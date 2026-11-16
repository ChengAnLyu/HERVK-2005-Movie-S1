"""Render the data-derived HERV-K 2005 review.

This is a scientific visualization of the supplied MRC arrays, not an imagegen
reconstruction. Geometry stays in the MRC physical coordinate system. Unknown
original scene settings are explicit parameters, never claimed as recovered.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import time

ROOT = Path(os.environ.get('HERVK_PROJECT_DIR', str(Path(__file__).resolve().parent)))
VTK_PRIVATE = Path(os.environ.get('HERVK_VTK_PATH', str(ROOT / 'vtk_egl')))
if VTK_PRIVATE.exists():
    sys.path.insert(0, str(VTK_PRIVATE))

import numpy as np
from scipy.sparse import coo_matrix
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from PIL import Image, ImageDraw, ImageFont

INPUT = next((p for p in (ROOT/'input', ROOT/'inputs') if p.is_dir()), ROOT)
CACHE = ROOT / 'render_cache'
OUT = ROOT / 'rendered'
PREVIEW = OUT / 'previews'
CACHE.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
PREVIEW.mkdir(parents=True, exist_ok=True)

WIDTH, HEIGHT, FPS = 1312, 844, 25
NX, NY, NZ = 660, 656, 185
CAPSID_COUNT = 4
SCAN_FRAMES = NZ-1
SCAN_END = 3*SCAN_FRAMES
OVERVIEW_END = SCAN_END+25
ZOOM_END = OVERVIEW_END+50
FIRST_END = ZOOM_END+200
TRANSITION_END = FIRST_END+50
TOTAL_FRAMES = TRANSITION_END+200
TOMOGRAM_CANDIDATES = (
    '_n2n_unet-medium_HervK_mature_2005_bin6_even(1).mrc',
    '_n2n_unet-medium_HervK_mature_2005_bin6_even.mrc',
)
# Fifty iterations were selected for a smooth review appearance. They are not
# recovered from the original session, which was not supplied.
SMOOTH_ITERATIONS = int(os.environ.get('HERVK_REVIEW_SMOOTHING', '50'))
SMOOTH_FACTOR = 0.3
LOW, HIGH = -0.55, 0.45
YELLOW, PALE_YELLOW, PURPLE = '#FFDB0F', '#FFFF7F', '#8B9CFF'
# The screenshot identifies cap3 as sushi-roll-like. Original movie zooms
# first to the multiple assemblies in cap4, then to the polygonal cap1.
SUSHI_CAPSIDS = (3,)
MULTIPLE_CAPSID = 4
FIRST_CAPSID, SECOND_CAPSID = 4, 1
FIRST_POSE = ((276, 232, 74), 230, -6, 5)
SECOND_POSE = ((249, 436, 93), 230, -14, 10)
# Brighter neutral surface appearance, selected against the supplied figure.
# The category RGB palette stays exact; lighting still conveys 3D shape.
SURFACE_AMBIENT, SURFACE_DIFFUSE, SURFACE_SPECULAR = 0.45, 0.60, 0.05
KEY_INTENSITY = 0.82


def capsid_color(number):
    return PURPLE if number == MULTIPLE_CAPSID else PALE_YELLOW if number in SUSHI_CAPSIDS else YELLOW


def resolve_tomogram():
    for name in TOMOGRAM_CANDIDATES:
        path = INPUT/name
        if path.is_file():
            return path
    raise FileNotFoundError('Missing denoised 2005 tomogram: '+', '.join(TOMOGRAM_CANDIDATES))


def ease(t):
    t = float(np.clip(t, 0, 1))
    return t*t*(3-2*t)


def hold_rock(frame):
    """Two seconds still, four seconds gently rocking, two seconds still."""
    if not 50 <= frame < 150:
        return 0.0
    wave = np.sin(2*np.pi*(frame-50)/99)
    return 10*wave*abs(wave)


def read_mrc(path):
    with open(path, 'rb') as f:
        h = f.read(1024)
    nx, ny, nz, mode = struct.unpack_from('<4i', h)
    grid = struct.unpack_from('<3i', h, 28)
    cell = struct.unpack_from('<3f', h, 40)
    axes = struct.unpack_from('<3i', h, 64)
    start = struct.unpack_from('<3i', h, 16)
    origin = struct.unpack_from('<3f', h, 196)
    nsymbt = struct.unpack_from('<i', h, 92)[0]
    assert axes == (1, 2, 3), axes
    assert start == (0, 0, 0), start
    dtype = {0: np.int8, 1: '<i2', 2: '<f4', 6: '<u2'}[mode]
    data = np.memmap(path, mode='r', dtype=dtype, offset=1024 + nsymbt,
                     shape=(nz, ny, nx))
    meta = dict(name=Path(path).name, shape_xyz=(nx, ny, nz), mode=mode,
                spacing=tuple(c / n for c, n in zip(cell, grid)),
                origin=origin, axes=axes, start=start)
    return data, meta


def smooth_equal_neighbors(poly):
    """The documented ChimeraX rule: move toward mean adjacent vertices."""
    xyz = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64)
    faces = vtk_to_numpy(poly.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    i = np.concatenate((faces[:, 0], faces[:, 1], faces[:, 2]))
    j = np.concatenate((faces[:, 1], faces[:, 2], faces[:, 0]))
    row, col = np.r_[i, j], np.r_[j, i]
    adj = coo_matrix((np.ones(len(row), dtype=np.float32), (row, col)),
                     shape=(len(xyz), len(xyz))).tocsr()
    adj.data[:] = 1
    degree = np.asarray(adj.sum(axis=1)).ravel()
    good = degree > 0
    for _ in range(SMOOTH_ITERATIONS):
        avg = adj @ xyz
        avg[good] /= degree[good, None]
        xyz[good] += SMOOTH_FACTOR * (avg[good] - xyz[good])
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(xyz.astype(np.float32), deep=True))
    poly.SetPoints(points)
    norm = vtk.vtkPolyDataNormals()
    norm.SetInputData(poly)
    norm.ComputePointNormalsOn()
    norm.ComputeCellNormalsOff()
    norm.SplittingOff()
    norm.ConsistencyOn()
    norm.AutoOrientNormalsOff()
    norm.Update()
    final = vtk.vtkPolyData()
    final.DeepCopy(norm.GetOutput())
    return final


def mesh_for_capsid(number):
    cache_path = CACHE / f'capsid_{number}_binary_union_s{SMOOTH_ITERATIONS}_f03.vtp'
    if cache_path.exists():
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(cache_path))
        reader.Update()
        poly = vtk.vtkPolyData()
        poly.DeepCopy(reader.GetOutput())
        return poly

    values, meta = read_mrc(INPUT / f'2005_capsid_{number}.mrc')
    # All positive segmentation IDs were checked; they are retained.
    # These integer region IDs are not physical density amplitudes.
    support = values > 0
    zz, yy, xx = np.where(support)
    # Take source padding where available. Do NOT invent zero voxels beyond
    # the original volume boundary: cropped capsids remain open at that edge.
    lo = np.maximum([int(xx.min())-1, int(yy.min())-1, int(zz.min())-1], 0)
    hi = np.minimum([int(xx.max())+1, int(yy.max())+1, int(zz.max())+1],
                    np.asarray(meta['shape_xyz'])-1)
    cropped = np.ascontiguousarray(support[lo[2]:hi[2]+1,
                                           lo[1]:hi[1]+1,
                                           lo[0]:hi[0]+1], dtype=np.uint8)
    del support, xx, yy, zz, values
    img = vtk.vtkImageData()
    img.SetDimensions(cropped.shape[2], cropped.shape[1], cropped.shape[0])
    img.SetSpacing(meta['spacing'])
    img.SetOrigin(np.asarray(meta['origin']) + lo * np.asarray(meta['spacing']))
    img.GetPointData().SetScalars(numpy_to_vtk(cropped.ravel(), deep=True,
                                              array_type=vtk.VTK_UNSIGNED_CHAR))
    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputData(img)
    contour.SetValue(0, 0.5)
    contour.ComputeScalarsOff()
    contour.ComputeNormalsOff()
    contour.Update()
    poly = vtk.vtkPolyData()
    poly.DeepCopy(contour.GetOutput())
    poly = smooth_equal_neighbors(poly)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(cache_path))
    writer.SetInputData(poly)
    writer.Write()
    print(f'Mesh {number}: {poly.GetNumberOfPoints():,} vertices, '
          f'{poly.GetNumberOfPolys():,} triangles', flush=True)
    return poly


def rgb(hexcolor):
    return tuple(int(hexcolor[i:i+2], 16)/255 for i in (1, 3, 5))


def font(size, bold=False):
    stem = 'DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf'
    return ImageFont.truetype('/usr/share/fonts/truetype/dejavu/' + stem, size)


class Scene:
    def __init__(self, width=WIDTH, height=HEIGHT):
        self.width, self.height = width, height
        self.scale = width/WIDTH
        self.map, self.meta = read_mrc(resolve_tomogram())
        assert tuple(self.meta['shape_xyz']) == (NX, NY, NZ)
        self.spacing = np.array(self.meta['spacing'])
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0, 0, 0)
        self.renderer.SetTwoSidedLighting(True)
        self.renderer.SetUseDepthPeeling(False)
        self.renderer.SetUseFXAA(True)
        self.window = vtk.vtkOSOpenGLRenderWindow()
        self.window.SetOffScreenRendering(1)
        self.window.SetSize(width, height)
        self.window.SetMultiSamples(0)
        self.window.AddRenderer(self.renderer)

        self.actors = []
        for i in range(1, CAPSID_COUNT+1):
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(mesh_for_capsid(i))
            mapper.ScalarVisibilityOff()
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            prop = actor.GetProperty()
            prop.SetColor(rgb(capsid_color(i)))
            prop.SetAmbient(SURFACE_AMBIENT)
            prop.SetDiffuse(SURFACE_DIFFUSE)
            prop.SetSpecular(SURFACE_SPECULAR)
            prop.SetSpecularPower(28)
            prop.SetInterpolationToPhong()
            prop.BackfaceCullingOff()
            actor.VisibilityOff()
            self.renderer.AddActor(actor)
            self.actors.append(actor)

        self.image = vtk.vtkImageData()
        self.image.SetDimensions(NX, NY, 1)
        self.image.SetSpacing(*self.spacing)
        self.image.SetOrigin(0, 0, 0)
        self.gray = np.empty((NY, NX), dtype=np.uint8)
        self.gray_vtk = numpy_to_vtk(self.gray.ravel(), deep=False,
                                      array_type=vtk.VTK_UNSIGNED_CHAR)
        self.image.GetPointData().SetScalars(self.gray_vtk)
        self.slice_actor = vtk.vtkImageActor()
        self.slice_actor.SetInputData(self.image)
        self.slice_actor.GetProperty().SetColorWindow(255)
        self.slice_actor.GetProperty().SetColorLevel(127.5)
        self.slice_actor.GetProperty().SetInterpolationTypeToLinear()
        self.slice_actor.GetProperty().SetOpacity(1)
        self.renderer.AddActor(self.slice_actor)

        self.camera = self.renderer.GetActiveCamera()
        self.camera.ParallelProjectionOn()
        self.camera.SetViewUp(0, 1, 0)
        self.target0 = np.array([(NX-1)/2, (NY-1)/2, (NZ-1)/2]) * self.spacing
        self.initial_scale = 344.4 * self.spacing[1]
        self.first_target = np.array(FIRST_POSE[0])*self.spacing
        self.second_target = np.array(SECOND_POSE[0])*self.spacing
        self.first_scale = FIRST_POSE[1]/2*self.spacing[1]
        self.second_scale = SECOND_POSE[1]/2*self.spacing[1]
        self.distance = 18000.0
        self.set_camera(self.target0, self.initial_scale, 0, 0)
        self.renderer.SetAmbient(1, 1, 1)
        self.renderer.RemoveAllLights()
        kit = vtk.vtkLightKit()
        kit.SetKeyLightIntensity(KEY_INTENSITY)
        kit.SetKeyLightWarmth(0.5)
        kit.SetFillLightWarmth(0.5)
        kit.SetBackLightWarmth(0.5)
        kit.SetHeadLightWarmth(0.5)
        kit.SetKeyToFillRatio(2.5)
        kit.SetKeyToBackRatio(4)
        kit.SetKeyToHeadRatio(5)
        kit.MaintainLuminanceOn()
        kit.AddLightsToRenderer(self.renderer)
        self.lightkit = kit

        self.capture = vtk.vtkWindowToImageFilter()
        self.capture.SetInput(self.window)
        self.capture.SetInputBufferTypeToRGB()
        self.capture.ReadFrontBufferOff()

    def set_camera(self, target, parallel_scale, yaw, pitch):
        y, p = np.radians([yaw, pitch])
        direction = np.array([np.sin(y)*np.cos(p), np.sin(p), np.cos(y)*np.cos(p)])
        self.camera.SetFocalPoint(target)
        self.camera.SetPosition(target + self.distance*direction)
        self.camera.SetViewUp(0, 1, 0)
        self.camera.SetParallelScale(parallel_scale)
        self.camera.SetClippingRange(100, 45000)

    def set_slice(self, z):
        self.gray[:] = np.clip((self.map[z] - LOW) * (255/(HIGH-LOW)), 0, 255).astype(np.uint8)
        self.gray_vtk.Modified()
        self.image.Modified()
        self.slice_actor.SetPosition(0, 0, z*self.spacing[2])

    def frame(self, i):
        if not 0 <= i < TOTAL_FRAMES:
            raise ValueError(f'Frame {i} outside 0..{TOTAL_FRAMES-1}')
        # Only visibility/opacity and the camera change; no model is moved.
        for actor in self.actors:
            actor.SetVisibility(False)
            actor.GetProperty().SetOpacity(1.0)
        if i < SCAN_END:
            z = SCAN_FRAMES-i if i < SCAN_FRAMES else (i-SCAN_FRAMES if i < 2*SCAN_FRAMES else SCAN_END-i)
            self.slice_actor.VisibilityOn()
            self.set_slice(z)
            self.set_camera(self.target0, self.initial_scale, 0, 0)
            for a in self.actors:
                a.SetVisibility(i >= 2*SCAN_FRAMES)
        else:
            z = None
            self.slice_actor.VisibilityOff()
            if i < OVERVIEW_END:
                for actor in self.actors:
                    actor.VisibilityOn()
                self.set_camera(self.target0, self.initial_scale, 0, 0)
            elif i < ZOOM_END:
                u = (i-OVERVIEW_END)/49
                e = ease(u)
                for number, actor in enumerate(self.actors, 1):
                    opacity = 1.0 if number == FIRST_CAPSID else 1-ease(u/0.4)
                    actor.SetVisibility(opacity > 0)
                    actor.GetProperty().SetOpacity(opacity)
                target = (1-e)*self.target0 + e*self.first_target
                scale = self.initial_scale*np.exp(np.log(self.first_scale/self.initial_scale)*e)
                self.set_camera(target, scale, FIRST_POSE[2]*e, FIRST_POSE[3]*e)
            elif i < FIRST_END:
                self.actors[FIRST_CAPSID-1].VisibilityOn()
                self.set_camera(self.first_target, self.first_scale,
                                FIRST_POSE[2]+hold_rock(i-ZOOM_END), FIRST_POSE[3])
            elif i < TRANSITION_END:
                u = (i-FIRST_END)/49
                e = ease(u)
                for number, opacity in ((FIRST_CAPSID, 1-ease(u/0.4)),
                                        (SECOND_CAPSID, ease((u-0.6)/0.4))):
                    actor = self.actors[number-1]
                    actor.SetVisibility(opacity > 0)
                    actor.GetProperty().SetOpacity(opacity)
                target = (1-e)*self.first_target + e*self.second_target
                scale = self.first_scale*np.exp(np.log(self.second_scale/self.first_scale)*e)
                self.set_camera(target, scale,
                                FIRST_POSE[2]+(SECOND_POSE[2]-FIRST_POSE[2])*e,
                                FIRST_POSE[3]+(SECOND_POSE[3]-FIRST_POSE[3])*e)
            else:
                self.actors[SECOND_CAPSID-1].VisibilityOn()
                self.set_camera(self.second_target, self.second_scale,
                                SECOND_POSE[2]+hold_rock(i-TRANSITION_END), SECOND_POSE[3])

        self.window.Render()
        self.capture.Modified()
        self.capture.Update()
        raw = vtk_to_numpy(self.capture.GetOutput().GetPointData().GetScalars())
        frame = Image.fromarray(np.flipud(raw.reshape(self.height, self.width, 3))).copy()
        self.add_annotations(frame, i, z)
        return frame

    def add_annotations(self, im, i, z):
        s=self.scale
        # Keep the scientific image area unobscured during the slice scans.
        d=ImageDraw.Draw(im)
        d.rectangle((0,0,round(238*s),self.height), fill='black')
        def xy(x,y):return (round(x*s),round(y*s))
        def txt(x,y,text,size=19,color='#D5D9E0',bold=False):
            d.multiline_text(xy(x,y),text,font=font(max(10,round(size*s)),bold),fill=color,spacing=round(7*s))
        txt(24,30,'HERV-K',24,'white',True)
        txt(24,67,'Tomogram 2005',17,'#A9B1BD')
        if i < 2*SCAN_FRAMES:
            txt(24,125,'Denoised\ntomographic slices',19)
        elif i < SCAN_END:
            txt(24,125,'Capsid\nsegmentations',19)
        else:
            txt(24,125,'Capsids in 3D',19)
        if 2*SCAN_FRAMES <= i < OVERVIEW_END:
            d.rounded_rectangle((*xy(24,222),*xy(40,238)),radius=max(2,round(3*s)),fill=YELLOW)
            txt(24,252,'Polygonal capsids',17,'#EDF0F5')
            d.rounded_rectangle((*xy(24,305),*xy(40,321)),radius=max(2,round(3*s)),fill=PALE_YELLOW)
            txt(24,335,'Sushi-roll-like\ncapsids',17,'#EDF0F5')
            d.rounded_rectangle((*xy(24,414),*xy(40,430)),radius=max(2,round(3*s)),fill=PURPLE)
            txt(24,444,'Multiple non-nested\ncapsid assemblies',17,'#EDF0F5')
        elif i >= OVERVIEW_END:
            first = i < FIRST_END+25
            color = capsid_color(FIRST_CAPSID if first else SECOND_CAPSID)
            label = 'Multiple non-nested\ncapsid assemblies' if first else 'Polygonal capsids'
            d.rounded_rectangle((*xy(24,222),*xy(44,242)),radius=max(2,round(3*s)),fill=color)
            txt(24,262,label,19,'#EDF0F5',True)
        if z is not None:
            txt(24,HEIGHT-66,f'Slice {z+1:03d} / {NZ}',15,'#98A2B2')
            # Scale bar on the black right margin, calibrated to the actual
            # orthographic world scale (500 A = 50 nm).
            bar=500/(2*self.camera.GetParallelScale()/self.height)
            x=round(1130*s);y=round((HEIGHT-65)*s)
            d.line((x,y,x+bar,y),fill='white',width=max(2,round(4*s)))
            txt(1128,HEIGHT-50,'50 nm',14,'#CED4DF')


def write_metadata():
    data=dict(
        title='Native mature HERV-K capsid morphology — tomogram 2005', revision=1,
        source_tomogram=resolve_tomogram().name,
        source_capsids=[f'2005_capsid_{i}.mrc' for i in range(1,CAPSID_COUNT+1)],
        reference_video='HERVK_video_2005(1).mp4',
        reference_video_frames=712, reference_video_duration_seconds=28.48,
        grid_xyz=[NX,NY,NZ], resolution=[WIDTH,HEIGHT], fps=FPS,
        frames=TOTAL_FRAMES, duration_seconds=TOTAL_FRAMES/FPS,
        source_slices={'pass1':[184,1], 'pass2':[0,183], 'pass3':[184,1]},
        frame_range_convention='start inclusive, end exclusive',
        stage_frames={'grayscale':[0,2*SCAN_FRAMES],
                      'capsid_reveal':[2*SCAN_FRAMES,SCAN_END],
                      'capsid_overview':[SCAN_END,OVERVIEW_END],
                      'zoom_to_multiple':[OVERVIEW_END,ZOOM_END],
                      'multiple_non_nested_assemblies':[ZOOM_END,FIRST_END],
                      'transition_to_polygonal':[FIRST_END,TRANSITION_END],
                      'polygonal_capsid':[TRANSITION_END,TOTAL_FRAMES]},
        isolated_hold_seconds=8,
        isolated_hold_motion='2 s still; 4 s gentle +/-10 degree yaw rock; 2 s still',
        closeup_order=['2005_capsid_4.mrc','2005_capsid_1.mrc'],
        closeup_poses={'2005_capsid_4.mrc':FIRST_POSE,'2005_capsid_1.mrc':SECOND_POSE},
        pose_units='target XYZ voxel indices; vertical span in Y voxels; yaw/pitch degrees',
        closeup_background='black; tomogram hidden; only selected capsid volume shown',
        membrane_surfaces=False,
        coords='Original MRC physical coordinates; no individual transformations',
        colors={f'2005_capsid_{i}.mrc':capsid_color(i) for i in range(1,CAPSID_COUNT+1)},
        terminology={'2005_capsid_1.mrc':'Polygonal capsids',
                     '2005_capsid_2.mrc':'Polygonal capsids',
                     '2005_capsid_3.mrc':'Sushi-roll-like capsids',
                     '2005_capsid_4.mrc':'Multiple non-nested capsid assemblies'},
        classification_provenance={'2005_capsid_3.mrc':'User explicitly classified screenshot image(3).png',
                                   'other_capsids':'Original 2005 video and supplied Figure 1A'},
        palette_reference='Figure category key, consistent with accepted 010 revision 3',
        image_levels=[LOW,HIGH],
        surface_appearance={'ambient':SURFACE_AMBIENT,'diffuse':SURFACE_DIFFUSE,
                            'specular':SURFACE_SPECULAR,'specular_power':28,
                            'key_intensity':KEY_INTENSITY,'light_colors':'neutral',
                            'key_to_fill_ratio':2.5,'key_to_back_ratio':4,'key_to_head_ratio':5,
                            'matches':'Accepted 010 revision 3'},
        visualization_mesh={'construction':'positive-label support, binary 0.5 contour; no outside-volume padding',
                            'selection':'All positive segmentation labels in all four capsid maps',
                            'labels':{'2005_capsid_1.mrc':[13],'2005_capsid_2.mrc':[12,15],
                                      '2005_capsid_3.mrc':[5,6,7],'2005_capsid_4.mrc':[16,19,22,23,24]},
                            'source_voxels_modified':False,'labels_excluded':[],
                            'smoothing':'uniform adjacent-vertex averaging',
                            'iterations':SMOOTH_ITERATIONS,'factor':SMOOTH_FACTOR,
                            'original_smoothing_recovered':False,
                            'smoothing_matches_010_review':True,'generic_component_filtering':False},
        renderer='VTK 9.5.2 / OSMesa; data-derived review render, not executed in ChimeraX',
        limitations=['No CXS session or original movie script was supplied.',
                     'Original numeric camera and material settings cannot be recovered exactly.',
                     'Original scan timing and close-up target order retained; camera path rebuilt with longer holds.',
                     'Cap1 and cap2 retain the broad polygonal category; no unsupported subtype is assigned.'])
    (OUT/'HERVK_2005_render_settings.json').write_text(json.dumps(data,indent=2))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--preview',action='store_true')
    ap.add_argument('--movie',action='store_true')
    ap.add_argument('--width',type=int,default=WIDTH)
    ap.add_argument('--height',type=int,default=HEIGHT)
    args=ap.parse_args()
    start=time.time()
    scene=Scene(args.width,args.height)
    write_metadata()
    if args.preview:
        frames=[100,500,552,600,627,702,877,952]
        thumbs=[]
        for i in frames:
            im=scene.frame(i)
            p=PREVIEW/f'preview_{i:04d}.png'
            im.save(p)
            print(f'Preview {i}: {p}',flush=True)
            thumb=im.resize((656,422),Image.Resampling.LANCZOS)
            thumbs.append((i,thumb))
        sheet=Image.new('RGB',(1312,4*454),(22,25,30))
        draw=ImageDraw.Draw(sheet)
        for k,(i,thumb) in enumerate(thumbs):
            x=(k%2)*656;y=(k//2)*454
            draw.text((x+12,y+7),f'{i/FPS:.2f} s',fill='white',font=font(17))
            sheet.paste(thumb,(x,y+32))
        sheet.save(PREVIEW/'review_storyboard.png')
    if args.movie:
        output=OUT/'HERVK_2005_denoised_capsids_review.mp4'
        cmd=['ffmpeg','-y','-v','error','-f','rawvideo','-pix_fmt','rgb24',
             '-s',f'{args.width}x{args.height}','-r',str(FPS),'-i','-',
             '-an','-c:v','libx264','-preset','medium','-crf','17',
             '-pix_fmt','yuv420p','-movflags','+faststart',
             '-metadata','title=Native mature HERV-K capsids — tomogram 2005',str(output)]
        process=subprocess.Popen(cmd,stdin=subprocess.PIPE)
        for i in range(TOTAL_FRAMES):
            im=scene.frame(i)
            process.stdin.write(im.tobytes())
            if i%50==0:
                print(f'Frame {i}/{TOTAL_FRAMES}, elapsed {time.time()-start:.1f}s',flush=True)
        process.stdin.close()
        rc=process.wait()
        if rc:raise RuntimeError(f'ffmpeg exit {rc}')
        print(f'Complete: {output} ({output.stat().st_size:,} bytes)',flush=True)
    scene.window.Finalize()
    print(f'Elapsed {time.time()-start:.1f}s',flush=True)


if __name__=='__main__':
    main()
