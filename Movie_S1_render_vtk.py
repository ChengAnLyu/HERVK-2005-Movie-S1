"""Render Movie S1 from the corrected HERV-K 2005 tomogram and capsid maps.

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

INPUT = Path(os.environ.get('HERVK_INPUT_DIR', str(next((p for p in (ROOT/'input', ROOT/'inputs') if p.is_dir()), ROOT))))
CACHE = ROOT / 'render_cache'
OUT = Path(os.environ.get('HERVK_OUTPUT_DIR', str(ROOT / 'rendered')))
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
TOMOGRAM_CANDIDATES = (
    'HervK_mature_2005_1_bin6_corrected.mrc',
    'HervK_mature_2005_1_bin6_corrected(1).mrc',
)
# Fifty iterations were selected for a smooth review appearance. They are not
# recovered from the original session, which was not supplied.
SMOOTH_ITERATIONS = int(os.environ.get('HERVK_REVIEW_SMOOTHING', '50'))
SMOOTH_FACTOR = 0.3
LOW, HIGH = -1.35, 0.30
POLYGONAL_COLOR, SUSHI_COLOR, MULTI_COLOR = '#dfbc0d', '#f9f900', '#7076de'
CLOSEUPS = [
    dict(index=1, name='Polygonal (nested capsids)', label='Polygonal\n(nested capsids)',
         color=POLYGONAL_COLOR, pose=((249,436,93),230,-14,10), key='nested_capsids'),
    dict(index=3, name='Polygonal (Sushi-roll-like)', label='Polygonal\n(Sushi-roll-like)',
         color=SUSHI_COLOR, pose=((85,377,81),210,0,0), key='sushi_roll_like'),
    dict(index=4, name='Multi-capsid assemblies', label='Multi-capsid\nassemblies',
         color=MULTI_COLOR, pose=((276,232,74),230,-6,5), key='multi_capsid_assemblies'),
]


def make_stages():
    stages = [dict(kind='overview', name='capsid_overview', start=SCAN_END, end=OVERVIEW_END)]
    cursor = OVERVIEW_END
    for j, closeup in enumerate(CLOSEUPS):
        stages.append(dict(kind='zoom' if j == 0 else 'transition',
                           name='zoom_to_'+closeup['key'], start=cursor, end=cursor+50,
                           target=j, previous=j-1))
        cursor += 50
        stages.append(dict(kind='hold', name=closeup['key'], start=cursor, end=cursor+200,
                           target=j))
        cursor += 200
    return stages, cursor


STAGES, TOTAL_FRAMES = make_stages()
# Brighter neutral surface appearance, selected against the supplied figure.
# The category RGB palette stays exact; lighting still conveys 3D shape.
SURFACE_AMBIENT, SURFACE_DIFFUSE, SURFACE_SPECULAR = 0.45, 0.60, 0.05
KEY_INTENSITY = 0.82


def capsid_color(number):
    return MULTI_COLOR if number == 4 else SUSHI_COLOR if number == 3 else POLYGONAL_COLOR


def resolve_capsid(number):
    for name in (f'2005_capsid_{number}(1).mrc', f'2005_capsid_{number}.mrc'):
        path = INPUT/name
        if path.is_file():
            return path
    raise FileNotFoundError(f'Missing capsid map {number} in {INPUT}')


def resolve_tomogram():
    for name in TOMOGRAM_CANDIDATES:
        path = INPUT/name
        if path.is_file():
            return path
    raise FileNotFoundError('Missing corrected 2005 tomogram: '+', '.join(TOMOGRAM_CANDIDATES))


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
    source_path = resolve_capsid(number)
    cache_path = CACHE / f'capsid_{number}_binary_union_s{SMOOTH_ITERATIONS}_f03.vtp'
    if cache_path.exists() and cache_path.stat().st_mtime >= source_path.stat().st_mtime:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(cache_path))
        reader.Update()
        poly = vtk.vtkPolyData()
        poly.DeepCopy(reader.GetOutput())
        return poly

    values, meta = read_mrc(source_path)
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
        self.closeup_targets = [np.array(c['pose'][0])*self.spacing for c in CLOSEUPS]
        self.closeup_scales = [c['pose'][1]/2*self.spacing[1] for c in CLOSEUPS]
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

    def interpolate_camera(self, previous, target, fraction):
        e = ease(fraction)
        if previous is None:
            a_target, a_scale, a_angles = self.target0, self.initial_scale, (0,0)
        else:
            a_target = self.closeup_targets[previous]
            a_scale = self.closeup_scales[previous]
            a_angles = CLOSEUPS[previous]['pose'][2:]
        b_target = self.closeup_targets[target]
        b_scale = self.closeup_scales[target]
        b_angles = CLOSEUPS[target]['pose'][2:]
        center = (1-e)*a_target + e*b_target
        scale = a_scale*np.exp(np.log(b_scale/a_scale)*e)
        self.set_camera(center, scale,
                        a_angles[0]+(b_angles[0]-a_angles[0])*e,
                        a_angles[1]+(b_angles[1]-a_angles[1])*e)

    def show_opacity(self, number, opacity):
        opacity = float(np.clip(opacity,0,1))
        actor = self.actors[number-1]
        actor.SetVisibility(opacity > 0)
        actor.GetProperty().SetOpacity(opacity)

    def frame(self, i):
        if not 0 <= i < TOTAL_FRAMES:
            raise ValueError(f'Frame {i} outside 0..{TOTAL_FRAMES-1}')
        for actor in self.actors:
            actor.SetVisibility(False)
            actor.GetProperty().SetOpacity(1.0)
        active = None
        if i < SCAN_END:
            z = SCAN_FRAMES-i if i < SCAN_FRAMES else (i-SCAN_FRAMES if i < 2*SCAN_FRAMES else SCAN_END-i)
            self.slice_actor.VisibilityOn()
            self.set_slice(z)
            self.set_camera(self.target0, self.initial_scale, 0, 0)
            for actor in self.actors:
                actor.SetVisibility(i >= 2*SCAN_FRAMES)
        else:
            z = None
            self.slice_actor.VisibilityOff()
            stage = next(stage for stage in STAGES if stage['start'] <= i < stage['end'])
            local_frame = i-stage['start']
            if stage['kind'] == 'overview':
                for actor in self.actors:
                    actor.VisibilityOn()
                self.set_camera(self.target0, self.initial_scale, 0, 0)
            else:
                active = stage['target']
                closeup = CLOSEUPS[active]
                if stage['kind'] == 'hold':
                    self.show_opacity(closeup['index'], 1)
                    self.set_camera(self.closeup_targets[active], self.closeup_scales[active],
                                    closeup['pose'][2]+hold_rock(local_frame),closeup['pose'][3])
                else:
                    u = local_frame/(stage['end']-stage['start']-1)
                    if stage['kind'] == 'zoom':
                        for number in range(1,CAPSID_COUNT+1):
                            self.show_opacity(number,1 if number == closeup['index'] else 1-ease(u/0.4))
                        self.interpolate_camera(None,active,u)
                    else:
                        previous = stage['previous']
                        self.show_opacity(CLOSEUPS[previous]['index'],1-ease(u/0.4))
                        self.show_opacity(closeup['index'],ease((u-0.6)/0.4))
                        self.interpolate_camera(previous,active,u)
                        if u < 0.5:
                            active = previous
        self.window.Render()
        self.capture.Modified()
        self.capture.Update()
        raw = vtk_to_numpy(self.capture.GetOutput().GetPointData().GetScalars())
        frame = Image.fromarray(np.flipud(raw.reshape(self.height,self.width,3))).copy()
        self.add_annotations(frame,i,z,active)
        return frame

    def add_annotations(self, im, i, z, active):
        s = self.scale
        d = ImageDraw.Draw(im)
        d.rectangle((0,0,round(238*s),self.height),fill='black')
        def xy(x,y): return (round(x*s),round(y*s))
        def txt(x,y,text,size=19,color='#D5D9E0',bold=False):
            d.multiline_text(xy(x,y),text,font=font(max(10,round(size*s)),bold),
                             fill=color,spacing=round(7*s))
        txt(24,30,'Movie S1',24,'white',True)
        txt(24,67,'Tomogram 2005',17,'#A9B1BD')
        if i < 2*SCAN_FRAMES:
            txt(24,125,'Tomographic slices',19)
        elif i < SCAN_END:
            txt(24,125,'Capsid\nsegmentations',19)
        else:
            txt(24,125,'Capsids in 3D',19)
        if 2*SCAN_FRAMES <= i < OVERVIEW_END:
            # The top gold context particle is generic Polygonal; the full
            # nested assignment is scoped to the middle particle's close-up.
            entries=[(222,POLYGONAL_COLOR,'Polygonal'),
                     (320,SUSHI_COLOR,'Polygonal\n(Sushi-roll-like)'),
                     (430,MULTI_COLOR,'Multi-capsid\nassemblies')]
            for y,color,label in entries:
                d.rounded_rectangle((*xy(24,y),*xy(40,y+16)),radius=max(2,round(3*s)),fill=color)
                txt(24,y+30,label,17,'#EDF0F5')
        elif active is not None:
            c = CLOSEUPS[active]
            d.rounded_rectangle((*xy(24,222),*xy(44,242)),radius=max(2,round(3*s)),fill=c['color'])
            txt(24,262,c['label'],19,'#EDF0F5',True)
        if z is not None:
            txt(24,HEIGHT-66,f'Slice {z+1:03d} / {NZ}',15,'#98A2B2')
            bar=500/(2*self.camera.GetParallelScale()/self.height)
            x=round(1130*s);y=round((HEIGHT-65)*s)
            d.line((x,y,x+bar,y),fill='white',width=max(2,round(4*s)))
            txt(1128,HEIGHT-50,'50 nm',14,'#CED4DF')


def write_metadata():
    stages={'tomographic_slices':[0,2*SCAN_FRAMES],
            'capsid_overlay':[2*SCAN_FRAMES,SCAN_END]}
    stages.update({stage['name']:[stage['start'],stage['end']] for stage in STAGES})
    classification={1:'Polygonal (nested capsids)',2:'Polygonal',
                    3:'Polygonal (Sushi-roll-like)',4:'Multi-capsid assemblies'}
    data=dict(
        title='Movie S1',dataset='HERV-K mature tomogram 2005',
        source_tomogram=resolve_tomogram().name,
        source_capsids=[resolve_capsid(i).name for i in range(1,CAPSID_COUNT+1)],
        visual_reference='2005 segmentation.tif (single rendered RGBA image, not a label volume)',
        grid_xyz=[NX,NY,NZ],resolution=[WIDTH,HEIGHT],fps=FPS,frames=TOTAL_FRAMES,
        duration_seconds=TOTAL_FRAMES/FPS,note='Tomographic slices',
        source_slices={'pass1':[184,1],'pass2':[0,183],'pass3':[184,1]},
        frame_range_convention='start inclusive, end exclusive',stage_frames=stages,
        stage_seconds={name:[a/FPS,b/FPS] for name,(a,b) in stages.items()},
        closeup_order=[c['index'] for c in CLOSEUPS],closeups=CLOSEUPS,
        isolated_hold_seconds=8,
        isolated_hold_motion='2 s still; 4 s gentle +/-10 degree yaw rock; 2 s still',
        pose_units='XYZ voxel indices; vertical span in Y voxels; yaw/pitch degrees',
        closeup_background='black; tomogram hidden; selected capsid volume only',
        membrane_surfaces=False,source_data_modified=False,
        coords='Original MRC physical coordinates; no individual transforms',
        colors={resolve_capsid(i).name:capsid_color(i) for i in range(1,CAPSID_COUNT+1)},
        terminology={resolve_capsid(i).name:classification[i] for i in range(1,CAPSID_COUNT+1)},
        top_particle_note='Cap2 remains gold Polygonal context, without a nested subtype assignment.',
        image_levels=[LOW,HIGH],
        contrast_reason='Fixed window selected for the newly supplied corrected tomogram intensity distribution.',
        surface_appearance={'ambient':SURFACE_AMBIENT,'diffuse':SURFACE_DIFFUSE,
                            'specular':SURFACE_SPECULAR,'specular_power':28,
                            'key_intensity':KEY_INTENSITY,'light_colors':'neutral',
                            'key_to_fill_ratio':2.5,'key_to_back_ratio':4,'key_to_head_ratio':5},
        visualization_mesh={'construction':'All positive region labels, binary 0.5 contour; no outside-volume padding',
                            'labels':{'cap1':[13],'cap2':[12,15],'cap3':[5,6,7],
                                      'cap4':[16,19,22,23,24]},
                            'labels_excluded':[],'generic_component_filtering':False,
                            'smoothing':'uniform adjacent-vertex averaging',
                            'iterations':SMOOTH_ITERATIONS,'factor':SMOOTH_FACTOR},
        renderer='VTK 9.5.2 / OSMesa; data-derived render, not executed in ChimeraX',
        limitations=['The original CXS session and numerical camera parameters were not supplied.',
                     'Camera path is reconstructed for three isolated highlights.',
                     'Native ChimeraX companion has not been executed in this environment.'])
    (OUT/'Movie_S1_render_settings.json').write_text(json.dumps(data,indent=2))


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
        frames=[100,500,552,627,702,877,1127,1202]
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
        output=OUT/'Movie_S1.mp4'
        cmd=['ffmpeg','-y','-v','error','-f','rawvideo','-pix_fmt','rgb24',
             '-s',f'{args.width}x{args.height}','-r',str(FPS),'-i','-',
             '-an','-c:v','libx264','-preset','medium','-crf','17',
             '-pix_fmt','yuv420p','-movflags','+faststart',
             '-metadata','title=Movie S1',str(output)]
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
