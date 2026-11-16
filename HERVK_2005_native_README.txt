HERVK 2005: native ChimeraX companion
===================================

This companion recreates the updated 2005 movie from the supplied denoised
tomogram and four capsid label maps. It uses the accepted figure palette,
surface smoothing, brighter neutral lighting, black 3D background, and longer
individual close-ups. The marked screenshot identifies capsid_3 as sushi-roll-
like. The two close-ups retain the original movie's subject sequence: capsid_4
(multiple non-nested assemblies), then capsid_1 (polygonal).

Status and provenance
---------------------
The script was checked for Python syntax, input resolution, generated frame
counts, color/category assignments, scan indices, visibility, camera calculations,
and consistency with official ChimeraX documentation/source. It has NOT been
executed in ChimeraX in the preparation environment.

The accompanying rendered review movie was produced separately from the same
data using VTK. It is not a native ChimeraX export. The native lighting below
is an intended visual counterpart, not a verified pixel match to VTK.

The original .cxs session was not supplied. The camera subjects and scan pattern
were compared with HERVK_video_2005(1).mp4; camera matrices, materials, and
smoothing are explicit recreation settings rather than recovered session state.
This script creates a new native .cxs session when run in ChimeraX.

Run
---
1. Put HERVK_2005_native_setup.py beside these five files, or set DATA_DIRECTORY
   near the top of the script:

   _n2n_unet-medium_HervK_mature_2005_bin6_even(1).mrc
   2005_capsid_1.mrc
   2005_capsid_2.mrc
   2005_capsid_3.mrc
   2005_capsid_4.mrc

   The denoised filename without "(1)" is also accepted:
   _n2n_unet-medium_HervK_mature_2005_bin6_even.mrc

   If both denoised filenames exist in the same directory, the "(1)" version is
   selected first. The script also checks an inputs subfolder beside itself.

2. Start an empty ChimeraX session. Save existing work before starting a new
   session. The script stops if models are already open.

3. Enter this in the ChimeraX command line:

   open "/your/path/HERVK_2005_native_setup.py"

   Do not run the script using a system Python interpreter.

4. By default, the script builds the scene, saves a native session, and records
   the movie. Outputs go in a native_output folder beside the script:

   HERVK_2005_native_updated.mp4
   HERVK_2005_native_ready.cxs
   HERVK_2005_native_record.cxc

   The saved session includes derived occupancy volumes and can be large.
   Keep the original MRC inputs accessible when restoring the native session.
   The five source MRC files are never edited.

To prepare without recording
---------------------------
Set RECORD_MOVIE = False, then open the Python file in ChimeraX. The scene and
native recording script will be prepared. Inspect the scene and open the
generated .cxc when ready. Its model IDs come from the scene built by the setup.

To repeat the recording, restore the generated .cxs in a fresh ChimeraX session
and open its corresponding .cxc. If outputs have moved, edit the output movie
path near the end of the .cxc. The recording script resets model visibility,
surface opacity, labels, and camera before recording. Its explicit camera
matrices do not require the setup Python script to remain loaded.

Colors, categories, and close-up subjects
----------------------------------------
2005_capsid_1.mrc:
  Polygonal capsids, yellow #FFDB0F; SECOND isolated close-up

2005_capsid_2.mrc:
  Polygonal capsids, yellow #FFDB0F

2005_capsid_3.mrc:
  Sushi-roll-like capsids, pale yellow #FFFF7F
  This is the object identified in the user's marked screenshot.

2005_capsid_4.mrc:
  Multiple non-nested capsid assemblies, purple #8B9CFF; FIRST isolated close-up

The overview key explicitly shows all three colors and category names. The
individual close-up labels are "Multiple non-nested capsid assemblies" and
"Polygonal capsids," matching the displayed subjects. Capsid_3 stays pale yellow
in the overview; its screenshot classification does not replace either of the
original movie's two close-up subjects. No nested subtype is assigned to either
of the generic yellow polygonal examples.

Label support and surface construction
--------------------------------------
These inputs are integer REGION LABEL maps, not density maps or binary masks.
For every map this companion forms:

   occupancy = (source labels > 0)

It then surfaces the derived binary 0/1 occupancy at 0.5. Numeric region IDs
are not interpreted as density values. No positive region is excluded from
the 2005 dataset: the label audit did not identify an unwanted membrane-like
region requiring removal. Capsid_2's internal sheet is retained, all five
labels of capsid_4 remain as separate assemblies, and the small positive
labels in capsid_3 are preserved.

There is no largest-component selection, dust removal, erosion, padding,
resampling, or external membrane surface. Boundary capping is disabled with
capFaces false, avoiding artificial closure at the limits of a source array.
RETAINED_LABELS is empty in this script; no label-selection rule from the
earlier 010 dataset is applied here.

Alignment, smoothing, and image contrast
---------------------------------------
All five inputs have matching 660 x 656 x 185 grids, zero origins/starts, and
XYZ axes. Supplied spacing is approximately 9.324 Angstrom per voxel. The script
checks interpreted grid size, origin, spacing, cell angles, and rotation before
constructing the surfaces. Those grid properties are preserved in the derived
volumes.

No capsid receives an independent translation, rotation, axis flip, or fit.
Isolation changes visibility only. Camera movement changes the view of the
fixed aligned scene; it does not move a capsid relative to the tomogram or to
another capsid.

surfaceSmoothing remains enabled with 50 iterations and factor 0.3, matching
the accepted previous recreation. These are recreation values rather than
parameters recovered from an original .cxs session. Smoothing changes displayed
surface geometry without changing any input label voxel.

Denoised-image contrast stays fixed at -0.55 (black) to 0.45 (white) throughout
all three scans. The image plane is hidden at 22.08 seconds, when the capsids-
only section begins. The overview and close-ups then have a fully black
background. Surfaces are opaque except for the brief transition fades.

Brighter native display
-----------------------
The material RGB palette remains exactly #FFDB0F, #FFFF7F, and #8B9CFF. Native
display brightness is adjusted with neutral white lighting. Both setup and the
generated recording apply the following as one command:

   lighting default color white fillColor white ambientColor white
     ambientIntensity 1.0 intensity 0.5 fillIntensity 0.25
     depthCue false shadows false multiShadow 0 moveWithCamera true

This preserves directional shading for surface shape while disabling strong
directional shadows, ambient shadows, and depth darkening. Lights follow the
camera and silhouettes remain off. The values match the accepted 010 revision 3
native companion.

The corresponding VTK review uses ambient/diffuse/specular coefficients of
0.45/0.60/0.05 and key intensity 0.82. ChimeraX lighting intensities and VTK
material coefficients control different rendering operations. The native
settings are a brighter visual counterpart; they have not been calibrated for
identical pixels in a ChimeraX runtime. LIGHTING_COMMAND is exposed near the top
of the script for later native visual adjustment.

Timing and camera
-----------------
Resolution: 1312 x 844
Playback:   25 frames per second
Frames:     1077
Duration:   43.08 seconds

Frames 0-183:      0.00-7.36 s    Denoised image, Z indices 184 down to 1
Frames 184-367:    7.36-14.72 s   Denoised image, Z indices 0 up to 183
Frames 368-551:    14.72-22.08 s  Capsid reveal, Z indices 184 down to 1
Frames 552-576:    22.08-23.08 s  All capsids on black, no image plane
Frames 577-626:    23.08-25.08 s  Isolate and zoom to capsid_4
Frames 627-826:    25.08-33.08 s  Multiple non-nested assemblies alone
Frames 827-876:    33.08-35.08 s  Transition from capsid_4 to capsid_1
Frames 877-1076:   35.08-43.08 s  Polygonal example alone

Frame numbers and Z indices are zero-based. Time intervals have an exclusive
end. The first three scans preserve the reference index progression. Each scan
index is written explicitly in the command script and followed by one frame
wait, so there is no implicit index interpolation.

Each eight-second close-up holds still for two seconds, gently rocks by up to
10 degrees either side of its base yaw over four seconds, then holds still for
two seconds. The rock returns to its starting orientation. During the two-second
camera transition, capsid_4 fades out over the first 0.8 seconds and capsid_1
fades in over the final 0.8 seconds, with a 0.4-second gap between.

The native camera remains orthographic with clipping disabled. Targets are
X,Y,Z indices in the original full tomogram. Vertical span is measured in Y
voxels, and yaw/pitch are degrees:

Overview: (329.5,327.5,92), span 688.8, yaw 0, pitch 0
Capsid_4: (276,232,74),     span 230,   yaw -6, pitch 5
Capsid_1: (249,436,93),     span 230,   yaw -14, pitch 10

These poses were rebuilt around the reference movie's subjects and the audited
capsid extents. They are not original camera matrices recovered from a session.

The initial wait occurs before recording. During recording there are exactly
1077 frame waits: 552 + 25 + 50 + 200 + 50 + 200. movie status is printed before
encoding so the native captured count can be checked. Open the .cxc as a file;
typing its commands individually can introduce additional recorded frames.

Official references checked
---------------------------
https://rbvi.github.io/chimerax-recipes/planes/planes.html
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/volume.html
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/lighting.html
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/view.html
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/zoom.html
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/movie.html
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/2dlabels.html
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/transparency.html
https://github.com/RBVI/ChimeraX/blob/develop/src/bundles/map/src/volume.py
https://github.com/RBVI/ChimeraX/blob/develop/src/bundles/graphics/src/camera.py
https://github.com/RBVI/ChimeraX/blob/develop/src/bundles/core/src/models.py
