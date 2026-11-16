Movie S1: native ChimeraX companion
==================================

Movie S1 uses the supplied corrected tomogram and four aligned capsid label
maps. The three isolated close-ups appear in the requested order: the middle
nested example (capsid_1), the sushi-roll-like example (capsid_3), and the
multi-capsid assemblies (capsid_4). Each has an eight-second close-up against
black. The title "Movie S1" stays visible throughout the recording; the slice
caption is "Tomographic slices."

Status and provenance
---------------------
The companion was checked for Python syntax, input-file selection, frame counts,
scan indices, category/color assignments, persistent title visibility, isolated
model visibility, and camera calculations. Its commands follow the official
ChimeraX documentation/source listed below. It has NOT been executed in ChimeraX
in the preparation environment.

The accompanying rendered movie was produced separately from the same data
using VTK. It is not a native ChimeraX export. The native lighting is an intended
visual counterpart, not a verified pixel match to the VTK rendering.

The original .cxs session was not supplied. Camera matrices, lighting, and
smoothing are explicit recreation settings rather than recovered session state.
Running this script in ChimeraX creates a new native .cxs session. The supplied
TIFF is a single reference image; the corrected MRC supplies the slice data and
the four segmentation MRCs supply the capsid geometry.

Run
---
1. Put Movie_S1_native_setup.py beside these five files, or set DATA_DIRECTORY
   near the top of the Python script:

   HervK_mature_2005_1_bin6_corrected.mrc
   2005_capsid_1(1).mrc
   2005_capsid_2(1).mrc
   2005_capsid_3(1).mrc
   2005_capsid_4(1).mrc

   For each capsid, its filename without "(1)" is also accepted. The corrected
   tomogram may also be named HervK_mature_2005_1_bin6_corrected(1).mrc. Where
   both variants exist, the "(1)" version is preferred. The script also checks
   an inputs subfolder beside itself. It selects only the corrected tomogram
   filename listed here; there is no fallback to a previous tomogram product.

2. Start an empty ChimeraX session. Save existing work before starting a new
   session. The script stops if models are already open.

3. Enter this in the ChimeraX command line:

   open "/your/path/Movie_S1_native_setup.py"

   Do not run the script using a system Python interpreter.

4. By default, the script builds the scene, saves a native session, and records
   the movie. Outputs go in a native_output folder beside the script:

   Movie_S1_native.mp4
   Movie_S1_native_ready.cxs
   Movie_S1_native_record.cxc

   The saved session includes derived occupancy volumes and can be large.
   Keep the input MRC files accessible when restoring the native session.
   The five source MRC files are never edited.

To prepare without recording
---------------------------
Set RECORD_MOVIE = False, then open the Python file in ChimeraX. The scene and
native recording script will be prepared. Inspect the scene and open the
generated .cxc when ready. Its model IDs come from the scene built by the setup.

To repeat a recording, restore the generated .cxs in a fresh ChimeraX session
and open its corresponding .cxc. If outputs have moved, edit the movie output
path near the end of the .cxc. The recording script resets model visibility,
surface opacity, labels, and camera before recording. Its explicit camera
matrices do not require the setup Python script to remain loaded.

Exact labels and material colors
--------------------------------
Capsid_1, the middle example:
  Polygonal (nested capsids)
  Gold #dfbc0d, RGB (223,188,13)
  FIRST isolated close-up

Capsid_3:
  Polygonal (Sushi-roll-like)
  Bright yellow #f9f900, RGB (249,249,0)
  SECOND isolated close-up

Capsid_4:
  Multi-capsid assemblies
  Purple #7076de, RGB (112,118,222)
  THIRD isolated close-up

Capsid_2, the upper context example:
  Polygonal
  Gold #dfbc0d, RGB (223,188,13)
  Retained in the overview; no nested/sushi subtype is assigned.

The upper context example stays gold, consistent with the reference image.
The overview key therefore uses the broad term "Polygonal" for gold, plus
"Polygonal (Sushi-roll-like)" for bright yellow and "Multi-capsid assemblies"
for purple. The full term "Polygonal (nested capsids)" appears when the middle
capsid_1 is isolated. This avoids assigning a nested subtype to the other gold
capsid while keeping the requested colors.

Label support and surface construction
--------------------------------------
The capsid inputs are integer REGION LABEL maps, not density maps or binary
masks. The new capsid headers, label counts, and extents match the previously
audited segmentations. Every positive label is retained in every file:

   occupancy = (source labels > 0)

The script surfaces that derived binary 0/1 occupancy at 0.5. Numeric region IDs
are not interpreted as density values. The audit found no additional region
requiring exclusion. Internal sheets, disconnected assemblies, and small
positive-label pieces remain included.

There is no largest-component selection, surface dust removal, erosion, padding,
resampling, or separate membrane surface. Boundary capping is disabled with
capFaces false, avoiding artificial closure at the limits of a source array.
RETAINED_LABELS is empty, so no label-specific exclusion is applied.

Alignment, smoothing, and corrected-map contrast
------------------------------------------------
All five MRCs have matching 660 x 656 x 185 grids, zero origins/starts, and XYZ
axes. Supplied spacing is approximately 9.324 Angstrom per voxel. The script
checks interpreted grid size, origin, spacing, cell angles, and rotation before
constructing the surfaces, and preserves those properties in the derived maps.

No capsid receives an independent translation, rotation, axis flip, or fit.
Isolation changes visibility only. Camera movement changes the view of the
fixed aligned scene; it does not move a capsid relative to the tomogram or to
another capsid.

surfaceSmoothing remains enabled with 50 iterations and factor 0.3. These are
the accepted recreation settings, not parameters recovered from an original
.cxs session. Smoothing changes displayed surface geometry without changing any
input label voxel.

The corrected tomogram has a different intensity distribution from the earlier
input. Its fixed display window is -1.35 (black) to 0.30 (white), applied to every
slice. There is no frame-by-frame contrast adjustment. These display settings
do not change the source MRC data.

The image plane is hidden at 22.08 seconds, when the capsids-only section begins.
The 3D overview and close-ups have a fully black background. Surfaces are opaque
except for the short fade-out/fade-in transitions between isolated examples.

Native lighting
---------------
The three material RGB colors remain exactly as specified above. Neutral white
lighting controls display brightness. Both setup and the generated recording
apply the following as one command:

   lighting default color white fillColor white ambientColor white
     ambientIntensity 1.0 intensity 0.5 fillIntensity 0.25
     depthCue false shadows false multiShadow 0 moveWithCamera true

Directional shading is retained to show surface shape. Strong directional
shadows, ambient shadows, and depth darkening are disabled. Lights follow the
camera and silhouettes remain off. These are the accepted native lighting
settings used for the previous recreation.

The corresponding VTK rendering uses ambient/diffuse/specular coefficients of
0.45/0.60/0.05 and key intensity 0.82. ChimeraX lighting intensities and VTK
material coefficients control different rendering operations. The native
settings have not been calibrated for identical pixels in a ChimeraX runtime.
LIGHTING_COMMAND is exposed near the top of the script for native adjustment.

Timing, persistent title, and camera
-----------------------------------
Resolution: 1312 x 844
Playback:   25 frames per second
Frames:     1327
Duration:   53.08 seconds

Frames 0-183:      0.00-7.36 s    Tomographic slices, Z indices 184 down to 1
Frames 184-367:    7.36-14.72 s   Tomographic slices, Z indices 0 up to 183
Frames 368-551:    14.72-22.08 s  Capsid reveal, Z indices 184 down to 1
Frames 552-576:    22.08-23.08 s  All capsids on black, no image plane
Frames 577-626:    23.08-25.08 s  Isolate and zoom to middle capsid_1
Frames 627-826:    25.08-33.08 s  Polygonal (nested capsids), alone
Frames 827-876:    33.08-35.08 s  Transition from capsid_1 to capsid_3
Frames 877-1076:   35.08-43.08 s  Polygonal (Sushi-roll-like), alone
Frames 1077-1126:  43.08-45.08 s  Transition from capsid_3 to capsid_4
Frames 1127-1326:  45.08-53.08 s  Multi-capsid assemblies, alone

Frame numbers and Z indices are zero-based. Time intervals have an exclusive
end. Each scan index is written explicitly in the command script and followed
by one frame wait, with no implicit index interpolation.

"Movie S1" remains visible during every stage, including transition gaps.
The title is rebuilt whenever other labels change, before the next recorded
frame. This adds no wait frames. "Tomographic slices" is shown while the image
plane is present and is removed when the capsids-only section begins.

Each eight-second close-up holds still for two seconds, gently rocks by up to
10 degrees either side of its base yaw over four seconds, then holds still for
two seconds. The rock returns to its starting orientation. During each two-
second camera transition, the previous example fades out over the first 0.8
seconds and the next fades in over the final 0.8 seconds, separated by a 0.4-
second gap. Model coordinates remain fixed throughout.

The camera stays orthographic with clipping disabled. Targets are X,Y,Z indices
in the original full tomogram. Vertical span is measured in Y voxels; yaw and
pitch are degrees:

Overview: (329.5,327.5,92), span 688.8, yaw 0,   pitch 0
Capsid_1: (249,436,93),     span 230,   yaw -14, pitch 10
Capsid_3: (85,377,81),      span 210,   yaw 0,   pitch 0
Capsid_4: (276,232,74),     span 230,   yaw -6,  pitch 5

These poses are reconstructed around the requested subjects and their audited
extents. They are not original camera matrices recovered from a session.

The initial setup wait occurs before recording. The recording contains exactly
1327 frame waits: 552 + 25 + 50 + 200 + 50 + 200 + 50 + 200. movie status is
printed before encoding so the native captured count can be checked. Open the
.cxc as a file; typing its commands individually can introduce extra frames.

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
