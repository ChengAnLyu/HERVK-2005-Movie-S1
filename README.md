# Movie S1 — Mature HERV-K Capsid Morphologies

Scripts and rendering metadata for **Movie S1**, visualizing mature HERV-K capsid morphologies in tomogram 2005. The updated sequence scans the corrected tomogram, reveals aligned capsid segmentations, and presents three isolated close-ups: **Polygonal (nested capsids)**, **Polygonal (Sushi-roll-like)**, and **Multi-capsid assemblies**. The title **Movie S1** remains visible throughout.

**Default movie settings:** 1312 × 844 pixels · 25 frames per second · 1,327 frames · 53.08 seconds · no audio.

## Files

| File | Purpose |
| --- | --- |
| [Movie_S1_render_vtk.py](Movie_S1_render_vtk.py) | VTK renderer for the movie, preview images, and storyboard. |
| [Movie_S1_native_setup.py](Movie_S1_native_setup.py) | Native UCSF ChimeraX scene setup, session saving, and recording commands. |
| [Movie_S1_render_settings.json](Movie_S1_render_settings.json) | Reference record of timing, colors, camera poses, and rendering parameters. |
| [Movie_S1_notes.txt](Movie_S1_notes.txt) | Detailed visualization notes, provenance, and suggested caption. |
| [README.txt](README.txt) | Detailed native ChimeraX instructions and technical notes. |

The five supplied files are preserved as provided. This Markdown README summarizes the updated package and provides links to the full documentation. Input MRC data, rendered movies, PNG previews, and native session files are not included in this repository.

## Required data

Supply the corrected tomogram and four capsid label maps:

```text
HervK_mature_2005_1_bin6_corrected.mrc
2005_capsid_1(1).mrc
2005_capsid_2(1).mrc
2005_capsid_3(1).mrc
2005_capsid_4(1).mrc
```

The scripts also accept the corrected-tomogram name with `(1)` before `.mrc` and capsid filenames without `(1)`. The VTK renderer prefers the plain corrected-tomogram filename; the native companion prefers its `(1)` variant when both are available in the same searched directory. Both prefer the `(1)` capsid filenames. The earlier denoised tomogram is superseded by the corrected input in this version.

The four capsid maps are integer region-label maps. The scripts retain all positive labels, form binary occupancy, and surface it at 0.5. The default surface smoothing uses 50 iterations with factor 0.3. The fixed tomogram display window is **−1.35 to 0.30**. Source MRC data are not edited.

## Render with VTK

The supplied renderer requires NumPy, SciPy, Pillow, VTK 9.5+ with `vtkOSOpenGLRenderWindow`, an OSMesa runtime, and `ffmpeg` on `PATH` with `libx264` encoding support. Its `font()` function uses the Linux DejaVu Sans font paths under `/usr/share/fonts/truetype/dejavu/`; adjust that function if your system uses a different font layout.

Put all five MRCs together beside the script, or together in an `input/` or `inputs/` subfolder. By default, the renderer chooses the first existing `input/` or `inputs/` directory, otherwise the script directory; it does not combine files across those locations.

```bash
python Movie_S1_render_vtk.py --preview
python Movie_S1_render_vtk.py --movie
```

Default outputs:

- `rendered/Movie_S1.mp4`
- `rendered/Movie_S1_render_settings.json`
- With `--preview`: eight preview PNGs and `rendered/previews/review_storyboard.png`

Mesh caches are stored in `render_cache/`. The environment variables `HERVK_PROJECT_DIR`, `HERVK_INPUT_DIR`, and `HERVK_OUTPUT_DIR` select project, input, and output folders. `HERVK_VTK_PATH` can select a private VTK installation, and `HERVK_REVIEW_SMOOTHING` controls the smoothing-iteration count.

**The JSON file is reference metadata, not a runtime configuration file.** The VTK renderer writes a settings record from its Python parameters. Neither script loads the supplied JSON to configure rendering.

## Native ChimeraX workflow

Place the five MRCs beside `Movie_S1_native_setup.py` or in its `inputs/` subfolder, or set `DATA_DIRECTORY` near the top of that script. In an empty UCSF ChimeraX session, run:

```text
open "/your/path/Movie_S1_native_setup.py"
```

Run this script through ChimeraX's Python environment. It checks for an empty session, builds the scene, saves a session by default, and generates recording commands. With the default `RECORD_MOVIE = True`, it also records the movie.

Outputs go to `native_output/`:

- `Movie_S1_native.mp4`
- `Movie_S1_native_ready.cxs`
- `Movie_S1_native_record.cxc`

Set `RECORD_MOVIE = False` to prepare the scene and commands before recording. `SAVE_NATIVE_SESSION` separately controls session saving. See [README.txt](README.txt) for repeat recordings, alignment checks, lighting, and camera details.

## Updated sequence and colors

Tomographic scans occupy the first 14.72 seconds, followed by the capsid reveal through 22.08 seconds and an all-capsid overview on black. The three isolated eight-second highlights follow in this order:

| Capsid map | Classification | Material color | Isolated interval |
| --- | --- | --- | --- |
| `2005_capsid_1(1).mrc` | Polygonal (nested capsids) | Gold `#dfbc0d` | 25.08–33.08 s |
| `2005_capsid_3(1).mrc` | Polygonal (Sushi-roll-like) | Yellow `#f9f900` | 35.08–43.08 s |
| `2005_capsid_4(1).mrc` | Multi-capsid assemblies | Purple `#7076de` | 45.08–53.08 s |

Capsid 2 remains gold **Polygonal** context in the overview. Each isolated highlight holds still for two seconds, gently rocks in yaw for four seconds, then holds still for two seconds. Camera movement, visibility, and opacity produce the close-ups while the capsids retain their original alignment. The full frame schedule and camera poses are recorded in the settings JSON and supplied notes.

## Provenance and validation

The supplied notes describe a separately produced VTK movie and a companion native ChimeraX workflow. The original CXS session was not supplied; camera and lighting settings are explicit recreations. Matching pixels between VTK and ChimeraX have not been established.

For this repository update, both Python files passed static syntax checks, the JSON parsed successfully, and the documented timeline was checked against the settings. The MRC inputs were not part of this upload, and neither renderer was executed here. Prior data-audit and movie-production statements in the supplied notes remain the author's provenance record.
