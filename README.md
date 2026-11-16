# Movie S1 — HERV-K Tomogram 2005

Scripts and rendering metadata for **Movie S1**, visualizing mature HERV-K capsid morphology in tomogram 2005. The animation scans denoised tomographic slices, reveals segmented capsids, and presents isolated close-ups of multiple non-nested capsid assemblies followed by a polygonal capsid.

**Movie settings:** 1312 × 844 pixels, 25 frames per second, 1,077 frames (43.08 seconds).

## Files

| File | Purpose |
| --- | --- |
| `render_2005_review_vtk.py` | VTK movie renderer and preview generation |
| `HERVK_2005_native_setup.py` | Native UCSF ChimeraX scene and recording setup |
| `HERVK_2005_render_settings.json` | Rendering parameters and provenance for reference; not loaded by the scripts |
| `HERVK_2005_native_README.txt` | Detailed ChimeraX instructions and visualization notes |

## Required data

Supply the denoised tomogram `_n2n_unet-medium_HervK_mature_2005_bin6_even(1).mrc` and four segmentation maps, `2005_capsid_1.mrc` through `2005_capsid_4.mrc`. These datasets are not included in this repository. The tomogram filename without `(1)` is also accepted; the `(1)` version is preferred when both are present.

Place all five files beside the scripts or in an `inputs/` subfolder.

## Render with VTK

Requires NumPy, SciPy, Pillow, VTK with OSMesa support, FFmpeg with `libx264`, and DejaVu Sans fonts. The script uses Linux font paths, which may need adjustment on other systems.

```bash
python render_2005_review_vtk.py --preview
python render_2005_review_vtk.py --movie
```

Outputs are written to `rendered/`.

## Native ChimeraX workflow

In an empty ChimeraX session, run:

```text
open "/path/to/HERVK_2005_native_setup.py"
```

Outputs are written to `native_output/`. See [the detailed instructions](HERVK_2005_native_README.txt) for input selection, camera settings, colors, and recording options.

The native companion was not executed in ChimeraX during original preparation. Camera and material settings were recreated; identical output between the ChimeraX and separate VTK workflows has not been verified.
