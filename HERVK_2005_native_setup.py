"""Recreate the HERVK 2005 update in UCSF ChimeraX.

Open this file from ChimeraX, not from a system Python interpreter.
Requires the denoised tomogram and all four uploaded capsid MRC files.
The original .cxs session was not supplied. Colors, smoothing, lighting, and
camera paths below are explicit recreation settings, not recovered settings.
This script has been checked against official ChimeraX documentation/source,
but has not been executed in ChimeraX in the preparation environment.
"""

from pathlib import Path
import math
import numpy as np

from chimerax.core.commands import run, quote_path_if_necessary
from chimerax.core.errors import UserError
from chimerax.map import Volume, volume_from_grid_data
from chimerax.map_data import ArrayGridData


# ------------------------------ USER SETTINGS ------------------------------
# Leave empty to look beside this script, in its inputs subfolder, or in the
# accompanying workspace's upload folder. Otherwise set an absolute directory.
DATA_DIRECTORY = ""
OUTPUT_DIRECTORY = ""  # Empty: a native_output folder beside this script.
RECORD_MOVIE = True     # False: prepare the scene and recording script only.
SAVE_NATIVE_SESSION = True

WIDTH, HEIGHT = 1312, 844
FPS = 25
SUPERSAMPLE = 1
CONTRAST_LOW, CONTRAST_HIGH = -0.55, 0.45  # Fixed for every tomogram plane.

# Chosen after visual review of the reference/recreation, not recovered from
# the missing original session. Surface smoothing remains enabled.
SMOOTHING_ITERATIONS = 50
SMOOTHING_FACTOR = 0.3

# Brighter neutral lighting for the figure palette. Ambient light is raised
# modestly from the full preset, with weaker directional lights and no shadow
# or depth darkening. This is a native visual counterpart to the rendered
# review; ChimeraX and VTK coefficients do not imply a verified pixel match.
LIGHTING_COMMAND = (
    "lighting default color white fillColor white ambientColor white "
    "ambientIntensity 1.0 intensity 0.5 fillIntensity 0.25 "
    "depthCue false shadows false multiShadow 0 moveWithCamera true")

POLYGONAL_COLOR = "#FFDB0F"
SUSHI_COLOR = "#FFFF7F"
MULTIPLE_COLOR = "#8B9CFF"
POLYGONAL_NAME = "Polygonal capsids"
SUSHI_NAME = "Sushi-roll-like capsids"
MULTIPLE_NAME = "Multiple non-nested capsid assemblies"
SHOW_CATEGORY_LEGEND = True  # Appears from the first surface reveal onward.

# Assignments follow comparison with the 2005 reference movie and the user's
# marked screenshot. No category assignments carry over from another dataset.
CAPSID_CATEGORIES = {1: "polygonal", 2: "polygonal", 3: "sushi", 4: "multiple"}
CATEGORY_STYLES = {
    "polygonal": (POLYGONAL_COLOR, POLYGONAL_NAME),
    "sushi": (SUSHI_COLOR, SUSHI_NAME),
    "multiple": (MULTIPLE_COLOR, MULTIPLE_NAME),
}

# These are REGION LABEL maps. The dataset audit found no membrane-like region
# to exclude: all positive labels remain, including capsid_2's internal sheet
# and every separate assembly in capsid_4. Empty selection dictionary = retain
# all positive labels. The displayed derived BINARY contour remains at 0.5.
RETAINED_LABELS = {}

# Recreated camera poses. Coordinates are X,Y,Z grid indices in the original
# full tomogram. Span is the vertical field of view, measured in Y voxels.
# All camera movement preserves the relative positions of every capsid.
OVERVIEW = ((329.5, 327.5, 92.0), 688.8, 0.0, 0.0)
CLOSEUPS = (
    {"index": 4, "pose": ((276.0, 232.0, 74.0), 230.0, -6.0, 5.0)},
    {"index": 1, "pose": ((249.0, 436.0, 93.0), 230.0, -14.0, 10.0)},
)
# Each scan is an explicit sequence of zero-based Z indices. This preserves
# the reference timing even if a source scan repeats or skips an index.
SCAN_Z_INDICES = (tuple(range(184, 0, -1)), tuple(range(184)), tuple(range(184, 0, -1)))
OVERVIEW_FRAMES = 25
ZOOM_FRAMES = 50
TRANSITION_FRAMES = 50
HOLD_FRAMES = 200
STATIC_HOLD_FRAMES = 50
ROCK_DEGREES = 10.0
FADE_FRAMES = 20

TOMOGRAM_FILES = (
    "_n2n_unet-medium_HervK_mature_2005_bin6_even(1).mrc",
    "_n2n_unet-medium_HervK_mature_2005_bin6_even.mrc",
)
CAPSID_FILES = [f"2005_capsid_{i}.mrc" for i in range(1, 5)]
EXPECTED_GRID_SIZE = (660, 656, 185)


def quote(value):
    return quote_path_if_necessary(str(value))


def model_spec(model):
    return "#" + model.id_string


def find_inputs(script_dir):
    candidates = ([Path(DATA_DIRECTORY).expanduser()] if DATA_DIRECTORY else
                  [script_dir, script_dir / "inputs",
                   script_dir.parent.parent / "upload"])
    for directory in candidates:
        if all((directory / name).is_file() for name in CAPSID_FILES):
            for name in TOMOGRAM_FILES:
                if (directory / name).is_file():
                    return directory.resolve(), (directory / name).resolve()
    locations = "\n".join(str(p) for p in candidates)
    raise UserError(
        "Could not find all five HERVK 2005 MRC files. Put them beside this "
        "script or set DATA_DIRECTORY near the top of the file. Searched:\n"
        + locations)


def open_one_volume(path):
    before = set(session.models.list(type=Volume))
    run(session, "open " + quote(path))
    opened = [m for m in session.models.list(type=Volume) if m not in before]
    if len(opened) != 1:
        raise UserError(f"Expected one volume from {path.name}; got {len(opened)}")
    opened[0].display = False
    return opened[0]


def check_alignment(reference, other):
    a, b = reference.data, other.data
    if tuple(a.size) != tuple(b.size):
        raise UserError(f"Grid size mismatch for {other.name}: {b.size} vs {a.size}")
    for field in ("origin", "step", "cell_angles", "rotation"):
        if not np.allclose(getattr(a, field), getattr(b, field), rtol=1e-6, atol=1e-5):
            raise UserError(f"Grid {field} mismatch for {other.name}; align the input maps first")


def make_capsid_surface(path, reference, color, category, index):
    source = open_one_volume(path)
    check_alignment(reference, source)
    values = source.full_matrix()
    if not np.issubdtype(values.dtype, np.integer):
        raise UserError(f"{path.name} must be an integer label map, got {values.dtype}")
    if values.min() < 0:
        raise UserError(f"{path.name} contains negative labels; check the label convention")

    # Select the source labels, then convert their union to explicit 0/1
    # occupancy. There is no component-size filtering, erosion, or padding.
    # All positive source labels are retained unless an explicit audited
    # selection is listed in RETAINED_LABELS for this dataset.
    source_positive = int(np.count_nonzero(values))
    occupancy = np.empty(values.shape, dtype=np.uint8)
    if index in RETAINED_LABELS:
        labels = tuple(RETAINED_LABELS[index])
        if not labels:
            raise UserError(f"The label selection for {path.name} is empty")
        occupancy[...] = np.isin(values, labels)
        selection_name = "labels " + ",".join(str(label) for label in labels)
    else:
        np.greater(values, 0, out=occupancy)
        selection_name = "all positive labels"
    kept_voxels = int(np.count_nonzero(occupancy))
    if kept_voxels == 0:
        raise UserError(f"{path.name} has no voxels matching {selection_name}")

    data = source.data
    grid = ArrayGridData(
        occupancy, origin=data.origin, step=data.step,
        cell_angles=data.cell_angles, rotation=data.rotation,
        name=f"{path.stem} - {category} ({selection_name})")
    volume = volume_from_grid_data(grid, session, style="surface", show_dialog=False)
    volume.position = source.scene_position
    run(session,
        f"volume {model_spec(volume)} region all style surface level 0.5 "
        f"step 1 limitVoxelCount false capFaces false color {color} "
        f"surfaceSmoothing true smoothingIterations {SMOOTHING_ITERATIONS} "
        f"smoothingFactor {SMOOTHING_FACTOR} calculateSurfaces true")

    # 0.5 is now the boundary of the derived BINARY occupancy, not a contour
    # applied to the arbitrary numeric magnitudes of the original label IDs.
    # capFaces false leaves cuts at the original array boundary open. This is
    # preserved for this dataset without adding new boundary faces.
    session.logger.info(
        f"{path.name}: {selection_name}; retained {kept_voxels:,} of "
        f"{source_positive:,} positive voxels; excluded "
        f"{source_positive - kept_voxels:,} voxels")
    session.models.close([source])
    return volume


def interpolate_pose(start, end, fraction):
    # Smooth acceleration and deceleration; geometric interpolation of zoom.
    t = fraction * fraction * (3.0 - 2.0 * fraction)
    target = tuple((1.0 - t) * a + t * b for a, b in zip(start[0], end[0]))
    span = math.exp((1.0 - t) * math.log(start[1]) + t * math.log(end[1]))
    yaw = (1.0 - t) * start[2] + t * end[2]
    pitch = (1.0 - t) * start[3] + t * end[3]
    return target, span, yaw, pitch


def camera_commands(reference, pose):
    # Camera local -Z is the viewing direction. Rotation affects only the
    # camera; no density map or segmentation is transformed independently.
    target_ijk, span_voxels, yaw_degrees, pitch_degrees = pose
    target = np.asarray(reference.data.ijk_to_xyz(target_ijk), dtype=float)
    yaw, pitch = math.radians(yaw_degrees), math.radians(pitch_degrees)
    cy, sy, cp, sp = math.cos(yaw), math.sin(yaw), math.cos(pitch), math.sin(pitch)
    ry = np.array(((cy, 0., sy), (0., 1., 0.), (-sy, 0., cy)))
    # Positive pitch means the camera moves toward +Y, matching the review.
    rx = np.array(((1., 0., 0.), (0., cp, sp), (0., -sp, cp)))
    rotation = ry @ rx
    extent = np.asarray(reference.data.size) * np.asarray(reference.data.step)
    distance = 3.0 * float(np.linalg.norm(extent))
    origin = target + rotation @ np.array((0., 0., distance))
    matrix = np.column_stack((rotation, origin)).reshape(-1)
    matrix_text = ",".join(f"{v:.12g}" for v in matrix)
    pixel_size = span_voxels * reference.data.step[1] / HEIGHT
    return f"view matrix camera {matrix_text}; zoom pixelSize {pixel_size:.12g}"


def category_for_capsid(index):
    return CATEGORY_STYLES[CAPSID_CATEGORIES[index]]


def legend_commands(stage, visible=True):
    # Labels are rebuilt at stage boundaries, so saved label IDs are not
    # required for replay. The setup requires an empty initial session.
    lines = ["2dlabels delete all"]
    if not SHOW_CATEGORY_LEGEND or stage is None:
        return lines
    visibility = "true" if visible else "false"

    def label(text, color, xpos, ypos, size=22):
        return (f"2dlabels text {quote(text)} color {color} bgColor black "
                f"margin 3 size {size} xpos {xpos} ypos {ypos} "
                f"visibility {visibility}")

    if stage == "overview":
        present = set(CAPSID_CATEGORIES.values())
        row = 0
        for category in ("polygonal", "sushi", "multiple"):
            if category not in present:
                continue
            color, text = CATEGORY_STYLES[category]
            ypos = 0.94 - row * 0.045
            lines.extend([label("\u25a0", color, 0.025, ypos),
                          label(text, "white", 0.057, ypos)])
            row += 1
    elif isinstance(stage, int):
        color, text = category_for_capsid(stage)
        lines.append(label(text, color, 0.025, 0.94, 26))
    else:
        raise ValueError(f"Unknown legend stage: {stage}")
    return lines


def isolate_commands(capsid_volumes, index):
    # Change visibility only. All model coordinates and parent positions
    # remain fixed. Parent groups stay visible while their children change.
    return [f"{'show' if number == index else 'hide'} {model_spec(volume)} models"
            for number, volume in enumerate(capsid_volumes, start=1)]


def held_pose(base, frame):
    rock_frames = HOLD_FRAMES - 2 * STATIC_HOLD_FRAMES
    if frame < STATIC_HOLD_FRAMES or frame >= HOLD_FRAMES - STATIC_HOLD_FRAMES:
        return base
    u = (frame - STATIC_HOLD_FRAMES) / (rock_frames - 1)
    sine = math.sin(2.0 * math.pi * u)
    # Signed sin-squared: a slow +/-10-degree rock with zero velocity at
    # each extremum and at the joins to the static portions of the hold.
    yaw_offset = ROCK_DEGREES * sine * abs(sine)
    return base[0], base[1], base[2] + yaw_offset, base[3]


def recording_commands(tomogram, capsids, capsid_volumes, movie_path):
    ts, cs = model_spec(tomogram), model_spec(capsids)
    first_index, second_index = (closeup["index"] for closeup in CLOSEUPS)
    first_pose, second_pose = (closeup["pose"] for closeup in CLOSEUPS)
    first = model_spec(capsid_volumes[first_index - 1])
    second = model_spec(capsid_volumes[second_index - 1])
    first_name, second_name = (category_for_capsid(i)[1] for i in (first_index, second_index))
    first_z = SCAN_Z_INDICES[0][0]
    lines = [
        "# Generated by HERVK_2005_native_setup.py. Run in the prepared scene.",
        "# All camera frames are explicit, so no original .cxs views are implied.",
        f"windowsize {WIDTH} {HEIGHT}",
        "set bgColor black", "clip off", "camera ortho", LIGHTING_COMMAND,
        "graphics silhouettes false",
        f"show {ts} models", f"show {cs} models",
        f"transparency {cs} 0 target s",
    ]
    # Restore every child after a prior recording, then hide their parent.
    lines.extend(f"show {model_spec(volume)} models" for volume in capsid_volumes)
    lines.extend([
        f"hide {cs} models",
        f"volume {ts} style image color white colorMode l8 step 1 planes z,{first_z}",
        f"volume {ts} level {CONTRAST_LOW},0 level {CONTRAST_HIGH},1 "
        "colormapExtendLeft true colormapExtendRight true showOutlineBox false",
        "2dlabels delete all",
        camera_commands(tomogram, OVERVIEW),
        "wait 1",  # Outside recording: flush setup and create display resources.
        f"movie record size {WIDTH},{HEIGHT} format png supersample {SUPERSAMPLE}",
    ])

    current_frame = 0
    for pass_index, z_indices in enumerate(SCAN_Z_INDICES):
        end_frame = current_frame + len(z_indices) - 1
        lines.append(f"# Frames {current_frame}-{end_frame}: scan pass {pass_index + 1}.")
        if pass_index == 2:
            lines.append(f"show {cs} models")
            lines.extend(legend_commands("overview"))
        for z_index in z_indices:
            lines.append(f"volume {ts} planes z,{z_index}; wait 1")
        current_frame += len(z_indices)

    lines.extend([
        f"# Frames {current_frame}-{current_frame + OVERVIEW_FRAMES - 1}: capsids on black.",
        f"hide {ts} models", "set bgColor black",
        f"transparency {cs} 0 target s", f"wait {OVERVIEW_FRAMES}",
    ])
    current_frame += OVERVIEW_FRAMES
    lines.append(f"# Frames {current_frame}-{current_frame + ZOOM_FRAMES - 1}: zoom to capsid_{first_index}.")
    lines.extend(isolate_commands(capsid_volumes, first_index))
    lines.extend(legend_commands(first_index))
    for frame in range(1, ZOOM_FRAMES + 1):
        pose = interpolate_pose(OVERVIEW, first_pose, frame / ZOOM_FRAMES)
        lines.append(camera_commands(tomogram, pose) + "; wait 1")
    current_frame += ZOOM_FRAMES

    lines.append(f"# Frames {current_frame}-{current_frame + HOLD_FRAMES - 1}: {first_name} alone.")
    for frame in range(HOLD_FRAMES):
        lines.append(camera_commands(tomogram, held_pose(first_pose, frame)) + "; wait 1")
    current_frame += HOLD_FRAMES

    lines.append(f"# Frames {current_frame}-{current_frame + TRANSITION_FRAMES - 1}: move between fixed positions.")
    lines.extend(legend_commands(None))
    for frame in range(TRANSITION_FRAMES):
        if frame < FADE_FRAMES:
            percent = 100.0 * (frame + 1) / FADE_FRAMES
            lines.append(f"transparency {first} {percent:.8g} target s")
        elif frame == FADE_FRAMES:
            lines.append(f"hide {first} models")
        if frame >= TRANSITION_FRAMES - FADE_FRAMES:
            if frame == TRANSITION_FRAMES - FADE_FRAMES:
                lines.append(f"show {second} models")
                lines.extend(legend_commands(second_index))
            percent = 100.0 * (TRANSITION_FRAMES - frame - 1) / FADE_FRAMES
            lines.append(f"transparency {second} {percent:.8g} target s")
        pose = interpolate_pose(first_pose, second_pose,
                                (frame + 1) / TRANSITION_FRAMES)
        lines.append(camera_commands(tomogram, pose) + "; wait 1")
    current_frame += TRANSITION_FRAMES

    lines.append(f"# Frames {current_frame}-{current_frame + HOLD_FRAMES - 1}: {second_name} alone.")
    lines.extend(isolate_commands(capsid_volumes, second_index))
    lines.append(f"transparency {second} 0 target s")
    for frame in range(HOLD_FRAMES):
        lines.append(camera_commands(tomogram, held_pose(second_pose, frame)) + "; wait 1")
    current_frame += HOLD_FRAMES
    lines.extend([
        "movie stop", "movie status",
        f"movie encode {quote(movie_path)} framerate {FPS} quality high wait true",
        f"# Expected: {current_frame} frames, {current_frame / FPS:.2f} s at {FPS} fps.",
        "",
    ])
    return lines


def main():
    if session.models.list():
        raise UserError(
            "Run this setup in an empty ChimeraX session. Save your current "
            "work and start a new session first; this script does not close existing work.")
    if CONTRAST_HIGH <= CONTRAST_LOW:
        raise UserError("CONTRAST_HIGH must exceed CONTRAST_LOW")
    if HOLD_FRAMES - 2 * STATIC_HOLD_FRAMES < 2:
        raise UserError("The hold must include at least two frames for the gentle rock")
    if not 0 < 2 * FADE_FRAMES <= TRANSITION_FRAMES:
        raise UserError("The two fades must fit within TRANSITION_FRAMES")
    expected_indices = set(range(1, len(CAPSID_FILES) + 1))
    if set(CAPSID_CATEGORIES) != expected_indices:
        raise UserError("Each capsid map must have exactly one category assignment")
    if not set(CAPSID_CATEGORIES.values()).issubset(CATEGORY_STYLES):
        raise UserError("A capsid category has no defined name/color")
    if len(CLOSEUPS) != 2 or any(c["index"] not in expected_indices for c in CLOSEUPS):
        raise UserError("Configure two close-ups using valid capsid map indices")
    if len(SCAN_Z_INDICES) != 3 or any(not scan for scan in SCAN_Z_INDICES):
        raise UserError("Configure the three nonempty scan-index sequences")
    if any(z < 0 or z >= EXPECTED_GRID_SIZE[2] for scan in SCAN_Z_INDICES for z in scan):
        raise UserError("A scan Z index is outside the expected tomogram grid")
    script_dir = Path(__file__).resolve().parent
    input_dir, tomogram_path = find_inputs(script_dir)
    output_dir = (Path(OUTPUT_DIRECTORY).expanduser() if OUTPUT_DIRECTORY else
                  script_dir / "native_output").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    tomogram = open_one_volume(tomogram_path)
    if tuple(tomogram.data.size) != EXPECTED_GRID_SIZE:
        raise UserError(
            f"This recording is configured for {EXPECTED_GRID_SIZE}; "
            f"the tomogram has grid size {tomogram.data.size}")

    capsid_volumes = []
    for index, filename in enumerate(CAPSID_FILES, start=1):
        color, category = category_for_capsid(index)
        capsid_volumes.append(make_capsid_surface(
            input_dir / filename, tomogram, color, category, index))
    groups = []
    for category in ("polygonal", "sushi", "multiple"):
        volumes = [volume for index, volume in enumerate(capsid_volumes, start=1)
                   if CAPSID_CATEGORIES[index] == category]
        if volumes:
            groups.append(session.models.add_group(volumes, name=CATEGORY_STYLES[category][1]))
    capsids = session.models.add_group(groups, name="HERVK 2005 capsids")

    run(session, f"windowsize {WIDTH} {HEIGHT}")
    run(session, "set bgColor black; clip off; camera ortho; graphics silhouettes false")
    run(session, LIGHTING_COMMAND)
    run(session, f"transparency {model_spec(capsids)} 0 target s")
    run(session, f"hide {model_spec(capsids)} models; show {model_spec(tomogram)} models")
    run(session,
        f"volume {model_spec(tomogram)} style image color white colorMode l8 "
        f"step 1 planes z,{SCAN_Z_INDICES[0][0]} level {CONTRAST_LOW},0 level {CONTRAST_HIGH},1 "
        "colormapExtendLeft true colormapExtendRight true showOutlineBox false")
    run(session, camera_commands(tomogram, OVERVIEW))

    for command in legend_commands("overview", visible=False):
        run(session, command)
    run(session, "wait 1")

    movie_path = output_dir / "HERVK_2005_native_updated.mp4"
    record_path = output_dir / "HERVK_2005_native_record.cxc"
    commands = recording_commands(tomogram, capsids, capsid_volumes, movie_path)
    record_path.write_text("\n".join(commands), encoding="utf-8")
    if SAVE_NATIVE_SESSION:
        run(session, "save " + quote(output_dir / "HERVK_2005_native_ready.cxs"))
    session.logger.info("Native recording script: " + str(record_path))
    frame_count = sum(map(len, SCAN_Z_INDICES)) + OVERVIEW_FRAMES + ZOOM_FRAMES + 2 * HOLD_FRAMES + TRANSITION_FRAMES
    session.logger.info(
        "HERVK 2005: fixed contrast -0.55 to 0.45; all positive region labels retained; "
        "original grid alignment retained; no membrane surface; "
        f"surface smoothing {SMOOTHING_ITERATIONS} iterations at {SMOOTHING_FACTOR}; "
        "capsids 1/2 polygonal yellow, capsid 3 sushi-roll-like pale yellow, "
        "capsid 4 multiple non-nested purple; "
        "brighter neutral lighting; black close-up background; "
        "eight-second individual holds of capsid 4 then capsid 1; "
        f"{frame_count} frames / {frame_count / FPS:.2f} seconds.")
    if RECORD_MOVIE:
        # The actual recording runs as a native command script, with explicit
        # waits; a Python tight loop is not used to capture graphics frames.
        run(session, "open " + quote(record_path))
    else:
        session.logger.info("To record, enter: open " + quote(record_path))


main()
