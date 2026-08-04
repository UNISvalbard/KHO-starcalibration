"""
star_calibration.py

Estimate all-sky camera calibration parameters for auroral/all-sky cameras
with fish-eye lenses, from a set of reference stars whose sky position
(altitude, azimuth) and pixel position (row, column) are both known.

There are always four "geometric" unknowns:
    - zenith_row, zenith_col : pixel location of the zenith
    - rot_angle              : rotation between north and the image's
                                column axis (radians)
plus one or more lens-model parameters (k, and for the extended equisolid
model also a and p) supplied by whichever LensModel is selected from
lens_models.py.

Library usage
-------------
    from star_calibration import calibrate_camera, predict_star_positions

    result = calibrate_camera(img.shape[:2], star_alt, star_az,
                               star_row, star_col, model="equisolid_ext")
    print(result["zenith_row"], result["zenith_col"], result["params"])

Run directly for a full worked example that loads an image and a star
table, calibrates, prints, and plots the result (edit the CONFIG section
in _run_example() below for your own files):

    python star_calibration.py
"""

from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from lens_models import MODELS, LensModel


# ============================================================================
# Core calibration engine
# ============================================================================

def _pack_params(zenith_row, zenith_col, extra_params, rot_angle):
    return np.concatenate(([zenith_row, zenith_col], extra_params, [rot_angle]))


def _unpack_params(x, n_extra):
    zenith_row, zenith_col = x[0], x[1]
    extra_params = x[2:2 + n_extra]
    rot_angle = x[2 + n_extra]
    return zenith_row, zenith_col, extra_params, rot_angle


def _predicted_positions(zenith_row, zenith_col, rot_angle, d, star_az_deg):
    theta = np.deg2rad(star_az_deg)
    row = zenith_row - d * np.cos(theta + rot_angle)
    col = zenith_col - d * np.sin(theta + rot_angle)
    return row, col


def _cost(x, model: LensModel, star_alt, star_az, star_row, star_col, n_extra):
    zenith_row, zenith_col, extra_params, rot_angle = _unpack_params(x, n_extra)
    z_deg = 90 - star_alt
    d = model.distance(z_deg, extra_params)
    row, col = _predicted_positions(zenith_row, zenith_col, rot_angle, d, star_az)
    return np.sum((star_row - row) ** 2 + (star_col - col) ** 2)


def calibrate_camera(imgsize, star_alt, star_az, star_row, star_col,
                      model="equisolid_ext", rot_angle0_deg=25.0):
    """
    Estimate all-sky camera calibration parameters by least-squares fitting
    of measured star pixel positions to their known sky positions.

    Parameters
    ----------
    imgsize : (n_rows, n_cols)
        Image size, e.g. img.shape[:2]
    star_alt, star_az : array-like
        Altitude and azimuth of the reference stars, in degrees
    star_row, star_col : array-like
        Pixel row / column of the reference stars in the image
    model : str or lens_models.LensModel
        One of 'equidistant', 'equisolid', 'equisolid_ext' (see
        lens_models.py), or a LensModel instance of your own
    rot_angle0_deg : float
        Initial guess for the rotation angle [degrees]. The default of 25
        degrees comes from the geodetic/geomagnetic north offset at KHO
        (Kjell Henriksen Observatory, Svalbard) in the original scripts --
        adjust it for your own site/camera.

    Returns
    -------
    dict with keys:
        zenith_row, zenith_col : float, pixel location of the zenith
        rot_angle               : float, radians
        rot_angle_deg           : float, degrees
        params                  : dict of the model's extra parameters,
                                   e.g. {'k': ...} or {'k': ..., 'a': ..., 'p': ...}
        rms_error                : float, pixels
        model                    : str, name of the lens model used
        success                  : bool, whether the optimiser converged
        fit                      : the raw scipy.optimize.OptimizeResult
    """
    star_alt = np.asarray(star_alt, dtype=float)
    star_az = np.asarray(star_az, dtype=float)
    star_row = np.asarray(star_row, dtype=float)
    star_col = np.asarray(star_col, dtype=float)

    lens_model = MODELS[model] if isinstance(model, str) else model
    n_extra = len(lens_model.param_names)

    # Rough initial guess: zenith at the image centre, k from a linear
    # regression through the origin of pixel distance vs. zenith angle,
    # rotation from the supplied default.
    zenith_row0 = imgsize[0] / 2
    zenith_col0 = imgsize[1] / 2

    star_zen_deg = 90 - star_alt
    star_distances = np.hypot(zenith_row0 - star_row, zenith_col0 - star_col)
    extra0 = lens_model.initial_guess(star_zen_deg, star_distances)
    rot0 = np.deg2rad(rot_angle0_deg)

    x0 = _pack_params(zenith_row0, zenith_col0, extra0, rot0)

    # Nelder-Mead is the direct equivalent of MATLAB's fminsearch.
    fit = minimize(
        _cost, x0,
        args=(lens_model, star_alt, star_az, star_row, star_col, n_extra),
        method="Nelder-Mead",
        options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 20000, "maxfev": 20000},
    )

    zenith_row, zenith_col, extra_params, rot_angle = _unpack_params(fit.x, n_extra)

    z_deg = 90 - star_alt
    d = lens_model.distance(z_deg, extra_params)
    row, col = _predicted_positions(zenith_row, zenith_col, rot_angle, d, star_az)
    rms_error = np.sqrt(np.mean((star_row - row) ** 2 + (star_col - col) ** 2))

    return {
        "zenith_row": zenith_row,
        "zenith_col": zenith_col,
        "rot_angle": rot_angle,
        "rot_angle_deg": np.rad2deg(rot_angle),
        "params": dict(zip(lens_model.param_names, extra_params)),
        "rms_error": rms_error,
        "model": lens_model.name,
        "success": bool(fit.success),
        "fit": fit,
    }


def predict_star_positions(result, star_alt, star_az):
    """
    Predicted pixel (row, col) positions of stars given a calibration
    result dict as returned by calibrate_camera().
    """
    lens_model = MODELS[result["model"]]
    extra_params = np.array([result["params"][name] for name in lens_model.param_names])
    star_alt = np.asarray(star_alt, dtype=float)
    star_az = np.asarray(star_az, dtype=float)
    z_deg = 90 - star_alt
    d = lens_model.distance(z_deg, extra_params)
    return _predicted_positions(result["zenith_row"], result["zenith_col"],
                                 result["rot_angle"], d, star_az)


# ============================================================================
# Worked example: load an image + star table, calibrate, print, plot
# (equivalent to WISC_test_sonyA7Scalibration3.m, generalised to any model)
# ============================================================================

def _run_example():
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    # ---- CONFIG: edit these for your own data --------------------------
    MODEL = "equisolid_ext"          # 'equidistant' | 'equisolid' | 'equisolid_ext'
    IMAGE_PATH = Path("Test_images") / "LYR-Sony-111220_200933.jpg"
    STAR_TABLE_PATH = Path("selected_stars_WISC") / "selected_stars_LYR_20201211.tsv"
    ROT_ANGLE0_DEG = 25.0             # initial guess (KHO geodetic/geomagnetic offset)
    # ----------------------------------------------------------------------

    img = mpimg.imread(IMAGE_PATH).astype(float)
    if img.max() > 1.0:
        img = img / 255.0

    stars = pd.read_csv(STAR_TABLE_PATH, sep="\t")
    star_names = stars["star"].to_numpy()
    star_az = stars["azimuth_deg"].to_numpy(dtype=float)
    star_alt = stars["altitude_deg"].to_numpy(dtype=float)
    star_col = stars["star_col_px_1based"].to_numpy(dtype=float)
    star_row = stars["star_row_px_1based"].to_numpy(dtype=float)

    result = calibrate_camera(img.shape[:2], star_alt, star_az, star_row, star_col,
                               model=MODEL, rot_angle0_deg=ROT_ANGLE0_DEG)

    print("-------------------------------------------------------")
    print(f"         Model = {result['model']}")
    print(f"        Zenith = ({result['zenith_row']:.1f}, {result['zenith_col']:.1f})")
    for name, val in result["params"].items():
        print(f"             {name} = {val:.3f}")
    print(f"Rotation angle = {result['rot_angle']:.3f} rad ({result['rot_angle_deg']:.1f} deg)")
    print(f"           RMS = {result['rms_error']:.2f} pixels")
    if not result["success"]:
        print("WARNING: optimiser did not report convergence -- inspect result['fit'].")

    new_row, new_col = predict_star_positions(result, star_alt, star_az)

    fig, ax = plt.subplots()
    ax.imshow(np.sqrt(img))
    ax.set_title(f"{IMAGE_PATH.name}\n{MODELS[MODEL].description}")

    ax.plot(result["zenith_col"], result["zenith_row"], "go")
    ax.plot(star_col, star_row, "ro", markersize=10)
    for name, c, r in zip(star_names, star_col, star_row):
        ax.text(c + 15, r, name, color="r")
    ax.plot(new_col, new_row, "g+")

    # North / east indicator lines, 45 degrees from zenith
    lens_model = MODELS[MODEL]
    extra_params = np.array([result["params"][n] for n in lens_model.param_names])
    d45 = lens_model.distance(np.array(45.0), extra_params)
    rot_angle = result["rot_angle"]

    north_row = result["zenith_row"] - d45 * np.cos(rot_angle)
    north_col = result["zenith_col"] - d45 * np.sin(rot_angle)
    ax.plot([result["zenith_col"], north_col], [result["zenith_row"], north_row], "b", linewidth=3)
    ax.text(north_col, north_row - 15, "NORTH", color="b")

    east_row = result["zenith_row"] - d45 * np.cos(np.pi / 2 + rot_angle)
    east_col = result["zenith_col"] - d45 * np.sin(np.pi / 2 + rot_angle)
    ax.plot([result["zenith_col"], east_col], [result["zenith_row"], east_row], "b", linewidth=3)
    ax.text(east_col, east_row + 15, "EAST", color="b")

    ax.axis("on")
    plt.show()


if __name__ == "__main__":
    _run_example()
