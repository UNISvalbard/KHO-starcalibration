# Star calibration of all-sky images

All-sky auroral imagers use fish-eye lenses that capture the whole sky in one frame. **Star calibration** finds the zenith pixel, the lens scale `k`, and the rotation angle `z` (offset from geodetic to geomagnetic north) — geometric camera calibration using stars as known targets (distinct from intensity calibration).

## Basic procedure

1. Pick an image with clearly visible stars
2. From date, time, and location, get each star’s altitude and azimuth (e.g. Astropy, Stellarium)
3. Record each star’s pixel position `(row, col)` in the image
4. Fit the camera parameters (zenith location, rotation, lens model)

More stars give a better fit; four can already be “pretty good” for many uses. See e.g. [A detailed study of auroral fragments](http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-388546) by Joshua Dreyer.

For a more automatic pipeline, see Juha’s [Widefield star calibrator](https://github.com/jvierine/widefield-star-calibrator).

---

## Python (recommended)

`Python/star_calibration.py` estimates zenith `(row, col)`, rotation about the optical axis, and lens parameters from known star positions.

### Setup

```bash
pip install -r requirements.txt
```

### Quick start (worked example)

Edit the `CONFIG` block at the top of `Python/star_calibration.py` (`MODEL`, `IMAGE_PATH`, `STAR_TABLE_PATH`), then from the repo root:

```bash
cd Python
python star_calibration.py
```

This loads an image + star table (TSV), fits, prints parameters / RMS, and plots measured vs predicted stars with north/east markers.

Example star tables live in `selected_stars_WISC/` (columns include `altitude_deg`, `azimuth_deg`, `star_row_px_1based`, `star_col_px_1based`). Sample images are under `Test_images/`.

### Library usage

```python
from star_calibration import calibrate_camera, predict_star_positions

result = calibrate_camera(
    img.shape[:2], star_alt, star_az, star_row, star_col,
    model="equisolid",          # or "equidistant" | "equisolid_ext"
)
# result: zenith_row/col, rot_angle(_deg), params, rms_error, ...

pred_row, pred_col = predict_star_positions(result, star_alt, star_az)
```

### Lens models

| Model | Mapping (`z` = zenith angle) | Parameters |
| ----- | ---------------------------- | ---------- |
| `equidistant` | `d = k · z` | `k` (px/deg) |
| `equisolid` | `d = k · sin(z/2)` | `k` (px) |
| `equisolid_ext` | `d = k · sin(a·z) + p` | `k`, `a`, `p` |

Defined in `Python/lens_models.py` (same models as the MATLAB `starcalibration*.m` scripts).

---

## MATLAB (legacy)

Original scripts and map-projection examples are in `Matlab/`. Prefer the Python tools above for new work; see [`Matlab/README.md`](Matlab/README.md) for details.

---

## Other notes

- Example calibrations / images for instruments at [Kjell Henriksen Observatory (KHO)](http://kho.unis.no)
- Narrow-field sCMOS at KHO: use Nick’s [calibration scripts](https://github.com/n-brindley/all-sky-calibration) (tracking mount → calibrate per pointing)
