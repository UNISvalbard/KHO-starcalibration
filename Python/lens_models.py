"""
lens_models.py

Fish-eye lens projection models for auroral all-sky camera star calibration.

Each model describes how the radial pixel distance `d` of a star from the
image zenith relates to that star's zenith angle `z` (its angular distance
from straight up, in degrees: 0 = zenith, 90 = horizon), via one or more
calibration parameters.

    equidistant     d = k * z                  (k in pixels/degree)
    equisolid       d = k * sin(z / 2)          (k in pixels)
    equisolid_ext   d = k * sin(a * z) + p      (k, p in pixels; a dimensionless)

These correspond to starcalibration.m, starcalibration2.m and
starcalibration3.m in the original MATLAB code, respectively.

Each LensModel bundles:
    - distance(z_deg, params) -> d
        the projection formula itself
    - initial_guess(star_zen_deg, star_distances) -> params0
        a rough starting estimate for its extra parameters, used to seed
        the optimiser in star_calibration.py

Add a new lens model by writing a distance()/initial_guess() pair and
registering it in MODELS at the bottom of this file.
"""

from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np


@dataclass(frozen=True)
class LensModel:
    name: str
    param_names: Tuple[str, ...]
    distance: Callable[[np.ndarray, np.ndarray], np.ndarray]
    initial_guess: Callable[[np.ndarray, np.ndarray], np.ndarray]
    description: str = ""


# ----------------------------------------------------------------------
# 1) Equidistant projection:  d = k * z
# ----------------------------------------------------------------------

def _equidistant_distance(z_deg: np.ndarray, params: np.ndarray) -> np.ndarray:
    (k,) = params
    return k * z_deg


def _equidistant_initial_guess(star_zen_deg: np.ndarray, star_distances: np.ndarray) -> np.ndarray:
    k0 = np.mean(star_distances / star_zen_deg)
    return np.array([k0])


EQUIDISTANT = LensModel(
    name="equidistant",
    param_names=("k",),
    distance=_equidistant_distance,
    initial_guess=_equidistant_initial_guess,
    description="d = k * z  [k in pixels/degree]",
)


# ----------------------------------------------------------------------
# 2) Equisolid projection:  d = k * sin(z / 2)
# ----------------------------------------------------------------------

def _equisolid_distance(z_deg: np.ndarray, params: np.ndarray) -> np.ndarray:
    (k,) = params
    z_rad = np.deg2rad(z_deg)
    return k * np.sin(z_rad / 2)


def _equisolid_initial_guess(star_zen_deg: np.ndarray, star_distances: np.ndarray) -> np.ndarray:
    z_rad = np.deg2rad(star_zen_deg)
    k0 = np.mean(star_distances / np.sin(z_rad / 2))
    return np.array([k0])


EQUISOLID = LensModel(
    name="equisolid",
    param_names=("k",),
    distance=_equisolid_distance,
    initial_guess=_equisolid_initial_guess,
    description="d = k * sin(z / 2)  [k in pixels]",
)


# ----------------------------------------------------------------------
# 3) Generalised equisolid projection:  d = k * sin(a * z) + p
# ----------------------------------------------------------------------

def _equisolid_ext_distance(z_deg: np.ndarray, params: np.ndarray) -> np.ndarray:
    k, a, p = params
    z_rad = np.deg2rad(z_deg)
    return k * np.sin(a * z_rad) + p


def _equisolid_ext_initial_guess(star_zen_deg: np.ndarray, star_distances: np.ndarray) -> np.ndarray:
    z_rad = np.deg2rad(star_zen_deg)
    k0 = np.mean(star_distances / np.sin(z_rad / 2))
    return np.array([k0, 0.5, 0.0])


EQUISOLID_EXT = LensModel(
    name="equisolid_ext",
    param_names=("k", "a", "p"),
    distance=_equisolid_ext_distance,
    initial_guess=_equisolid_ext_initial_guess,
    description="d = k * sin(a * z) + p  [k, p in pixels; a dimensionless]",
)


MODELS = {
    "equidistant": EQUIDISTANT,
    "equisolid": EQUISOLID,
    "equisolid_ext": EQUISOLID_EXT,
}
