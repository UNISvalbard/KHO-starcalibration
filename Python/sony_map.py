"""
Map an all-sky image onto a fixed-altitude latitude/longitude grid.

Python translation of the MATLAB script.
Requires: numpy, matplotlib, scipy, Pillow
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.ndimage import map_coordinates
from scipy.interpolate import RegularGridInterpolator


# -----------------------------------------------------------
# Configuration
Ngrid = 300
max_zenith_angle = 70.0  # degrees
auroral_h = 140.0        # km

filename = Path("..") / "Test_images" / '2020-01-03' / "LYR-Sony-20200103_085735.jpg"
# filename = Path("..") / "Test_images" / "LYR-Sony-211119_061942.jpg"

# -----------------------------------------------------------
# Read image as floating-point RGB values in [0, 1]
img = np.asarray(Image.open(filename).convert("RGB"), dtype=np.float64) / 255.0
image_height, image_width, _ = img.shape

# -----------------------------------------------------------
# Calibration coefficients

MODEL = "equisolid"
zenithRow = 247.2
zenithCol = 243.6
k = 330.8        # pixels
rotAngle = 0.50       # radians
p = 0 # needed for EQUISOLID_EXT lens model

glat = 78.20            # imager latitude, degrees
glon = 15.70            # imager longitude, degrees

# -----------------------------------------------------------
# Geometry
Re = 6378.1             # Earth radius, km

lats = np.linspace(glat - 2.55, glat + 2.55, Ngrid)
lons = np.linspace(glon - 12.6, glon + 12.6, Ngrid)

theta_0 = np.deg2rad(glat)
fii_0 = np.deg2rad(glon)

r_0 = Re * np.array([
    np.cos(theta_0) * np.cos(fii_0),
    np.cos(theta_0) * np.sin(fii_0),
    np.sin(theta_0),
])

n_0 = r_0 / np.linalg.norm(r_0)

e_theta = np.array([
    -np.sin(theta_0) * np.cos(fii_0),
    -np.sin(theta_0) * np.sin(fii_0),
    np.cos(theta_0),
])

e_fii = np.array([
    -np.sin(fii_0),
    np.cos(fii_0),
    0.0,
])

# MATLAB meshgrid(lats, lons) gives arrays shaped (len(lons), len(lats)).
# Keep the same convention so the final transpose matches the MATLAB script.
theta_mesh, fii_mesh = np.meshgrid(np.deg2rad(lats), np.deg2rad(lons))
theta = theta_mesh.ravel(order="F")
fii = fii_mesh.ravel(order="F")
n_points = theta.size

# Position vectors, shape (3, n_points)
r = (Re + auroral_h) * np.vstack((
    np.cos(theta) * np.cos(fii),
    np.cos(theta) * np.sin(fii),
    np.sin(theta),
))

# -----------------------------------------------------------
# Azimuth using broadcasting
dr = r - r_0[:, None]

X = np.sum(dr * e_theta[:, None], axis=0)
Y = np.sum(dr * e_fii[:, None], axis=0)
azimuth = np.arctan2(X, Y)

# -----------------------------------------------------------
# Zenith angle
v_norm = np.linalg.norm(dr, axis=0)
v = np.divide(
    dr,
    v_norm[None, :],
    out=np.full_like(dr, np.nan),
    where=v_norm[None, :] > 0,
)

cos_z = np.sum(v * n_0[:, None], axis=0)
z = np.arccos(np.clip(cos_z, -1.0, 1.0))

# Mask points outside the allowed zenith angle
z[z > np.deg2rad(max_zenith_angle)] = np.nan

# Image pixel coordinates (MATLAB uses 1-based coordinates; Python uses 0-based)
if MODEL == 'equidistant':
    d = k * z
elif MODEL == 'equisolid':
    d = k * np.sin(z / 2)  
elif MODEL == 'equisolid_ext':
     d = k * np.sin(a * z) + p
else:
    print('Need to choose valid lens model')
    
i = (zenithCol - d * np.sin(azimuth - rotAngle)) - 1.0
j = (zenithRow - d * np.cos(azimuth - rotAngle)) - 1.0

# -----------------------------------------------------------
# Interpolate each RGB channel at (row=i, col=j).
# scipy.ndimage.map_coordinates uses zero-based pixel coordinates.
coords = np.vstack((i, j))
mapped_channels = []

for channel in range(3):
    sampled = map_coordinates(
        img[:, :, channel],
        coords,
        order=1,             # linear interpolation
        mode="constant",
        cval=np.nan,
        prefilter=False,
    )
    sampled[~np.isfinite(z)] = np.nan
    mapped = sampled.reshape(theta_mesh.shape, order="F").T
    mapped_channels.append(mapped)

mapped_asc = np.stack(mapped_channels, axis=-1)

np.save("mapped_asc.npy", mapped_asc)

# -----------------------------------------------------------
# Plots
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Original image with sample points
axes[0].imshow(img, origin="upper")
axes[0].plot(j, i, ".", color="red", markersize=1)
axes[0].set_title("ASC image and sample points")
axes[0].set_xlabel("Pixel coordinates")
axes[0].set_ylabel("Pixel coordinates")
axes[0].set_aspect("equal")

# Mapped image
axes[1].imshow(
    mapped_asc,
    extent=[lons.min(), lons.max(), lats.min(), lats.max()],
    origin="lower",
    aspect="auto",
)
axes[1].plot(glon, glat, "ro", markersize=8)
axes[1].set_xlabel("Longitude (degrees)")
axes[1].set_ylabel("Latitude (degrees)")
axes[1].set_title("Mapped image")

# Reference image rotated by -rotAngle, cropped to original dimensions
from scipy.ndimage import rotate
imgrot = rotate(
    img,
    angle=-np.rad2deg(rotAngle),
    axes=(1, 0),
    reshape=False,
    order=1,
    mode="constant",
    cval=0.0,
)
axes[2].imshow(np.fliplr(imgrot), origin="upper")
axes[2].set_title("Reference image")
axes[2].set_aspect("equal")

plt.tight_layout()
plt.show()
