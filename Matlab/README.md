#  Estimation of camera calibration parameters with Matlab

## Equidistant projection

For equidistant projection, we have a simple formula `d=k*z` where
- `d` is the distance from the centre
- `z` is the zenith angle
- `k` is a calibration coefficient

A well-aligned all-sky imager is pointing straight up towards the
zenith. However, the image may not be aligned in north-south
direction, so there is usually an offset in azimuth (i.e. rotation about the optical axis) that needs to be estimated as well.

For the all-sky images, the upper left corner row and column
coordinates are (1,1). Also, we assume that the image is oriented
similarly to most sky charts without any "left-right" mirroring. So,
if north isup, then east is left. Or if north is to the right, then east
is up.

The Matlab-function for estimating the parameters is
```
 [zenithRow, zenithCol, k, rotAngle]= ...
    starcalibration(img, starAlt, starAz, starRow, starCol)
```

| Input parameter | Description | Type |
| --------------- | ------------ | --- |
| `imgsize` | The size of the all-sky image | *double* (2x1) |
| `starAlt` |  The altitudes (above horizon) of the selected stars | *double* (Nx1) |
| `starAz` | The azimuths of the stars | *double* (Nx1) |
| `starRow` | The pixel rows for selected stars | *double* (Nx1) |
| `starCol` | The pixel columns for the stars |  *double* (Nx1) |

| Output parameter | Description | Type |
| --------------- | ------------ | --- |
| `zenithRow` | The pixel row of the zenith (location) | *double* |
| `zenithCol` | the pixel column for the zenith | *double* |
| `k` | The coefficient [pixels/degree] | *double* |
| `rotAngle` | The rotation angle (in radians) | *double* |

The rotation angle is  about the optical axis between the
north and the vertical axis (pixel columns) of the image.
A positive value means that if you rotate the image right,
then the image will be aligned north-south with north being
up.

## How to do the calibrations

There are two example calibration for an airglow imager and an all-sky
camera.You can use one as a template for a new camera and then 
update and/or add more stars for the calibration.


