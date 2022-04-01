# Star calibration of all-sky images

All-sky auroral imagers are used in studying the physical processes that take place in the upper atmosphere. The name "all-sky" refers to an ultra wide-angle field-of-view: one image captures the whole sky. Typically one uses so called circular fish-eye lenses where you get a smaller circular image of the whole sky inside the image frame.

For more detailed analysis, one needs to know which way the instrument is "looking". In other words, which pixel in the image corresponds to the zenith (straight up), which way is north and east etc. In computer vision terminology, this is called camera calibration. However, in the context of auroral research, one usually uses the term "star calibration" to separate this procedure from *intensity* calibration, where the recorded intensitites are mapped to physical units.

## Basic procedure

Given an all-sky image, we have metadata about the date and time as well as the location of the imager. Based on this, we can calculate the positions of stars in the sky and this provides us with known targets.

1. Select an image where stars are clearly visible
2. Use the date, time and location to calculate star positions to produce a sky chart (e.g. astropy, Stellarium)
3. Identify the stars visible in the image
   - take a note of the azimuth and altitude of the stars apparent position
   - determine the pixel location (row, column) of the star in the image
5. Estimate the camera parameters based on the stars' locations
   - note that there are several different lens projections

## What you will find in this repo

- short scripts that produce a "pretty good" calibration
  - the lens is assumed to be perfect i.e. we use the ideal [mapping function](https://en.wikipedia.org/wiki/Fisheye_lens#Mapping_function)
  - calibration parameters that you can use to e.g. project data on maps for further analysis such as satellite trajectories across the field-of-views of the imagers
- example calibrations for selected instruments at [Kjell Henriksen Observatory](http://kho.unis.no)

## Notes

The short test scripts are there to do a quick functional test, but for an
actual scientific study, one should use a larger number of reference points (stars). However, even with only four stars, we can obtain a pretty good results that may be sufficient.

See e.g. [A detailed study of auroral fragments](http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-388546) by Joshua Dreyer.


## To do

- [ ] Convert Matlab-scripts to python
- [ ] Add a map projection example
