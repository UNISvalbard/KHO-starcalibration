% function definition to be added here

% For equidistant projection, we have a simple formula d=k*z where
%   d is the distance from the centre
%   z is the zenith angle
%   k is a calibration coefficient
%
% The calibration can be formulated as a minimisation problem where
% there are four unknowns
% - the pixel row and column coordinates of the zenith in the image
% - the calibration coefficient k
% - a rotation angle about the zenith direction (clockwise/anticlockwise)
%
% Estimate rough zenith location by assuming the imager is pointing
% straight up

zenithRow=size(img,1)/2;
zenithCol=size(img,2)/2;

% and then do a rough estimate of the calibration coefficient
% using the zenith angle and (pixel) distances from the zenith to
% the stars

starZen=90-starAlt;

starDistances= sqrt( (zenithRow-starRow).^2 + ...
    (zenithCol-starCol).^2);

k=mean(starDistances./starZen);

% Start with a rough rotation angle based on the offset between
% the geodetic and geomagnetic north direction at KHO

rotAngle=25*pi/180;