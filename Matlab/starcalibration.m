function [zenithRow, zenithCol, k, rotAngle]= ...
    starcalibration(img, starAlt, starAz, starRow, starCol)
%
% Estimation of camera calibration parameters for auroral all-sky
% cameras using fish-eye lenses
%
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
% A well-aligned all-sky imager is pointing straight up towards the
% zenith. So, a good estimate for the zenith pixel coordinates is
% in the centre of the image.
%
% For the all-sky images, the upper left corner row and column
% coordinates are (1,1). Also, we assume that the image is oriented
% similarly to most sky charts without any "left-right" mirroring. So,
% if north is up, then east is left. Or if north is to the right, then east
% is up.
%
% [zenithRow, zenithCol, k, rotAngle]= ...
%    starcalibration(img, starAlt, starAz, starRow, starCol)
%
% The input parameters are
%    img      all-sky image
%
%    starAlt  altitudes (above horizon) of the selected stars
%    starAz   azimuths
%
%    starRow  pixel rows for selected stars
%    starCol  pixel columns
%
% The outputs are
%    zenithRow the pixel row of the zenith (location)
%    zenithCol the pixel column for the zenith
%
%    k         the coefficient [pixels/degree]
%
%    rotAngle  the rotation about the optical axis between the
%              north and the horizontal axis (pixel columns) of the image.
%              A positive value means that if you rotate the image right,
%              then the image will be aligned north-south with north being
%              up (if the upper left corner is at (1,1) in pixel
%              coordinates)
%

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

% Given the initial estimates, find the values that
% result in smallest sum of position errors

mystartparams=[zenithRow zenithCol k rotAngle];
myconstants=[starAlt starAz starRow starCol];

% Create a helper function to be able to supply contants to
% fminsearch, see Matlab documentation for details.

fun=@(x)checkParameters(x,starAlt,starAz,starRow,starCol);

x=fminsearch(fun,mystartparams);

zenithRow=x(1);
zenithCol=x(2);
k=x(3);
rotAngle=x(4);



end

%===============================================================
% We minimise the error between the star positions (row, column)
% in the image and the expected locations computed from the
% apparent positions of the stars in the sky. The return
% value is simply the sum of squared distances.

function x=checkParameters(myparams,starAlt,starAz,starRow,starCol)
zenithRow=myparams(1);
zenithCol=myparams(2);
k=myparams(3);
rotAngle=myparams(4);

theta=starAz*pi/180;
d=k*(90-starAlt);
newStarRow=zenithRow-d.*cos(theta+rotAngle);
newStarCol=zenithRow-d.*sin(theta+rotAngle);
x=sum(((starRow-newStarRow).^2+(starCol-newStarCol).^2));

end