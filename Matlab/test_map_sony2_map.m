% Map an all-sky image into a fixed altitude for map projections
%
% We use equisolid lens projection
%
% For later overlaying on top of geographical maps etc. use
% pcolor(lons,lats,mapped_asc)
%
% Note: this script uses NaN's to handle out of range etc. values

clear

%-----------------------------------------------------------
% The mapping is done from the all-sky image onto a regular
% latitude and longitude grid. 
%
% The auroral altitude is depatable depending on activity
% and these are just suggestions used in literature
% - nighttime 110km or 125km
% - daytime 220km

Ngrid=300; % Number of grid points along lat and lon
max_zenith_angle=80; % The limit for the zenith angle
auroral_h= 110; % Auroral altitude in kilometers 

% Equal number of grid points both in lat and lon
% - the extents in latitude and longitude were manually
%   selected to cover as much as possible within the field-of-view
%   while still having nice integer numbers for latitude and longitude

% Nightside aurora at 110km
lats=linspace(73, 83, Ngrid);
lons=linspace(-8,39,Ngrid);

% Nightside aurora at 125km
%lats=linspace(73, 84, Ngrid);
%lons=linspace(-10,42,Ngrid);

% Dayside aurora at 220km
%lats=linspace(70, 87,Ngrid);
%lons=linspace(-28, 59,Ngrid);


%----------------------
% Read the image
filename=fullfile('..','Test_images/','LYR-Sony-211121_170751.jpg');
%filename=fullfile('..','Test_images/','LYR-Sony-211119_061942.jpg');
%filename=fullfile('..','Test_images/','LYR-Sony-20260110_083148.jpg'); %
%filename=fullfile('..','Test_images/','LYR-Sony-20251222_201712.jpg');

img=im2double(imread(filename));
%img=sqrt(img(:,:,1));

%-----------------------------------------------
% Calibration parameters
% - these are from test_sonyA7Scalibration2.m
% d = k sin(z/2)

zenithRow=1452.2;
zenithCol=1424.8;
k=1949.4; 
rotAngle=0.4900; 

glat = 78.20; % Imager latitude
glon = 15.70; % Imager longitude



%-----------------------------------------------------------------------
% For each point (lat, lon) at the auroral altitude, we determine
% the azimuth and zenith angle as seen from the station. Based on these,
% the corresponding pixels in the image are determined using the output
% from geometric calibration.

Re=6378.1; % Radius of the Earth (km)


%----------------------------------------------------------------
% First, determine a tangent plane at the station.
%  - the position vector to the station (r_0)
%  - the plane's normal (unit) vector (n_0)
%  - unit vectors in latitude and longitude (e_theta,e_fii)
%  - position vectors for all grid points (r)

theta_0=glat*pi/180;
fii_0=glon*pi/180;

r_0=Re * [cos(theta_0)*cos(fii_0); ...
    cos(theta_0)*sin(fii_0); ...
    sin(theta_0)];

n_0=r_0/norm(r_0);

e_theta=[-sin(theta_0)*cos(fii_0); ...
    -sin(theta_0)*sin(fii_0); ...
    cos(theta_0)];

e_fii=[-sin(fii_0); cos(fii_0); 0];

% Form the vectors to the lat,lon grid points
% - first mesh points from the lat,lon grids
% - then mesh point matrices into vectors

[theta_mesh,fii_mesh]=meshgrid(lats*pi/180,lons*pi/180);

theta=reshape(theta_mesh,[],1); % Convert to vectors for following
fii=reshape(fii_mesh,[],1);     % calculations
n_points=length(theta);

r=(Re+auroral_h) * [cos(theta).*cos(fii), ... % Position vectors (columns)
    cos(theta).*sin(fii), ...                 % to all grid points
    sin(theta)]';

% Project the locations onto a plane (e_theta,e_fii)
% for azimuth calculations
% - the "repmat" operations allow vectorised execution
%   of later calculations

r_0_rep=repmat(r_0,1,n_points);
e_theta_rep=repmat(e_theta,1,n_points);
e_fii_rep=repmat(e_fii,1,n_points);
n_0_rep=repmat(n_0,1,n_points);

% - plane coordinates using the dot product

X=dot(r-r_0_rep,e_theta_rep,1);
Y=dot(r-r_0_rep,e_fii_rep,1);
azimuth=atan2(X,Y);

% Obtain the zenith angles for each point
% - the angle is the one between the vector v=r-r_0 from
%   the station to the auroral point and the plane normal vector
%
% - the position vector to the station would work as well as
%   the plane normal vector as they are in the same direction

% Angle between vectors
% - use the dot product equation:  v (dot) w = |v| |w| cos(angle)

v=r-r_0_rep;
v=bsxfun(@rdivide,v,sqrt(sum(v.^2,1))); % Each column into unit vector

% The zenith angles between unit vectors and plane normal
z=acos(dot(v,n_0_rep)); 

% Ignore pixels with too large a zenith angle
z(z>max_zenith_angle*pi/180)=NaN;

% Find the corresponding pixel locations (i,j) in the all-sky image
% Use equisolid projection

d=k*sin(z/2);

i=zenithCol-d.*sin(azimuth-rotAngle);
j=zenithRow-d.*cos(azimuth-rotAngle);

%------------------------------------------------------------
% Interpolate intensitites from the image
% - convert the vector to image matrix
% - out of matrix values will be NaN
% - reorient to match the lat,lon grid

mapped_r=reshape(interp2(img(:,:,1),j,i,'linear',NaN),size(theta_mesh));
mapped_r=mapped_r';

mapped_g=reshape(interp2(img(:,:,2),j,i,'linear',NaN),size(theta_mesh));
mapped_g=mapped_g';

mapped_b=reshape(interp2(img(:,:,3),j,i,'linear',NaN),size(theta_mesh));
mapped_b=mapped_b';

mapped_asc=zeros(Ngrid,Ngrid,3);
mapped_asc(:,:,1)=mapped_r;
mapped_asc(:,:,2)=mapped_g;
mapped_asc(:,:,3)=mapped_b;


% The end of the map projection part of the script...

%*******************************************************************
% Sanity check: show the original image with sample points
% and the resulting interpolated mapped image with the station
% location overlayed.

% 

subplot(1,2,1)
image(img)
axis image
hold on
plot(j,i,'.','color','red','markersize',1)
hold off
title('ASC image and sample points')
xlabel('Pixel coordinates')
ylabel('Pixel coordinates')


subplot(1,2,2)
image(lons,lats,mapped_asc)
axis xy
xlabel('Longitude')
ylabel('Latitude')

