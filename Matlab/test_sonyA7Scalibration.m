% Testing the star calibration with the Sony images
% 1) Read the image and to define the locations of stars in the image
% 2) Run the calibration parameter optimisation routine
% 3) Plot results

clear
filename=fullfile('..','Test_images/','LYR-Sony-111220_200933.jpg');
img=imread(filename);

imshow(img)
axis on
title('Sony A7S 2020/12/11 20:09:33 UT')

% Since the image is practically taken at the same time as the airglow
% image, we are using the same stars and just update the pixel coordinates

starNames={'Vega','Capella','Dubhe','Deneb'};

starAz=[306+42/60+19.3/3600, ... % Stellarium shows the azimuth in degrees, 
    129+20/60+26.5/3600, ...     % minutes and seconds
    39+22/60+51.8/3600, ...
    280+15/60+56.1/3600];

starAlt=[32+21/60+35.9/3600, ...
    54+19/60+2.7/3600, ...
    53+24/60+34.4/3600, ...
    44+27.5/60+57.6/3600];

starRow=[599, 2005, 1219, 982];
starCol=[1822, 1196, 861, 2015];

hold on
plot(starCol, starRow, 'ro','markersize',10)
text(starCol+15, starRow, starNames,'color','r')
hold off

%---------------------------------------------------

[zenithRow, zenithCol, k, rotAngle]= ...
    starcalibration(img,starAlt,starAz, starRow, starCol);

fprintf('-------------------------------------------------------\n')
fprintf('        Zenith = (%.1f,%.1f)\n',zenithRow,zenithCol);
fprintf('             k = %.1f [pixel/deg]\n',k);
fprintf('Rotation angle = %.1f (%.1f deg)\n', rotAngle, rotAngle*180/pi); 


% Do a quick plot of the final results.
% Note that 'test_readairglow' has already plotted the image with 
% a number of identified stars, so we just need to add lines, markers etc.

hold on
plot(zenithCol,zenithRow,'go')

theta=starAz*pi/180;
d=k*(90-starAlt);
newStarRow=zenithRow-d.*cos(theta+rotAngle);
newStarCol=zenithRow-d.*sin(theta+rotAngle);

plot(zenithCol,zenithRow,'go')
for i=1:length(starAz)
    plot([zenithCol newStarCol(i)],[zenithRow newStarRow(i)],'g')
end

% Plot a thick long line towards north and east

d=k*45; % 45 degrees from the horizon

northRow=zenithRow-d*cos(rotAngle);
northCol=zenithCol-d*sin(rotAngle);
plot([zenithCol northCol],[zenithRow northRow],'b','linewidth',3)
text(northCol,northRow-15,'NORTH','color','b')

eastRow=zenithRow-d*cos(pi/2+rotAngle);
eastCol=zenithCol-d*sin(pi/2+rotAngle);
plot([zenithCol eastCol],[zenithRow eastRow],'b','linewidth',3)
text(eastCol,eastRow+15,'EAST','color','b')

hold off