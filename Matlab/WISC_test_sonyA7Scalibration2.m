% Testing the star calibration with the Sony images when selecting the
% stars in WISC
% 1) Read the image and to define the locations of stars in the image
% 2) Run the calibration parameter optimisation routine
% 3) Plot results

% This routine uses equisolid projection instead of equidistant

clear

filename=fullfile('..','Test_images/','LYR-Sony-111220_200933.jpg');
img=im2double(imread(filename));

[~, name, ext] = fileparts(filename);

imshow(sqrt(img))
axis on
title({[name ext],'Equisolid projection'})

%--------------------------------------------------
% Import the selected stars from WISC

T = readtable('/Users/guillaume/Code/KHO-starcalibration/selected_stars_WISC/selected_stars_LYR_20201211.tsv', 'FileType', 'text', ...
    'Delimiter', '\t');
starNames = T.star';
starAz    = T.azimuth_deg';
starAlt   = T.altitude_deg';
starCol   = T.star_col_px_1based';
starRow   = T.star_row_px_1based';


%---------------------------------------------------
% Do the star calibration

[zenithRow, zenithCol, k, rotAngle]= ...
    starcalibration2(size(img),starAlt,starAz, starRow, starCol);

fprintf('-------------------------------------------------------\n')
fprintf('        Zenith = (%.1f,%.1f)\n',zenithRow,zenithCol);
fprintf('             k = %.1f \n',k);
fprintf('           k/2 = %.1f \n',k/2);
fprintf('Rotation angle = %.2f (%.1f deg)\n', rotAngle, rotAngle*180/pi); 


%------------------------------------------------------------------
% Do a quick plot of the final results.

hold on
plot(zenithCol,zenithRow,'go')
plot(starCol, starRow, 'ro','markersize',10)
text(starCol+15, starRow, starNames,'color','r')

% Compute where the stars are in the image based on the calibration
% parameters to visualise the calibration results

theta=starAz*pi/180;
z=(90-starAlt)*pi/180;
d=k*sin(z/2);
newStarRow=zenithRow-d.*cos(theta+rotAngle);
newStarCol=zenithCol-d.*sin(theta+rotAngle);

plot(zenithCol,zenithRow,'go')
for i=1:length(starAz)
    plot(newStarCol(i),newStarRow(i),'g+')
end

% Plot a thick long line towards north and east

d=k*sin(45*pi/180/2); % 45 degrees from the horizon

northRow=zenithRow-d*cos(rotAngle);
northCol=zenithCol-d*sin(rotAngle);
plot([zenithCol northCol],[zenithRow northRow],'b','linewidth',3)
text(northCol,northRow-15,'NORTH','color','b')

eastRow=zenithRow-d*cos(pi/2+rotAngle);
eastCol=zenithCol-d*sin(pi/2+rotAngle);
plot([zenithCol eastCol],[zenithRow eastRow],'b','linewidth',3)
text(eastCol,eastRow+15,'EAST','color','b')

hold off

% Quick error analysis based on true and predicted star locations

rmsError=sqrt(sum((starRow-newStarRow).^2+(starCol-newStarCol).^2)...
        /length(starCol));

fprintf('Root Mean Square Error = %.2f pixels\n', rmsError);