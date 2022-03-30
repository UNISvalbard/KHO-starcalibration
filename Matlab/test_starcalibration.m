% Testing the star calibration
% 1) Run the file reading test to read the image and to define the
%    locations of stars in the image
% 2) Run the calibration parameter optimisation routine
% 3) Plot results

test_readairglow;

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