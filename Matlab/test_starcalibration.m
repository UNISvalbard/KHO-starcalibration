% Testing the star calibration
% 1) Run the file reading test
% 2) Find the centre of the image (which should be close the zenith)
% 3) Use the known stars/planets to estimate the calibration parameters

test_readairglow;

starcalibration; % Replace with a proper function call

fprintf('Zenith = (%.1f,%.1f)\n',zenithRow,zenithCol);
fprintf('     k = %.1f\n',k);
fprintf('offset = %.1f (%.1f deg)\n', rotAngle, rotAngle*180/pi); 

% Do a plot of the final results

hold on
plot(zenithCol,zenithRow,'go')

for i=1:length(starAz)
    d=k*(90-starAlt(i));
    myAzimuth=starAz(i)*pi/180;
    i=-d*cos(myAzimuth+rotAngle);
    j=-d*sin(myAzimuth+rotAngle);
    plot(zenithCol+[0 j],zenithRow+[0 i],'g')
    plot(zenithCol+j, zenithRow+i,'go')
end

hold off
