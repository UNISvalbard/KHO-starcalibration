% Testing the airglow image reading routine
%
% The expected result is a plot of an airglow image
% where you can clearly see stars. As the original image uses
% 16-bit intensity values, we use the square root to compress the range to
% 8 bits for visualisation purposes.
%


filename=fullfile('..','Test_images/','W8446B18_2020-12-11_20-09-36.img');
img=readairglow(filename);

imagesc(sqrt(img))
axis image
colormap('gray')
title('KHO airglow imager 2020/12/11 20:09:36UT')

minvalue=min(img,[],'all');
maxvalue=max(img,[],'all');
fprintf("Min value read %d vs. expected value 51\n",minvalue);
fprintf("Max value read %d vs. expected value 64941\n",maxvalue);

%--------------------------------------------
% Label a couple of identified stars/planets

starNames={'Vega','Capella','Mars','Dubhe','Deneb'};

starAz=[306+42/60+19.3/3600, ... % Stellarium shows the azimuth in degrees, 
    129+20/60+26.5/3600, ...     % minutes and seconds
    201+54/60+13.9/3600, ...
    39+22/60+51.8/3600, ...
    280+15/60+56.1/3600];

starAlt=[32+21/60+35.9/3600, ...
    54+19/60+2.7/3600, ...
    19+6/60+1.5/3600, ...
    53+24/60+34.4/3600, ...
    44+27.5/60+57.6/3600];

starRow=[71, 361, 394,195, 150];
starCol=[340, 206, 428, 138, 380];

hold on
plot(starCol, starRow, 'ro','markersize',10)
text(starCol+15, starRow, starNames,'color','r')
hold off
