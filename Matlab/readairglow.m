function img=readairglow(filename)
%
% img=readairglow(filename)
%
% The airglow files are raw binary files comprising
%   - a header (60 bytes)
%   - Dimension 1 as a 16-bit unsigned integer
%   - Dimension 2
%   - Size of each value in bytes (i.e. 2 for a 16-bit unsigned integer)
%   - Array (Dimension1 x Dimension2)
%

fid=fopen(filename);
if fid==-1
    error('readairglow: Cannot open file %s',filename);
end
[~,headercount]=fread(fid,60,'uint8');
if headercount~=60
    error('readairglow: could not read the header in %s',filename)
end

[nx, nxcount]=fread(fid,1,'uint16');
if nxcount~=1
    error('readairglow: could not read the first dimension in %s',filename)
end

[ny, nycount]=fread(fid,1,'uint16');
if nycount~=1
    error('readairglow: could not read the second dimension in %s',filename)
end

% The size
[~, ncount]=fread(fid,1,'uint16');
if ncount~=1
    error('readairglow: could not read the size of data in %s',filename)
end

[img, imgcount]=fread(fid,[nx ny],'uint16');
if imgcount~=nx*ny
    error('readairglow: read %d instead of %d values from %s', ...
        filename, imgcount, nx*ny)
end

fclose(fid);

img=img'; % Rotate the image to get north up, east to the left
end

