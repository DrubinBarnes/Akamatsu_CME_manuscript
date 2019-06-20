function indx=locmin1d(x,winSize)
%LOCMIN1D returns a list of local minima in vector x within a window.
%
% SYNOPSIS indx=locmin1d(x, n)
%
% INPUT 
%   x:  data vector
%
%   winSize:  (optional) length of the window. Must be odd. Default value
%   is 3
%
% OUTPUT indx : index list to local minima in x

% STARTED GD 29-Nov-1999
% MODIFIED SB 6-Apr-2010

if nargin > 1 && ~isempty(winSize)
    if mod(winSize,2) == 0
        error('winSize must be odd.');
    end
else
    winSize = 3;
end

hside = floor(winSize/2);
xPad = padarray(x(:), hside);
nPad = length(xPad);

xMin = arrayfun(@(pos) min(xPad(pos-hside:pos+hside)), hside+1:nPad-hside);

indx = find(x(:) == xMin(:));
