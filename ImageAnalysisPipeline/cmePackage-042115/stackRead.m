function [stack, stackinfo] = stackRead(stackpath)
%
% [stack, stackinfo] = stackRead(stackpath)
%
% Wrapper function for François Nedelec's tiffread.m
% This works on STK files as well as multipage TIFFs.
% Outputs:
%
%   stack     : 3D data array
%   stackinfo : Any other information included in original file
%
% Francois Aguet, 01/2010

if (nargin == 0 || isempty(stackpath))
    stackinfo = tiffread();
else
    stackinfo = tiffread(stackpath);
end

stack = cat(3,stackinfo(:).data);
stackinfo = rmfield(stackinfo, 'data');