function isit = isposint(in)
%ISPOSINT returns true if all elements of the input are positive integers
%
% tf = isposint(in)
%
% This function duplicates the functionality of isposintscalar and
% isposintmat, but those functions are in the system identification
% toolbox, so this avoids having your code depend on that toolbox.
%
% Hunter Elliott
% 4/2011
%
isit = isnumeric(in) & in > 0 & abs(round(in)) == in;