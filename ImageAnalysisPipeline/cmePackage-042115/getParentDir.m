% [cpath] = getParentDir(cpath) Returns the parent directory's path from the input path

% Francois Aguet, 11/05/2010

function cpath = getParentDir(cpath, level)
if nargin<2
    level = 1;
end

fsIdx = regexp(cpath, filesep);

if strcmp(cpath(end), filesep)
    cpath = cpath(1:fsIdx(end-level));
else
    cpath = cpath(1:fsIdx(end-level+1));
end