function C = pcombs(v,keepDiag)
% This function generate all pairwise combinations of a vector v of length
% n. The number of possible pairwise combinations is n x (n - 1) / 2 with
% keepDiag == false (default), n x (n + 1) / 2 otherwise. The elements of v
% must be all different.
%
% e.g.:
% pcombs([1 2 4 3], false)
% ans = 
%   1  2
%   1  4
%   1  3
%   2  4
%   2  3
%   4  3
%
% pcombs([1 2 4 3], true)
% ans =
%   1  1
%   1  2
%   1  4
%   1  3
%   2  2
%   2  4
%   2  3
%   4  4
%   4  3
%   3  3

if nargin < 2 || isempty(keepDiag)
    keepDiag = false;
end

assert(islogical(keepDiag));

assert(length(v) == length(unique(v)));

[I J] = meshgrid(v, v');
k = ~keepDiag;
I = triu(I, k);
J = triu(J, k);
C = horzcat(nonzeros(J), nonzeros(I));

end