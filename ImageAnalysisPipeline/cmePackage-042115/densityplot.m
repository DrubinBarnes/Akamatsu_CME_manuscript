%[ha, dRange] = densityplot(x, y, varargin) generates a 2D histogram / density plot
% from the input vectors
%
% Inputs:
%          x : sample vector for the 1st dimension
%          y : sample vector for the 2nd dimension
%
% Optional inputs:
%         xv : bin centers for the 1st dimension
%         yv : bin centers for the 2nd dimension
%
% Parameters ('specifier', value)
%           'Parent' : handle of the parent axes
%  'DisplayFunction' : intensity mapping, i.e., @sqrt

% Francois Aguet, 02/19/2013

function [ha, dRange] = densityplot(x, y, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('xv', linspace(min(x), max(x), 100));
ip.addOptional('yv', linspace(min(y), max(y), 100));
ip.addParamValue('Parent', gca, @ishandle);
ip.addParamValue('Div', 1, @isscalar);
ip.addParamValue('NormX', false, @islogical);
ip.addParamValue('DisplayFunction', @(x) x, @(x) isa(x, 'function_handle'));
ip.parse(varargin{:});
xv = ip.Results.xv(:)';
yv = ip.Results.yv(:)';
ha = ip.Results.Parent;

M = hist3([y(:) x(:)], {yv, xv});
if ip.Results.NormX
    M = M./repmat(sum(M,1), [size(M,1) 1]);
end

M = ip.Results.DisplayFunction(M/ip.Results.Div);

imagesc(xv, yv, M, 'Parent', ha);
axis(ha, 'xy', 'tight');
dRange = [min(M(:)) max(M(:))];
