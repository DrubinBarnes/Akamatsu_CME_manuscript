%setErrorbarStyle(he, pos, de) modifies the width of the error bars

% Francois Aguet, Feb 22 2011 (last modif. 01/21/2012)

function setErrorbarStyle(he, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('he');
ip.addOptional('de', 0.2, @isscalar);
ip.addParamValue('Position', 'both', @(x) any(strcmpi(x, {'both', 'top', 'bottom'})));
ip.parse(he, varargin{:});
de = ip.Results.de;

he = get(he, 'Children');
xd = get(he(2), 'XData');
if strcmpi(ip.Results.Position, 'bottom')
    xd(4:9:end) = xd(1:9:end);
    xd(5:9:end) = xd(1:9:end);
else
    xd(4:9:end) = xd(1:9:end) - de;
    xd(5:9:end) = xd(1:9:end) + de;
end
if strcmpi(ip.Results.Position, 'top')
    xd(7:9:end) = xd(1:9:end);
    xd(8:9:end) = xd(1:9:end);
else
    xd(7:9:end) = xd(1:9:end) - de;
    xd(8:9:end) = xd(1:9:end) + de;
end
set(he(2), 'XData', xd);