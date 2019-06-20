
% Francois Aguet, 06/30/2011

function lambda = name2wavelength(name)

s = getFluorPropStruct();
lv = [s.lambda_em];

if ischar(name)
    name = {name};
end

lambda = cellfun(@(x) lv(strcmpi(x, {s.name})), name, 'UniformOutput', false);

invalid = cellfun(@isempty, lambda);
if any(invalid)
    invalid = name(invalid);
    invalid = cellfun(@(x) [x ' '], invalid, 'UniformOutput', false);
    error(['Unsupported fluorophores: ' invalid{:}]);
end

lambda = [lambda{:}];
