function stk2tiffDirs(varargin)
% stktiffDirs splits TIFF stacks (*.tif or *.stk) in input directory into folders with TIFF files.
% 
% Synopsis:    stk2tiffDirs(path)
%              stk2tiffDirs(path,'Crop','on')
%
% Input:
%      path : optional - path to the directory containing STK files. Can be
%      a string or a cell array of strings. If not input, the user will be
%      asked to select a folder.
%
%      Optional parameter/value pairs
%           Crop ('on'/'off') :  if 'on', a window will open asking the
%           user to select the region to crop.
%
% Francois Aguet, 09/01/2010


ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('path', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParamValue('Crop', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('Channels', 1, @isscalar);
ip.addParamValue('ChannelOrder', 'interleaved', @(x) any(strcmpi(x, {'interleaved', 'consecutive'})));
ip.addParamValue('ChannelNames', cell(1));
ip.addParamValue('WriteMode', 'overwrite', @(x) any(strcmpi(x, {'overwrite', 'append'})));
ip.addParamValue('Compression', 'none', @(x) any(strcmpi(x, {'none', 'lzw'})));
ip.parse(varargin{:});
stkpath = ip.Results.path;
crop = ip.Results.Crop;
nc = ip.Results.Channels;
cdir = ip.Results.ChannelNames;
if nc>1
    if isempty(cdir{1})
        cdir = arrayfun(@(i) ['c' num2str(i) filesep], 1:nc, 'UniformOutput', false);
    else
        cdir = cellfun(@(i) [i filesep], cdir, 'UniformOutput', false);
    end
end

if isempty(stkpath)
   stkpath = uigetdir('Select directory containing the STK files:'); 
   if (stkpath == 0)
       return;
   end
end

% Recursive call if input is cell array
if iscell(stkpath), 
    cellfun(@(x)stk2tiffDirs(x,'Crop',crop),stkpath);
    return
end

stkpath = [stkpath filesep];
stkList = dir(stkpath);
stkList = {stkList(~[stkList.isdir]).name};
idx = ~cellfun(@isempty, regexpi(stkList, '(\.tiff?|\.stk)$'));
stkList = stkList(idx);

N = length(stkList);
if N==0
    fprintf('No TIFF files found in input directory.\n');
end

for k = 1:N
    fprintf('Converting: %s\n', stkList{k});
    [~,stkname] = fileparts(stkList{k});
    stack = stackRead([stkpath stkList{k}]);
    
    if strcmpi(ip.Results.Crop, 'on')
       h = figure;
       imagesc(stack(:,:,1)); colormap(gray(256)); axis image;
       hr = imrect();
       pos = round(wait(hr));
       close(h);
       stack = stack(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
    end
    
    nf = size(stack,3) / nc;
    if strcmpi(ip.Results.ChannelOrder, 'interleaved')
        idx = @(c,f) c+(f-1)*nc;
    else % consecutive channels
        idx = @(c,f) f+(c-1)*nf;
    end
    
    fmt = ['%0' num2str(length(num2str(nf))) 'd'];
    for c = 1:nc
        destdir = [stkpath stkname filesep cdir{c}];
        [~,~] = mkdir(destdir);
        if strcmpi(ip.Results.WriteMode, 'Overwrite')
            for f = 1:nf
                imwrite(stack(:,:,idx(c,f)), [destdir stkname '_' num2str(f, fmt) '.tif'],...
                    'tif', 'Compression', ip.Results.Compression);
            end
        else
            for f = 1:nf
                imwrite(stack(:,:,idx(c,f)), [destdir stkname '.tif'],...
                    'tif', 'WriteMode', 'append', 'Compression', ip.Results.Compression);
            end
        end
    end
end