function makeCMEanalysisFolders_twocolor(varargin)
%

% This program takes movies in a single folder and puts them in separate
% folders, with subfolders named in a format that cmeAnalysis can read.

% If you use Sun Hong's imageRegistration.py program, this program will
% separate into one expt folder with separate folders for each field (and
% each channel)

% Matt Akamatsu, code borrowed from stk2tiffDirs.m

% from .nd2 fild

% split into 2+ channels

% save as GFP RFP etc with right name. 

%stktiffDirs splits TIFF stacks (*.tif or *.stk) in input directory into folders with TIFF files.
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

timestep = input('What the time interval in seconds? ');
nbChannels = input('How many channels? ');

% nbFields = input('How many fields? '); % could automate later


% 
% if length(timestep) == 1
%     
%     timesteps = repmat(timestep, nbFields, 1);
% 
% else
%     timesteps = timestep;
% end

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('path', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParamValue('Crop', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('Channels', 1, @isscalar);
ip.addParamValue('ChannelOrder', 'interleaved', @(x) any(strcmpi(x, {'interleaved', 'consecutive'})));
ip.addParamValue('ChannelNames', cell(1));
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
if iscell(stkpath) 
    cellfun(@(x)stk2tiffDirs(x,'Crop',crop),stkpath);
    return
end

stkpath = [stkpath filesep];
stkList = dir(stkpath);
stkList = {stkList(~[stkList.isdir]).name};
idx = ~cellfun(@isempty, regexpi(stkList, '(\.tiff?|\.stk|\.nd2)$'));
stkList = stkList(idx);

N = length(stkList);
if N==0
    fprintf('No TIFF files found in input directory.\n');
end

% for k = 1:N
%     fprintf('Converting: %s\n', stkList{k});
    [~,stkname] = fileparts(stkList{1});
%     stackTemp = stackread(
%     stack = stackRead([stkpath stkList{k}]);
    
%     if strcmpi(ip.Results.Crop, 'on')
%        h = figure;
%        imagesc(stack(:,:,1)); colormap(gray(256)); axis image;
%        hr = imrect();
%        pos = round(wait(hr));
%        close(h);
%        stack = stack(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
%     end
%     
%     nf = size(stack,3) / nc;
%     if strcmpi(ip.Results.ChannelOrder, 'interleaved')
%         idx = @(c,f) c+(f-1)*nc;
%     else % consecutive channels
%         idx = @(c,f) f+(c-1)*nf;
%     end
%     
%     fmt = ['%0' num2str(length(num2str(nf))) 'd'];
    
    stknameTree = stkname(1:end-8)

%     colors = {'left' 'right'};
    
    %%% here make it flexible how many channels are being added. read in from metadata in nd2 file.
    colorCode = {'GFP' 'RFP' 'farRed'};
    cd(stkpath)
    allFiles= dir('*tif');
    nbFields=ceil(length(allFiles)/nbChannels)
    kk = 1;
    
    if length(timestep) == 1
    
    timesteps = repmat(timestep, nbFields, 1);

    else
        timesteps = timestep;
    end
    
    for m = 1:nbFields
 
        destdirCell = [stkpath stknameTree filesep 'cell' num2str(m) '_' num2str(timesteps(m)) 's'];
        [~,~] = mkdir(destdirCell);

        for n = 1:nbChannels
        
            fileEnding = [colorCode{n} '_0' num2str(m)];
        
            %disp(fileEnding);
            
            destdirColor = [destdirCell filesep colorCode{n}];
            
            disp(destdirColor);
            [~,~] = mkdir(destdirColor);
            movefile([stkpath stkList{kk}], destdirColor)

            kk = kk+1;
        end
    end
            
%     destdir = [stkpath stkname filesep 'cell1_' num2str(timestep) 's'];
%        % for f = 1:nf
%             %imwrite(stack(:,:,idx(c,f)), [destdir stkname '_' num2str(f, fmt) '.tif'], 'tif');
%             movefile([stkpath stkList{k}], destdir)
        %end
    
%     for c = 1:nc
%         destdir = [stkpath stkname filesep cdir{c}];
%         [~,~] = mkdir(destdir);
%         for f = 1:nf
%             %imwrite(stack(:,:,idx(c,f)), [destdir stkname '_' num2str(f, fmt) '.tif'], 'tif');
%             movefile
%         end
% end
    
disp('Done!')
end