%% Extract XY coordinates and plot

% This program extracts XY coordinates from u-track and saves and plots
% them for each track. Note you need to use Matlab 2013b or later. 

% Matt Akamatsu, Drubin lab, December 2015
% modified 10/2016 to incorporate a second track rejection step from Josh
% Ferguson (OSU, Kural lab). This step removes tracks that have erratic intensity jumps between timepoints. 

% Select TrackFeatures file
function extractXYcoordinates_useCellMask(data)
disp('running extractXYcoords');
%clear;


%% second track rejection step from Josh Ferguson (OSU, Kural lab). 
% This step removes tracks that have erratic intensity jumps between timepoints.

rejectErraticTracks = true;

%% set "makefigs" to false if you don't want those figures

makefigs = false;
saveallgraphs = false;

% number of tracks to plot, for track positions (fig 1) and track durations (fig 5). set to -1 if you want to plot all tracks.
% Plots of lifetime distributions do not change based on this number.

nbTracksToPlot = 500; %-1

%% Choose cutoff for lifetime distribution. 
% The program will discard lifetimes less than this value (in seconds) if "cutLifetimes = true".

cutLifetimes = true;

lifetimesCutoffSeconds = 3


% If you ran cmeAnalysis and chose a region of the movie (cell mask), turn
% this on to keep only the tracks within the mask.

useCellMask = true;


% Later this program could scan through pre-defined folder structures...

%% Ask for file name
if ~exist('data', 'var')
    curfilename = input('What do you want to call this file? ', 's');
    
    pixel_size = input('What was the pixel size of the camera in nm? ');
    
    frame_time = input('How many seconds between each frame? ');
    
    load('previous_file.mat')
    [utrack_file, utrack_folder_name] = uigetfile('*.mat','Select the tracking_result.mat or trackedFeatures.mat file.', previous_file);
    if isdir(utrack_folder_name)
        previous_file = fullfile(utrack_folder_name, utrack_file);
        save('previous_file', 'previous_file')
    else
        disp('You did not pick a folder.')
        return
    end
else
    curfilename = data.curfilename;
    pixel_size = data.pixelSize * (1e9)/data.M;
    frame_time = data.framerate;
    utrack_file = 'trackedFeatures.mat';
    utrack_folder_name = strcat(data.source, 'Tracking/');
end

lifetimesCutoffFrames = ceil(lifetimesCutoffSeconds/frame_time)

utrack_full = fullfile(utrack_folder_name, utrack_file);
uTrackParent = open(utrack_full);

uTracks = uTrackParent.tracksFinal;

if useCellMask
    

    load('previous_file.mat')
    
    if (exist('data', 'var'))
        cellmask_file = 'cellmask.tif';
        cellmask_folder_name = strcat(data.source, 'Detection/');
    else
        [cellmask_file, cellmask_folder_name] = uigetfile('*.tif','Select the cellmask file.', previous_file);
    end
    if isdir(cellmask_folder_name)
        previous_file = fullfile(cellmask_folder_name, cellmask_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
    cellmask_full = fullfile(cellmask_folder_name, cellmask_file);
    cellmask = readtiff(cellmask_full); %or imread

% % flip Y 

% cellmaskUnflipped = cellmask;
% 
% cellmask = flipud(cellmask);
% % 
%     
%     separatorsMask = utrack_folder_name==pathsep;
%     [~, ]
%     
%    curCellMask = utrack_folder_name
    
end

% if isdir(ref_img_folder_name)
%     previous_file = fullfile(ref_img_folder_name, ref_img_file);
%     save('previous_file', 'previous_file')
% else
%     display('You did not pick a folder.')
%     return
% end

alltracks = {};
alltracksNoCutoff = {};
trackLengths = [];
trackLengthsNoCutoff = [];

% n is the current track number.

% calculate number of tracks (no cutoff)

nbTracksNoCutoff = length(uTracks);


% Identify the number of XY data points for each track

% NOTE will fail currently f you have splitting and merging tracks. it
% assumes you only have nomerging tracks.

% Find the x and y values for the first time point of each track, and save
% as xs and ys

% Keep ONLY if the first time point of the track is within the cell mask.

nn = 1;

for n = 1:nbTracksNoCutoff
    trackLengthsNoCutoff(n) = length(uTracks(n).tracksCoordAmpCG(1,:))/8;
    
    xs = [];
    ys = [];
    Is = [];
    Istd = [];
    
    for m = 1:trackLengthsNoCutoff(n)
        
        curIndex = (m-1)*8+1; %Find where X and Y are (indecies)
        
        % Just for the primary track. ignores track splits and merges.
        
        xs(m,1) = uTracks(n).tracksCoordAmpCG(1, curIndex);
        ys(m,1) = uTracks(n).tracksCoordAmpCG(1, curIndex+1);
        Is(m,1) = uTracks(n).tracksCoordAmpCG(1, curIndex+3);
        Istd(m,1) = uTracks(n).tracksCoordAmpCG(1, curIndex+7);
        
        
    end
    
    % This treatment of tracks ignores merging and splitting events. This
    % might cause problems later.
    
    trackend = uTracks(n).seqOfEvents(end,1);
    trackbeg = uTracks(n).seqOfEvents(1,1);
    
    firstX = round(xs(1));
    firstY = round(ys(1));
    
    firstTimePoint = [round(xs(1)), round(ys(1))];
    
    firstTimePoints(n,:) = firstTimePoint;
    curMaskState(n) = sum(cellmask(firstTimePoint));    
    timepointsAndMask(n,:) = [firstTimePoints(n) curMaskState(n)];
    %disp(cellmask(firstTimePoint));
    
    if useCellMask
    
%         if cellmask(firstTimePoint)
        if cellmask(firstY, firstX)

            trackduration = transpose(trackbeg:trackend);

            alltracksNoCutoff(nn).xs = xs;
            alltracksNoCutoff(nn).ys = ys;
            alltracksNoCutoff(nn).Is = Is;
            alltracksNoCutoff(nn).Istd = Istd;
            alltracksNoCutoff(nn).trackduration = trackduration;
            alltracksNoCutoff(nn).tracklength = trackLengthsNoCutoff(n);
            alltracksNoCutoff(nn).lifetime = alltracksNoCutoff(nn).tracklength * frame_time;

            nn = nn+1;
        end
    else
        trackduration = transpose(trackbeg:trackend);
        
        alltracksNoCutoff(n).xs = xs;
        alltracksNoCutoff(n).ys = ys;
        alltracksNoCutoff(n).Is = Is;
        alltracksNoCutoff(n).Istd = Istd;
        alltracksNoCutoff(n).trackduration = trackduration;
        alltracksNoCutoff(n).tracklength = trackLengthsNoCutoff(n);
        alltracksNoCutoff(n).lifetime = alltracksNoCutoff(n).tracklength * frame_time;
        
    end
end

% disp('halfway there');

%% Added by Josh Ferguson, Kural Group, The Ohio State University, Sep 21 2016. Modified by Matt Akamatsu

if rejectErraticTracks

    % record number of original tracks

    lengthAlltracksIncludingErratic = length(alltracksNoCutoff);
    
    disp('Rejecting erratic tracks...')
    
    % remove "erratic" tracks (ones that don't follow a linear progression
    % of intensity change over a, say, 4 frame window)
    
    alltracksNoCutoff = trace_rejection_function(alltracksNoCutoff, lifetimesCutoffFrames);
    lengthAllTracks = length(alltracksNoCutoff);
    
    % display how many tracks you rejected.
    
%     fprintf('%d erratic tracks rejected. \n', lengthAlltracksIncludingErratic - lengthAllTracks)
    
end



if cutLifetimes

    nn = 1;

    for n = 1:length(alltracksNoCutoff)



        % Remove tracks with lifetimes < cutoff value listed at top of program

        if alltracksNoCutoff(n).lifetime >= lifetimesCutoffSeconds

            trackLengths(nn) = alltracksNoCutoff(n).tracklength;
            
            alltracks(nn).xs = alltracksNoCutoff(n).xs;
            alltracks(nn).ys = alltracksNoCutoff(n).ys;
            alltracks(nn).Is = alltracksNoCutoff(n).Is;
            alltracks(nn).Istd = alltracksNoCutoff(n).Istd;
            alltracks(nn).trackduration = alltracksNoCutoff(n).trackduration;
            alltracks(nn).tracklength = alltracksNoCutoff(n).tracklength;
            alltracks(nn).lifetime = alltracksNoCutoff(n).lifetime;

            nn = nn + 1;

        end

    end

else
    
    alltracks = alltracksNoCutoff;
    
    nn = 1;

    for n = 1:length(alltracksNoCutoff)
    
        trackLengths(nn) = alltracksNoCutoff(n).tracklength;
        nn = nn+1;
    end
    
    
end

nbTracks = length(alltracks);

% If nbTracksToPlot is -1, then nbTracksToPlot becomes same value as nbTracks
% And if nbTracksToPlot is > nbTracks then set it to nbTracks
if nbTracksToPlot == -1
    
   nbTracksToPlot = nbTracks;
   
elseif nbTracksToPlot > nbTracks
    
    nbTracksToPlot = nbTracks;
   
end
   


LifetimesNoCutoff = (trackLengthsNoCutoff*frame_time)';
Lifetimes         = (trackLengths*frame_time)';




%% Make a csv file for Associate_tracks program
% Format is parent, time, x, y

assocarray = [];

for n = 1:nbTracks
    
    nLarge = n+10^9; % Format Imaris gives for track number.
    
    %     curassocfile = [repmat(nLarge, trackLength(n), 1) alltracks(n).trackduration alltracks(n).xs alltracks(n).ys];
    
    % Edited 12/21/15 by Julian to convert pixels to nm for input into
    % 'Associate_tracks.m'
    curassocarray = [alltracks(n).xs*pixel_size, alltracks(n).ys*pixel_size, NaN(trackLengths(n), 1), NaN(trackLengths(n), 1), NaN(trackLengths(n), 1), alltracks(n).trackduration, repmat(nLarge, trackLengths(n), 1), NaN(trackLengths(n), 1)];
    
    % Edited 12/21/15 by Julian to avoid adding NaNs from tracking to the
    % output file which breaks 'Associate_tracks.m'
    
    % MA edit further to remove just tracks with NaN in first or last
    % element
    
  %  if ~isnan(alltracks(n).xs)
     if isnan(alltracks(n).xs(1))||isnan(alltracks(n).xs(end))
         
         nanTrackNumbers(n) = 1; % list of tracks that start or end with NaN (probably due to track splitting)
        % nanTracks{n} = curassocarray;
         
     else
         
        assocarray = [assocarray; curassocarray];
     end
    
end


% convert to table

assoctable = array2table(assocarray, 'VariableNames', {'PositionX', 'PositionY', 'Unit', 'Category', 'Collection', 'Time', 'TrackID', 'ID'});

% Rearrange by time not track

assoctablesort = sortrows(assoctable, 'Time');

% Save as csv file

assocfilename = fullfile(utrack_folder_name, [curfilename '_utrack_Position.csv']);

%writetable(assoctablesort, 'utrack_Positions.csv')
writetable(assoctablesort, assocfilename);

% Save a csv file with raw lifetimes of tracks
lifetimeTable = array2table(Lifetimes);
lifetimeFilename = fullfile(utrack_folder_name, [curfilename '_utrack_Lifetimes.csv']);
writetable(lifetimeTable, lifetimeFilename);

meanLifetime = mean(Lifetimes);
stdLifetime = std(Lifetimes);

%% Plot track durations

trackdurations = {};

for n = 1:nbTracks

    trackdurations{n} = alltracks(n).trackduration;
    
end

if makefigs

trackPlottingOrder = randperm(nbTracks, nbTracksToPlot);

    
    figure(1); clf;
   
    
    for n = trackPlottingOrder
        
        xs = alltracks(n).xs;
        ys = alltracks(n).ys;
        
        xsNoNan = xs(~isnan(xs));
        ysNoNan = ys(~isnan(ys));
        
        plot(xsNoNan,-ysNoNan), hold all;
    end
    
    title(sprintf('%2.f tracks', nbTracksToPlot));
    xlabel('x');
    ylabel('y');

    
    figure(4); clf;
    
    colormap('jet');

    for n = trackPlottingOrder

%      plot(tracksall{n}, n:length(tracksall{n}):n, '.'); hold all;
      scatter(trackdurations{n}*frame_time, repmat(n,length(trackdurations{n}),1)', [], repmat(Lifetimes(n),length(trackdurations{n}),1),'.'); hold all;

    end
    
    title(sprintf('Track durations for %2.f tracks', nbTracksToPlot));
    xlabel('Time (seconds)');
    ylabel('Track number');
 
    set(gca, 'box', 'on');
    h = colorbar;
    h.Label.String = 'Lifetime(s)';
    
    figure(2); clf;
    
    
    for i = 1:20
        
        subplot(5,4,i);
        
        figure(2);
        errorbar(alltracks(i).Is, alltracks(i).Istd); hold on;
        
        xlabel('Time (frames)');
        ylabel('Intensity');
        
    end
    
    title('Intensity vs frames for 20  tracks');

    
    figure(3); clf;
    histogram(Lifetimes);
    xlabel('Lifetime (seconds)');
    title(sprintf('Lifetimes \n %2.f � %2.f s', meanLifetime, stdLifetime));

    
      
    %% Plot cumulative distribution of lifetimes
    
    [yy, xx] = ecdf(Lifetimes);
    
    t50mask = yy<0.55&yy>0.45;
    t50 = mean(xx(t50mask));
    
    figure(5); clf;
    ecdf(Lifetimes);

    title(sprintf('Track lifetimes cumulative distribution \n t50 = %2.f � %2.f s', t50, stdLifetime));
    xlabel('Lifetime (s)');
    ylabel('Proportion of population');
    
end

%% Save all variables
% this is off for now. comment out to save all the data.
%save(fullfile(utrack_folder_name, curfilename))

fprintf('saved %s \n in %s \n', [curfilename '_utrack_Position.csv'], utrack_folder_name) 

%% Save graphs

% if saveallgraphs
% 
%     saveCurGraphsInLocation([curfilename '_utrack'], utrack_folder_name, 1:5); 
%     
% end

