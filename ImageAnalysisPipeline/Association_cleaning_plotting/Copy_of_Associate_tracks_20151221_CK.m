% Read in imaris result file exported as csv.
% Uses only "Position" results.
% Determine associated reference and corresponding tracks.
%
%
%      Sun Hae Hong
%      David Drubin lab
%      University of  California, Berkeley
%
%      Copyright 2014
%
% Log
% 2014/12/18
% Local BG correction

% 4/30/2015
% Changed
% from '[cor_file, cor_folder_name] = uigetfile('*.Position.csv','Pick corresponding position file for plotting.',previous_file);'
% to '[cor_file, cor_folder_name] = uigetfile('*.csv','Pick corresponding position file for plotting.',previous_file);'
% 4/30/2015
% updated 'function_plot_ref_cor_tracks' to add two more input variables.

% 6/8/2016
% This program now does the following:
% Allow tracks to end at the end of the movie and categorize
% separately
%  Categorize tracks that start at the beginning separately if
% lifetime > certain amount
% To Do: Make input variable that determines whether to throw out or not
% To Do: Test whether results would have changed and keep an index
% To Do: If all works, Call "throw tracks at edges" after track association step


function [ref_assoc_persistant_nb, cor_assoc_persistant_nb, complete_associated_nb] = Associate_tracks_20151221(data)
%% Parameter to use
% align_grid is the size of the grid to be used for alignment. An image
% with a pixel size of align_grid is made.
% number of extra time points before and after the track.

%set a default value for data, then just use all input to get the data so
%this can still be run standalone
if (nargin == 0)
    data.mode = 'manual';
else
    data.mode = 'sequence';
end
num_thrown_tsrrta=0; %tsrrta means "track stat ref . ref tracks assoc"
num_thrown_tsrrpt=0; %tsrrpt means "track stat ref. ref primary track"
num_thrown_roiart=0; %ROI art means "region of interest assoc ref tracks"
num_thrown_ref=0;
num_thrown_cor=0;
times_flagged_ref=0;
times_flagged_cor=0;
extra_tp = 5;
threshold = 1/350^2; % distance btw two spots are less than 350nm.
montage_size = 3;
peak_region_mask = fspecial('disk', montage_size) > 0.01;
bg_size = 5;
bg_region_mask = (fspecial('disk', bg_size) > 0.001);
peak_index_within_bg = (bg_size-montage_size+1):(bg_size+montage_size+1);
bg_region_mask(peak_index_within_bg, peak_index_within_bg)=bg_region_mask(peak_index_within_bg, peak_index_within_bg).*(1-peak_region_mask);

%% info for debugging.
% pixel_size = 108;
%
% analysis_info.ref_channel = 'RFP';
% analysis_info.cor_channel = 'GFP';
%
% ref_folder_name = '/Users/SH/Documents/Matlab/Machine_analysis/Training_Set/AP_R_CL_G_2013-10-25/Field1_1a';
% ref_file = 'Field1_1a_RFP.Position.csv';
% cor_folder_name = '/Users/SH/Documents/Matlab/Machine_analysis/Training_Set/AP_R_CL_G_2013-10-25/Field1_1a';
% cor_file = 'Field1_1a_GFP.Position.csv';
%
% ref_img_folder_name = '/Users/SH/Documents/Matlab/Machine_analysis/Training_Set/AP_R_CL_G_2013-10-25/Field1_1a';
% ref_img_file = 'Field1_1a_RFP.tif';
% cor_img_folder_name =  '/Users/SH/Documents/Matlab/Machine_analysis/Training_Set/AP_R_CL_G_2013-10-25/Field1_1a';
% cor_img_file = 'Field1_1a_GFP.tif';
%
% shift_x = 2;
% shift_y = 2;
% shift_y_nm = 162.9878;
% shift_x_nm = 162.9878;

% 9/9/2015
% Now the code can read in position file even if the file format is not
% part of the name.

%% If on manual mode, allow input
if strcmp(data.mode,'manual')
    % Let users input pixel size
    pixel_size = input('What was the pixel size of the camera in nm? ');
    %pixel_size = 108;
    %% Let user put in channel info
    channel_info = input('Is reference tracks RFP? 1 for yes, 0 for no.');
    while channel_info ~= 1  && channel_info ~= 0
        channel_info = input('The answer needs to be either 0 or 1.\n 1 is for yes and 0 is for no.');
    end
    if channel_info
        analysis_info.ref_channel = 'RFP';
        analysis_info.cor_channel = 'GFP';
        display('You chose RFP channel to be the reference channel, and GFP to be the corresponding channel.')
    else
        analysis_info.ref_channel = 'GFP';
        analysis_info.cor_channel = 'RFP';
        display('You chose GFP channel to be the reference channel, and RFP to be the corresponding channel.')
    end
    %% choose reference file
    load('previous_file.mat')
    [ref_file, ref_folder_name] = uigetfile('*.csv','Pick reference position file for plotting.',previous_file);
    if isdir(ref_folder_name)
        previous_file = fullfile(ref_folder_name, ref_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
    
    %% choose corresponding file
    load('previous_file.mat')
    [cor_file, cor_folder_name] = uigetfile('*.csv','Pick corresponding position file for plotting.',previous_file);
    if isdir(cor_folder_name)
        previous_file = fullfile(cor_folder_name, cor_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
    
    %% choose reference images
    load('previous_file.mat')
    [ref_img_file, ref_img_folder_name] = uigetfile('*.tif','Pick reference images plotting.',previous_file);
    if isdir(ref_img_folder_name)
        previous_file = fullfile(ref_img_folder_name, ref_img_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
    
    %% choose corresponding images
    load('previous_file.mat')
    [cor_img_file, cor_img_folder_name] = uigetfile('*.tif','Pick corresponding images plotting.',previous_file);
    if isdir(cor_img_folder_name)
        previous_file = fullfile(cor_img_folder_name, cor_img_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
    data.ref_img = fillfile(ref_img_folder_name, ref_img_file);
    data.cor_img = fullfile(cor_img_folder_name, cor_img_file);
end

%% if on sequence mode, then get all of the data from the struct
if (data.mode == 'sequence')
    pixel_size = ceil(data.pixelSize * 10e8 / data.M)
    
    %% import choice of channel info
    channel_info = data.rfpRefChannel;
    
%     channel_info = input('Is reference tracks RFP? 1 for yes, 0 for no.');
    while channel_info ~= 1  && channel_info ~= 0
        channel_info = input(' Is reference tracks RFP? 1 for yes, 0 for no. \n The answer needs to be either 0 or 1.\n 1 is for yes and 0 is for no.');
    end
    if channel_info
        analysis_info.ref_channel = 'RFP';
        analysis_info.cor_channel = 'GFP';       
%         display('You chose RFP channel to be the reference channel, and GFP to be the corresponding channel.')
    else
        analysis_info.ref_channel = 'GFP';
        analysis_info.cor_channel = 'RFP';
%         display('You chose GFP channel to be the reference channel, and RFP to be the corresponding channel.')
    end
    home = fullfile(data.source(1:end-4));
    ref_file = strcat(data.curfilename, '_utrack_Position.csv');
    ref_folder_name = fullfile(home, analysis_info.ref_channel, 'Tracking/');
    cor_file = strcat(data.curfilename, '_utrack_Position.csv');
    cor_folder_name = fullfile(home, analysis_info.cor_channel, 'Tracking/');
    temp = dir(strcat(home, analysis_info.ref_channel, '/*.tif'));
    ref_img_file = temp.name;
    ref_img_folder_name = fullfile(home, strcat(analysis_info.ref_channel, '/'));
    temp = dir(strcat(home, analysis_info.cor_channel, '/*.tif'));
    cor_img_file = temp.name;
    cor_img_folder_name = fullfile(home, strcat(analysis_info.cor_channel, '/'));
    
end

%% read in the csv file for ref.
% the columns are for Position X, Position Y, Unit, Category, Collection, Time, Parent (Track num), ID
text_scan_format = '%f%f%s%s%s%f%f%f';
csv_file_ref = fullfile(ref_folder_name, ref_file);
fid = fopen(csv_file_ref);

% Edited 12/21/15 by Julian to conform with output from
% 'extractXYcoordinates.m'
if strcmpi(ref_file(end-19:end), '_utrack_Position.csv')
    headerlines = 1;
elseif strcmpi(ref_file(end-12:end), '_Position.csv')
    headerlines = 4;
elseif strcmpi(ref_file(end-12:end), '.Position.csv')
    headerlines = 2;
else
    display('The file name should end with .Position.csv or _Position.csv.')
    return
end
ref_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',headerlines);
if numel(ref_data{8}) == 0
    text_scan_format = '%f%f%f%s%s%s%f%f%f';
    ref_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',headerlines);
end
fclose(fid);

%% copy data in different format.
% position_data: [Parent, Time, Position X, Position Y].
if numel(ref_data) == 8
    position_data_ref = [ref_data{1,7} - 10^9+1, ref_data{1,6}, ref_data{1,1}, ref_data{1,2}];
elseif numel(ref_data) == 9
    position_data_ref = [ref_data{1,8} - 10^9+1, ref_data{1,7}, ref_data{1,1}, ref_data{1,2}];
else
    error('Dimension mismatch.')
end
clear ref_data


%% read in the csv file for cor.
% the columns are for Position X, Position Y, Unit, Category, Collection, Time, Parent (Track num), ID
csv_file_cor= fullfile(cor_folder_name, cor_file);
fid = fopen(csv_file_cor);
cor_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',headerlines);
fclose(fid);

%% copy data in different format.
% position_data: [Parent, Time, Position X, Position Y].
if numel(cor_data) == 8
    position_data_cor= [cor_data{1,7} - 10^9+1, cor_data{1,6}, cor_data{1,1}, cor_data{1,2}];
elseif numel(cor_data) == 9
    position_data_cor= [cor_data{1,8} - 10^9+1, cor_data{1,7}, cor_data{1,1}, cor_data{1,2}];
else
    error('Dimension mismatch.')
end
clear cor_data


%% Plot all tracks together.
% h = function_plot_ref_cor_tracks(position_data_ref, position_data_cor, 'w', false);
% figure(h); title('Original track positions.')

%% Throw away tracks that touch the edges.
position_data_ref_original = position_data_ref;
position_data_cor_original = position_data_cor;

position_data_cor= throw_tracks_at_edges(position_data_cor, (bg_size+montage_size)*pixel_size);
position_data_ref = throw_tracks_at_edges(position_data_ref, (bg_size+montage_size)*pixel_size);

%% Make lists of tracks
position_data_ref = sortrows(position_data_ref, [1,2]);
track_list_ref = sort(unique(position_data_ref(:,1)));

tp_max = max(max(position_data_ref(:, 2)));

%% Make lists of tracks
position_data_cor= sortrows(position_data_cor, [1,2]);
track_list_cor= sort(unique(position_data_cor(:,1)));

%% Align two channel manually.
shift_y_nm = 0;
shift_x_nm = 0;

% Now we assume that the two channels are aligned in XY. comment this next
% line out if you want to pop up a flag to align the channels.

align_flag = 0;
%align_flag = input('Do you want to align the tracks manually?\n Type 1 for yes and 0 for no.');
while align_flag ~= 1  && align_flag ~= 0
    align_flag = input('The answer needs to be either 0 or 1.\n 1 is for yes and 0 is for no.');
end
if align_flag
    display('Click a blue spot first and then click the corresponding red spot.');
    figure(h)
    [x1,y1] = ginput(1);  % Select a point with the mouse
    text(x1,y1,['(', num2str(x1), ', ', num2str(y1), ')']) % show the clicked coordinate
    display(['First locus (cor):',  num2str(x1), ', ', num2str(y1)]);
    [x2,y2] = ginput(1);  % Select a point with the mouse
    text(x2,y2,['(', num2str(x2), ', ', num2str(y2), ')']) % show the clicked coordinate
    display(['Second locus (ref):',  num2str(x2), ', ', num2str(y2)]);
    display(['First locus - second locus =',  num2str(x1 - x2), ', ', num2str(y1 - y2)]);
    %% Align two channel manually.
    %% Let users input shift in x and y until satisfied.
    position_ref_aligned = position_data_ref;
    position_cor_aligned = position_data_cor;
    %     display('Based on the figure, "Original track positions", please give shift in x and y.');
    
    satisfied = false;
    while ~satisfied
        close(h);
        %% ask users for shift
        
        %         shift_y_nm = input('How much shift toward upper direction to reference tracks do you want to apply?');
        %         shift_x_nm = input('How much shift towards right to reference tracks do you want to apply?');
        
        shift_y_nm = shift_y_nm + y1-y2;
        shift_x_nm = shift_x_nm + x1-x2;
        %%
        position_ref_aligned(:,3) = position_data_ref(:,3) + shift_x_nm;
        position_ref_aligned(:,4) = position_data_ref(:,4) + shift_y_nm;
        
        %% show changed tracks
        h = function_plot_ref_cor_tracks(position_ref_aligned, position_data_cor);
        figure(h); title('Manual shift is applied.');
        
        satisfied = input('Are you satisfied with your alignment?\n Type 1 for yes and 0 for no.');
        while satisfied ~= 1  && satisfied ~= 0
            satisfied = input('The answer needs to be either 0 or 1.\n 1 is for yes and 0 is for no.');
        end
        if ~satisfied
            [x1,y1] = ginput(1);  % Select a point with the mouse
            text(x1,y1,['(', num2str(x1), ', ', num2str(y1), ')']) % show the clicked coordinate
            display(['First locus:',  num2str(x1), ', ', num2str(y1)]);
            [x2,y2] = ginput(1);  % Select a point with the mouse
            text(x2,y2,['(', num2str(x2), ', ', num2str(y2), ')']) % show the clicked coordinate
            display(['Second locus:',  num2str(x2), ', ', num2str(y2)]);
            display(['First locus - second locus =',  num2str(x1 - x2), ', ', num2str(y1 - y2)]);
        end
    end
end

position_ref_aligned = position_data_ref;
position_cor_aligned = position_data_cor;

position_ref_aligned(:,3) = position_data_ref(:,3) + shift_x_nm;
position_ref_aligned(:,4) = position_data_ref(:,4) + shift_y_nm;

%% convert shifts in nm to pixels
shift_y = round(shift_y_nm / pixel_size);
shift_x = round(shift_x_nm / pixel_size);

%% make copies of tracks that start at beginning of movie (low) or end at end of movie (high)
[ROI_track_ref, ROI_track_ref_low, ROI_track_ref_high, ...
    track_stat_ref, track_stat_ref_low, track_stat_ref_high] = ...
    throw_tracks_at_timepoint(track_list_ref, position_data_ref, extra_tp);

[ROI_track_cor, ROI_track_cor_low, ROI_track_cor_high, ...
    track_stat_cor, track_stat_cor_low, track_stat_cor_high] = ...
    throw_tracks_at_timepoint(track_list_cor, position_data_cor, extra_tp);


%% Read in images and manually go through them.
cor_stack_3D = double(openTiffStack(fullfile(cor_img_folder_name, cor_img_file)));
ref_stack_3D = double(openTiffStack(fullfile(ref_img_folder_name, ref_img_file)));
[im_size1, im_size2, im_size3] = size(ref_stack_3D);

%% Determine background intensity of each time point as median
cor_stack_bg = NaN*ones(1, im_size3);
ref_stack_bg = NaN*ones(1, im_size3);

for tp = 1:im_size3
    cor_stack_bg(tp) = median(reshape(cor_stack_3D(:,:,tp), [1, im_size1*im_size2]));
    ref_stack_bg(tp) = median(reshape(ref_stack_3D(:,:,tp), [1, im_size1*im_size2]));
end
%% Determine Associated Tracks
display('Associating corresponding and reference tracks.')
ref_tracks_taken = [];
ref_primary_track = [];
cor_tracks_taken = [];

ref_matrix_size2 = size(track_stat_ref.xcoord, 2);
cor_matrix_size2 = size(track_stat_cor.xcoord, 2);

% ref_track_info = [ref track id, number of reference tracks, number of corresponding trcaks].
ref_track_info = zeros*ones(numel(track_stat_ref.track_list), 3);
count = 0;
for ref_track = track_stat_ref.track_list'
    if sum(ref_track == ref_tracks_taken)
        continue
    end
    ref_primary_track = [ref_primary_track; ref_track];
    ref_tracks_taken = [ref_tracks_taken; ref_track];
    
    dist_score_mat = (track_stat_cor.xcoord- repmat(track_stat_ref.xcoord(:, ref_track), [1, cor_matrix_size2]) - shift_x_nm).^2 +...
        (track_stat_cor.ycoord- repmat(track_stat_ref.ycoord(:, ref_track), [1, cor_matrix_size2]) - shift_y_nm).^2;
    dist_score_mat = 1./dist_score_mat;
    dist_score = nansum(dist_score_mat,1);
    overlapping_time = max(10, sum(~isnan(dist_score_mat),1));
    candidates = find(dist_score > threshold*overlapping_time/2);
    
    ROI_x = [round(min(track_stat_ref.xcoord(:,ref_track))/pixel_size - 5), round(max(track_stat_ref.xcoord(:,ref_track))/pixel_size + 5)];
    ROI_y = [round(min(track_stat_ref.ycoord(:,ref_track))/pixel_size - 5), round(max(track_stat_ref.ycoord(:,ref_track))/pixel_size + 5)];
    
    tp_prev = ROI_track_ref(ref_track, 1).tp_prev;
    tp_post = ROI_track_ref(ref_track, 1).tp_post;
    
    ROI_track_ref(ref_track, 1).associated_cor_track = [];
    ROI_track_ref(ref_track, 1).associated_ref_track = [ref_track];
    
    ROI_track_ref(ref_track, 1).ref_tp_first = ROI_track_ref(ref_track, 1).tp_first;
    ROI_track_ref(ref_track, 1).ref_tp_last = ROI_track_ref(ref_track, 1).tp_last;
    
    ROI_track_ref(ref_track, 1).cor_tp_first = Inf;
    ROI_track_ref(ref_track, 1).cor_tp_last = -Inf;
    
    %%
    for cor_track = candidates;
        if sum(cor_track == cor_tracks_taken)
            continue
        end
        ROI_track_ref(ref_track, 1).associated_cor_track = [ROI_track_ref(ref_track, 1).associated_cor_track, cor_track];
        
        tp_prev = min(tp_prev, ROI_track_cor(cor_track, 1).tp_prev);
        tp_post = max(tp_post, ROI_track_cor(cor_track, 1).tp_post);
        
        ROI_track_ref(ref_track, 1).cor_tp_first = min(ROI_track_ref(ref_track, 1).cor_tp_first, ROI_track_cor(cor_track, 1).tp_first);
        ROI_track_ref(ref_track, 1).cor_tp_last = max(ROI_track_ref(ref_track, 1).cor_tp_last, ROI_track_cor(cor_track, 1).tp_last);
        
        cor_tracks_taken = [cor_tracks_taken, cor_track];
        
        %% search for reference tracks associated with corresponding tracks.
        dist_score_mat = (track_stat_ref.xcoord +shift_y_nm - repmat(track_stat_cor.xcoord(:, cor_track), [1, ref_matrix_size2])).^2 +...
            (track_stat_ref.ycoord+shift_x_nm - repmat(track_stat_cor.ycoord(:, cor_track), [1, ref_matrix_size2])).^2;
        dist_score_mat = 1./dist_score_mat;
        dist_score = nansum(dist_score_mat,1);
        
        overlapping_time = max(10, sum(~isnan(dist_score_mat),1));
        ref_candidates = find(dist_score > threshold*overlapping_time/2);
        
        n_associated_ref = numel(ref_candidates);
        for ref_track2 = ref_candidates
            if sum(ref_track2 == ref_tracks_taken)
                continue
            end
            % if the time of reference peak overlaps, it is not used.
            if intersect(ROI_track_ref(ref_track2, 1).tp_first:ROI_track_ref(ref_track2, 1).tp_last,...
                    ROI_track_ref(ref_track, 1).tp_first:ROI_track_ref(ref_track, 1).tp_last)
                continue
            end
            
            ROI_track_ref(ref_track, 1).associated_ref_track = [ROI_track_ref(ref_track, 1).associated_ref_track, ref_track2];
            ref_tracks_taken = [ref_tracks_taken; ref_track2];
            tp_prev = min(tp_prev, ROI_track_ref(ref_track2, 1).tp_prev);
            tp_post = max(tp_post, ROI_track_ref(ref_track2, 1).tp_post);
            
            ROI_track_ref(ref_track, 1).ref_tp_first = min(ROI_track_ref(ref_track, 1).ref_tp_first, ROI_track_ref(ref_track2, 1).tp_first);
            ROI_track_ref(ref_track, 1).ref_tp_last = max(ROI_track_ref(ref_track, 1).ref_tp_last, ROI_track_ref(ref_track2, 1).tp_last);
            
        end
    end
    
    %% if corresponding tracks overlap in time, the track is not analyzed
    tp_sequence = zeros(1,tp_max);
    for cor_track = ROI_track_ref(ref_track, 1).associated_cor_track
        tp_sequence(ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last) = ...
            tp_sequence(ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last) + 1;
    end
    if any(tp_sequence > 1)
        ref_primary_track = ref_primary_track(1:end-1);
        ref_tracks_taken = [ref_tracks_taken; ref_track];
        continue
    end
    
    
    ROI_track_ref(ref_track,1).tp_prev_all = tp_prev;
    ROI_track_ref(ref_track,1).tp_post_all = tp_post;
    
    count = count + 1;
    ref_track_info(count, 1) = ref_track;
    ref_track_info(count, 2) = numel(ROI_track_ref(ref_track, 1).associated_ref_track);
    ref_track_info(count, 3) = numel(ROI_track_ref(ref_track, 1).associated_cor_track);
end

%% create list of associated and unassociated tracks

track_stat_ref.ref_primary_track = ref_primary_track;
track_stat_ref.ref_tracks_associated = ref_track_info(and(ref_track_info(:,2)>0, ref_track_info(:,3)>0), 1);
track_stat_ref.ref_tracks_unassociated = ref_track_info(and(ref_track_info(:,2)>0, ref_track_info(:,3)==0), 1);
track_stat_cor.cor_tracks_unassociated = setdiff(track_stat_cor.track_list, cor_tracks_taken);

clear ref_track_info
%% Determine inter- and extra-polated track
for ref_track = ref_primary_track'
    
    tp_prev_all = ROI_track_ref(ref_track,1).tp_prev_all;
    tp_post_all = ROI_track_ref(ref_track,1).tp_post_all;
    
    % initialize extrapolated coord.
    ROI_track_ref(ref_track).ref_coord_extrapol = NaN*ones(tp_post_all - tp_prev_all + 1, 2);
    ROI_track_ref(ref_track).cor_coord_extrapol = NaN*ones(tp_post_all - tp_prev_all + 1, 2);
    
    %% reference coordinates
    coord_extrapol = ROI_track_ref(ref_track).ref_coord_extrapol;
    for ref_track2 = ROI_track_ref(ref_track, 1).associated_ref_track
        coord_extrapol((ROI_track_ref(ref_track2, 1).tp_first - tp_prev_all + 1):(ROI_track_ref(ref_track2, 1).tp_last - tp_prev_all + 1), :) = ...
            ROI_track_ref(ref_track2, 1).center_coord;
    end
    % fill in the coordinates before, after and in between
    filled_index = find(~isnan(coord_extrapol(:,1)));
    if filled_index(1) > 1
        coord_extrapol(1:(filled_index(1)-1),1) = coord_extrapol(filled_index(1),1);
        coord_extrapol(1:(filled_index(1)-1),2) = coord_extrapol(filled_index(1),2);
    end
    if filled_index(end) < numel(coord_extrapol(:,1))
        coord_extrapol((filled_index(end)+1):end,1) = coord_extrapol(filled_index(end),1);
        coord_extrapol((filled_index(end)+1):end,2) = coord_extrapol(filled_index(end),2);
    end
    % Fill in the gaps by averaging the positions before and after the gap
    coord_extrapol = fill_gap(coord_extrapol);
    ROI_track_ref(ref_track).ref_coord_extrapol = coord_extrapol;
    
    %% corresponding coordinates
    coord_extrapol = ROI_track_ref(ref_track).cor_coord_extrapol;
    if numel(ROI_track_ref(ref_track, 1).associated_cor_track) == 0;
        ROI_track_ref(ref_track).cor_coord_extrapol = NaN*ones(size(ROI_track_ref(ref_track).ref_coord_extrapol));
        % if the coordinate is not at the boundary assign the known
        % coordinate to the unknown coordinate
        if min(ROI_track_ref(ref_track).ref_coord_extrapol(:,1) + shift_x_nm) >= (1+bg_size)*pixel_size &&...
                max(ROI_track_ref(ref_track).ref_coord_extrapol(:,1) + shift_x_nm) <= (track_stat_ref.maxx - bg_size*pixel_size) &&...
                min(ROI_track_ref(ref_track).ref_coord_extrapol(:,2) + shift_y_nm) >= (1+bg_size)*pixel_size &&...
                max(ROI_track_ref(ref_track).ref_coord_extrapol(:,2) + shift_y_nm) <= (track_stat_ref.maxy - bg_size*pixel_size)
            ROI_track_ref(ref_track).cor_coord_extrapol(:,1) = ROI_track_ref(ref_track).ref_coord_extrapol(:,1) + shift_x_nm;
            ROI_track_ref(ref_track).cor_coord_extrapol(:,2) = ROI_track_ref(ref_track).ref_coord_extrapol(:,2) + shift_y_nm;
        end
    else
        for cor_track = ROI_track_ref(ref_track, 1).associated_cor_track
            coord_extrapol((ROI_track_cor(cor_track, 1).tp_first - tp_prev_all + 1):(ROI_track_cor(cor_track, 1).tp_last - tp_prev_all + 1), :) = ...
                ROI_track_cor(cor_track, 1).center_coord;
        end
        % fill in the coordinates before, after and in between
        filled_index = find(~isnan(coord_extrapol(:,1)));
        if filled_index(1) > 1
            coord_extrapol(1:(filled_index(1)-1),1) = coord_extrapol(filled_index(1),1);
            coord_extrapol(1:(filled_index(1)-1),2) = coord_extrapol(filled_index(1),2);
        end
        if filled_index(end) < numel(coord_extrapol(:,1))
            coord_extrapol((filled_index(end)+1):end,1) = coord_extrapol(filled_index(end),1);
            coord_extrapol((filled_index(end)+1):end,2) = coord_extrapol(filled_index(end),2);
        end
        % Fill in the gaps by averaging the positions before and after the gap
        prev_coord = coord_extrapol;
        coord_extrapol = fill_gap(coord_extrapol);
        ROI_track_ref(ref_track).cor_coord_extrapol = coord_extrapol;
    end
end

%% Based on inter- and extra-polated track, determine region of interest (ROI)
for ref_track = ref_primary_track'
    % ROI to include all of ref_coord_extrapol
    ROI_x(1,1) = round(min(ROI_track_ref(ref_track).ref_coord_extrapol(:,1)/pixel_size)) - 5;
    ROI_x(1,2) = round(max(ROI_track_ref(ref_track).ref_coord_extrapol(:,1)/pixel_size)) + 5;
    ROI_y(1,1) = round(min(ROI_track_ref(ref_track).ref_coord_extrapol(:,2)/pixel_size)) - 5;
    ROI_y(1,2) = round(max(ROI_track_ref(ref_track).ref_coord_extrapol(:,2)/pixel_size)) + 5;
    
    % ROI to include all of cor_coord_extrapol
    if ~any(isnan(ROI_track_ref(ref_track).cor_coord_extrapol))
        ROI_x(1,1) = min(ROI_x(1,1), round(min(ROI_track_ref(ref_track).cor_coord_extrapol(:,1)/pixel_size) - shift_x) - 5);
        ROI_x(1,2) = max(ROI_x(1,2), round(max(ROI_track_ref(ref_track).cor_coord_extrapol(:,1)/pixel_size) - shift_x) + 5);
        ROI_y(1,1) = min(ROI_y(1,1), round(min(ROI_track_ref(ref_track).cor_coord_extrapol(:,2)/pixel_size) - shift_y) - 5);
        ROI_y(1,2) = max(ROI_y(1,2), round(min(ROI_track_ref(ref_track).cor_coord_extrapol(:,2)/pixel_size) - shift_y) - 5);
    end
    ROI_track_ref(ref_track,1).ROI_x(1,1) = max([ROI_x(1,1),1, 1 - shift_x]);
    ROI_track_ref(ref_track,1).ROI_x(1,2) = min([ROI_x(1,2), im_size2, im_size2 - shift_x]);
    ROI_track_ref(ref_track,1).ROI_y(1,1) = max([ROI_y(1,1),1, 1 - shift_y]);
    ROI_track_ref(ref_track,1).ROI_y(1,2) = min([ROI_y(1,2), im_size1, im_size1 - shift_y]);
end
%% Initialize reference montage
display('saving montage.')
associated_tracks = struct('ref_movie', {},'cor_movie', {});
for ref_track = ref_primary_track'
    % for the time frame of interest, determine the integrated intensity of the
    % region of interest.
    ROI_x = ROI_track_ref(ref_track,1).ROI_x;
    ROI_y = ROI_track_ref(ref_track,1).ROI_y;
    
    ROI_delx = ROI_x(1,2) - ROI_x(1,1) + 1;
    ROI_dely = ROI_y(1,2) - ROI_y(1,1) + 1;
    
    tp_prev = ROI_track_ref(ref_track,1).tp_prev_all;
    tp_post = ROI_track_ref(ref_track,1).tp_post_all;
    
    %     associated_tracks(ref_track, 1).ref_movie = NaN*ones(ROI_dely, ROI_delx*(tp_post - tp_prev+ 1));
    %     associated_tracks(ref_track, 1).cor_movie = NaN*ones(ROI_dely, ROI_delx*(tp_post - tp_prev+ 1));
    if ~any(isnan(ROI_track_ref(ref_track).ref_coord_extrapol))
        associated_tracks(ref_track, 1).cor_montage = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
        associated_tracks(ref_track, 1).ref_montage = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
        associated_tracks(ref_track, 1).cor_intensity = NaN*ones(2, (tp_post - tp_prev+ 1));
        associated_tracks(ref_track, 1).ref_intensity = NaN*ones(2, (tp_post - tp_prev+ 1));
        associated_tracks(ref_track, 1).cor_intensity_bg_corrected = NaN*ones(2, (tp_post - tp_prev+ 1));
        associated_tracks(ref_track, 1).ref_intensity_bg_corrected = NaN*ones(2, (tp_post - tp_prev+ 1));
        
        associated_tracks(ref_track, 1).cor_montage_bg = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
        associated_tracks(ref_track, 1).ref_montage_bg = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
        
    end
end

for ref_track = ref_primary_track'
    
    ROI_x = ROI_track_ref(ref_track,1).ROI_x;
    ROI_y = ROI_track_ref(ref_track,1).ROI_y;
    
    ROI_delx = ROI_x(1,2) - ROI_x(1,1) + 1;
    ROI_dely = ROI_y(1,2) - ROI_y(1,1) + 1;
    
    tp_prev = ROI_track_ref(ref_track,1).tp_prev_all;
    tp_post = ROI_track_ref(ref_track,1).tp_post_all;
    
    ROI_x_range = ROI_x(1):ROI_x(2);
    ROI_y_range = ROI_y(1):ROI_y(2);
    associated_tracks(ref_track, 1).ref_movie = reshape(ref_stack_3D(ROI_y_range, ROI_x_range, tp_prev:tp_post),...
        [ROI_dely, ROI_delx*(tp_post - tp_prev+ 1), 1]);
    associated_tracks(ref_track, 1).cor_movie = reshape(cor_stack_3D(ROI_y_range + shift_y, ROI_x_range + shift_x,...
        tp_prev:tp_post), [ROI_dely, ROI_delx*(tp_post - tp_prev+ 1), 1]);
    % Determine the integrated intensity
    if ~any(isnan(ROI_track_ref(ref_track).ref_coord_extrapol))
        for tp = ROI_track_ref(ref_track,1).tp_prev_all:ROI_track_ref(ref_track,1).tp_post_all
            rel_tp = (tp - ROI_track_ref(ref_track,1).tp_prev_all + 1);
            
            % save reference montage
            peak_x = round(ROI_track_ref(ref_track).ref_coord_extrapol(rel_tp, 1)/pixel_size);
            peak_y = round(ROI_track_ref(ref_track).ref_coord_extrapol(rel_tp, 2)/pixel_size);
            
            peak_region_image = ref_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp).*peak_region_mask;
            associated_tracks(ref_track, 1).ref_intensity(1, rel_tp) = sum(sum(peak_region_image));
            associated_tracks(ref_track, 1).ref_intensity(2, rel_tp) = std2(peak_region_image);
            associated_tracks(ref_track, 1).ref_montage(:,(rel_tp*7-6):(rel_tp*7)) = peak_region_image;
            
            bg_region_image = ref_stack_3D(peak_y-bg_size:peak_y+bg_size, peak_x-bg_size:peak_x+bg_size, tp).*bg_region_mask;
            bg_value = median(bg_region_image(bg_region_mask));
            peak_region_image_bg_corrected = (ref_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp)...
                - bg_value).*peak_region_mask;
            associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, rel_tp) = sum(sum(peak_region_image_bg_corrected));
            associated_tracks(ref_track, 1).ref_intensity_bg_corrected(2, rel_tp) = std2(peak_region_image_bg_corrected);
            
            associated_tracks(ref_track, 1).ref_montage_bg(:,(rel_tp*7-6):(rel_tp*7)) = peak_region_image_bg_corrected;
            
            % save corresponding montage
            if ~any(isnan(ROI_track_ref(ref_track).cor_coord_extrapol))
                
                peak_x = round(ROI_track_ref(ref_track).cor_coord_extrapol(rel_tp, 1)/pixel_size);
                peak_y = round(ROI_track_ref(ref_track).cor_coord_extrapol(rel_tp, 2)/pixel_size);
                
                peak_region_image = cor_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp).*peak_region_mask;
                associated_tracks(ref_track, 1).cor_intensity(1, rel_tp) = sum(sum(peak_region_image));
                associated_tracks(ref_track, 1).cor_intensity(2, rel_tp) = std2(peak_region_image);
                associated_tracks(ref_track, 1).cor_montage(:,(rel_tp*7-6):(rel_tp*7)) = peak_region_image;
                
                bg_region_image = cor_stack_3D(peak_y-bg_size:peak_y+bg_size, peak_x-bg_size:peak_x+bg_size, tp).*bg_region_mask;
                bg_value = median(bg_region_image(bg_region_mask));
                peak_region_image_bg_corrected = (cor_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp)...
                    - bg_value).*peak_region_mask;
                associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, rel_tp) = sum(sum(peak_region_image_bg_corrected));
                associated_tracks(ref_track, 1).cor_intensity_bg_corrected(2, rel_tp) = std2(peak_region_image_bg_corrected);
                
                associated_tracks(ref_track, 1).cor_montage_bg(:,(rel_tp*7-6):(rel_tp*7)) = peak_region_image_bg_corrected;
            end
        end
    end
    
end

%%
if size(track_stat_cor.cor_tracks_unassociated, 1) > size(track_stat_cor.cor_tracks_unassociated, 2)
    track_stat_cor.cor_tracks_unassociated = track_stat_cor.cor_tracks_unassociated';
end

if ~isempty(track_stat_cor.cor_tracks_unassociated)

for cor_track = track_stat_cor.cor_tracks_unassociated
    
    ROI_x = [round(min(track_stat_cor.xcoord(:,cor_track))/pixel_size - 5), round(max(track_stat_cor.xcoord(:,cor_track))/pixel_size + 5)];
    ROI_y = [round(min(track_stat_cor.ycoord(:,cor_track))/pixel_size - 5), round(max(track_stat_cor.ycoord(:,cor_track))/pixel_size + 5)];
    
    ROI_track_cor(cor_track, 1).ROI_x(1,1) = max([ROI_x(1,1), 1, 1 + shift_x]);
    ROI_track_cor(cor_track, 1).ROI_x(1,2) = min([ROI_x(1,2), im_size2, im_size2 + shift_x]);
    ROI_track_cor(cor_track, 1).ROI_y(1,1) = max([ROI_y(1,1), 1, 1 + shift_y]);
    ROI_track_cor(cor_track, 1).ROI_y(1,2) = min([ROI_y(1,2), im_size1, im_size1 + shift_y]);
    
    
    tp_prev = ROI_track_cor(cor_track, 1).tp_prev;
    tp_post = ROI_track_cor(cor_track, 1).tp_post;
    
    %% Determine inter- and extra-polated track
    % initialize extrapolated coord.
    coord_extrapol = NaN*ones(tp_post - tp_prev + 1, 2);
    coord_extrapol((ROI_track_cor(cor_track, 1).tp_first - tp_prev + 1):(ROI_track_cor(cor_track, 1).tp_last - tp_prev + 1), :) = ...
        ROI_track_cor(cor_track, 1).center_coord;
    
    % fill in the coordinates before, after and in between
    filled_index = find(~isnan(coord_extrapol(:,1)));
    if filled_index(1) > 1
        coord_extrapol(1:(filled_index(1)-1),1) = coord_extrapol(filled_index(1),1);
        coord_extrapol(1:(filled_index(1)-1),2) = coord_extrapol(filled_index(1),2);
    end
    if filled_index(end) < numel(coord_extrapol(:,1))
        coord_extrapol((filled_index(end)+1):end,1) = coord_extrapol(filled_index(end),1);
        coord_extrapol((filled_index(end)+1):end,2) = coord_extrapol(filled_index(end),2);
    end
    % Fill in the gaps by averaging the positions before and after the gap
    coord_extrapol = fill_gap(coord_extrapol);
    ROI_track_cor(cor_track).cor_coord_extrapol = coord_extrapol;
    ROI_track_cor(cor_track).ref_coord_extrapol = NaN*ones(size(ROI_track_cor(cor_track).cor_coord_extrapol));
    % if the coordinate is not at the boundary assign the known
    % coordinate to the unknown coordinate
    if min(ROI_track_cor(cor_track).cor_coord_extrapol(:,1) - shift_x_nm) >= (1+bg_size)*pixel_size &&...
            max(ROI_track_cor(cor_track).cor_coord_extrapol(:,1) - shift_x_nm) <= (track_stat_ref.maxx - bg_size*pixel_size) &&...
            min(ROI_track_cor(cor_track).cor_coord_extrapol(:,2) - shift_y_nm) >= (1+bg_size)*pixel_size &&...
            max(ROI_track_cor(cor_track).cor_coord_extrapol(:,2) - shift_y_nm) <= (track_stat_ref.maxy - bg_size*pixel_size)
        ROI_track_cor(cor_track).ref_coord_extrapol(:,1) = ROI_track_cor(cor_track).cor_coord_extrapol(:,1) - shift_x_nm;
        ROI_track_cor(cor_track).ref_coord_extrapol(:,2) = ROI_track_cor(cor_track).cor_coord_extrapol(:,2) - shift_y_nm;
    end
    
end

end
%% Remove any tracks that went out of bounds in time

% first do this for unassociated ref tracks, and categorize those tracks
% that went out of bounds in 2 new vectors so they can later be analyzed
% independently (high and low violators)
num_thrown=0;
lowindex=0;
highindex=0;

track_stat_ref.ref_tracks_unassociated_high_breach = [];
track_stat_ref.ref_tracks_unassociated_low_breach  = [];

track_stat_cor.cor_tracks_unassociated_high_breach = [];
track_stat_cor.cor_tracks_unassociated_low_breach  = [];

track_stat_ref.ref_tracks_associated_high_breach = [];
track_stat_ref.ref_tracks_associated_low_breach = [];

track_stat_cor.cor_tracks_associated_high_breach  = [];
track_stat_cor.cor_tracks_associated_low_breach  = [];

if ~isempty(track_stat_ref.ref_tracks_unassociated)

    for i=1:size(track_stat_ref.ref_tracks_unassociated)
        aa=track_stat_ref.ref_tracks_unassociated(i-num_thrown);
        
        I=find(track_stat_ref_low.track_list==aa);
        if isempty(I)==1
            I=0;
        end
        
        if I>0
            lowindex=lowindex+1;
            track_stat_ref.ref_tracks_unassociated_low_breach(lowindex)=aa;
        end
        
        J=find(track_stat_ref_high.track_list==aa);
        if isempty(J)==1
            J=0;
        end
        
        if J>0
            highindex=highindex+1;
            track_stat_ref.ref_tracks_unassociated_high_breach(highindex)=aa;
        end
        
        if or(I>0, J>0) % delete appropriate vector entry to keep lists pure
            track_stat_ref.ref_tracks_unassociated(i-num_thrown,:)=[];
            num_thrown=num_thrown+1;
        end
    end
    
% else
%     
%     track_stat_ref.ref_tracks_unassociated_low_breach  = []
%     track_stat_ref.ref_tracks_unassociated_high_breach = []
    
end
% next do this for unassociated cor tracks
% tscctu stands for track stat cor cor tracks unassociated.

undeleted_tscctu = track_stat_cor.cor_tracks_unassociated;

num_thrown=0;
lowindex=0;
highindex=0;

% initialize track_stat_cor.cor_tracks_unassociated_low_breach and
% _high_breach (they'll be empty if none are added)

track_stat_cor.cor_tracks_unassociated_low_breach  = [];
track_stat_cor.cor_tracks_unassociated_high_breach = [];

if ~isempty(track_stat_cor.cor_tracks_unassociated)
    
    for i=1:size(track_stat_cor.cor_tracks_unassociated,2)
        aa=track_stat_cor.cor_tracks_unassociated(i-num_thrown);
        
        I=find(track_stat_cor_low.track_list==aa);
        if isempty(I)==1
            I=0;
        end
        
        if I>0
            lowindex=lowindex+1;
            track_stat_cor.cor_tracks_unassociated_low_breach(lowindex)=aa;
        end
        
        J=find(track_stat_cor_high.track_list==aa);
        if isempty(J)==1
            J=0;
        end
        
        if J>0
            highindex=highindex+1;
            track_stat_cor.cor_tracks_unassociated_high_breach(highindex)=aa;
        end
        
        if or(I>0, J>0) % delete appropriate vector entry to keep lists pure
            track_stat_cor.cor_tracks_unassociated(i-num_thrown)=[];
            num_thrown=num_thrown+1;
        end
    end
    
% else
%     
%     track_stat_cor.cor_tracks_unassociated_low_breach  = []
%     track_stat_cor.cor_tracks_unassociated_high_breach = []
    
end
% next do this for associated ref tracks and associated cor tracks

lowindexa=0;
lowindexb=0;
highindexa=0;
highindexb=0;


for i=1:size(track_stat_ref.ref_tracks_associated)
    
    aa=track_stat_ref.ref_tracks_associated(i-num_thrown_tsrrta);
    
    bb=ROI_track_ref(track_stat_ref.ref_tracks_associated(i-num_thrown_tsrrta)).associated_cor_track;
    
    
    I=find(track_stat_ref_low.track_list==aa);
    if isempty(I)==1
        I=0;
    end
    
    % this seems to populate a list
    % "track_stat_ref.ref_tracks_associated_low_breach" which lists the
    % track numbers that are in this category(associated reference tracks
    % that start at the beginning of the movie)
    
    if I>0
        lowindexa=lowindexa+1;
        track_stat_ref.ref_tracks_associated_low_breach(lowindexa)=aa;
    end
    
    J=find(track_stat_ref_high.track_list==aa);
    if isempty(J)==1
        J=0;
    end
    
     % this seems to populate a list
    % "track_stat_ref.ref_tracks_associated_high_breach" which lists the
    % track numbers that are in this category(associated reference tracks
    % that end at the end of the movie)
    
    
    if J>0
        highindexa=highindexa+1;
        track_stat_ref.ref_tracks_associated_high_breach(highindexa)=aa;
    end
    
    
    K=0;
    L=0;
    
    for j=1:size(bb,2)
        Ktemp=find(track_stat_cor_low.track_list==bb(j));
        if isempty(Ktemp)==1
            Ktemp=0;
        end
        
        if Ktemp>0
            lowindexb=lowindexb+1;
            track_stat_cor.cor_tracks_associated_low_breach(lowindexb)=bb(j);
        end
        
        Ltemp=find(track_stat_cor_high.track_list==bb(j));
        if isempty(Ltemp)==1
            Ltemp=0;
        end
        
        if Ltemp>0
            highindexb=highindexb+1;
            track_stat_cor.cor_tracks_associated_high_breach(highindexb)=bb(j);
        end
        
        K=K+Ktemp;
        L=L+Ltemp;
    end
    
    
    ref_unassoc_persistant_tracks = [];
    cor_unassoc_persistant_tracks = [];
    ref_assoc_persistant_tracks = [];
    cor_assoc_persistant_tracks = [];
    ref_unassoc_persistant_status = ismember(track_stat_ref.ref_tracks_unassociated_low_breach, track_stat_ref.ref_tracks_unassociated_high_breach);
    cor_unassoc_persistant_status = ismember(track_stat_cor.cor_tracks_unassociated_low_breach, track_stat_cor.cor_tracks_unassociated_high_breach);
    ref_assoc_persistant_status = ismember(track_stat_ref.ref_tracks_associated_low_breach, track_stat_ref.ref_tracks_associated_high_breach);
    cor_assoc_persistant_status = ismember(track_stat_cor.cor_tracks_associated_low_breach, track_stat_cor.cor_tracks_associated_high_breach);
    
    for ii = 1:numel(ref_unassoc_persistant_status)
        if ref_unassoc_persistant_status(ii) == 1
            ref_unassoc_persistant_tracks = [ref_unassoc_persistant_tracks ref_unassoc_persistant_status(ii)];
        end
    end
    for ii = 1:numel(cor_unassoc_persistant_status)
        if cor_unassoc_persistant_status(ii) == 1
            cor_unassoc_persistant_tracks = [cor_unassoc_persistant_tracks cor_unassoc_persistant_status(ii)];
        end
    end
    for ii = 1:numel(ref_assoc_persistant_status)
        if ref_assoc_persistant_status(ii) == 1
            ref_assoc_persistant_tracks = [ref_assoc_persistant_tracks ref_assoc_persistant_status(ii)];
        end
    end
    for ii = 1:numel(cor_assoc_persistant_status)
        if cor_assoc_persistant_status(ii) == 1
            cor_assoc_persistant_tracks = [cor_assoc_persistant_tracks cor_assoc_persistant_status(ii)];
        end
    end
    
   % ref_unassoc_persistant_nb = numel(ref_unassoc_persistant_tracks);
   % ref_assoc_persistant_nb = numel(ref_assoc_persistant_tracks);
    %cor_unassoc_persistant_nb = numel(cor_unassoc_persistant_tracks);
    %cor_assoc_persistant_nb = numel(cor_assoc_persistant_tracks);
    
   % ref_unassoc_persistant = numel(ref_unassoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_unassociated);
    %ref_assoc_persistant = numel(ref_assoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_associated);
    %cor_unassoc_persistant = numel(cor_unassoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_unassociated);
    %cor_assoc_persistant = numel(cor_assoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_unassociated);
    
    % change by CK 180802, reference to the wrong variable
    %num_thrown_ref_associated = highindexa + lowindexa - numel(ref_assoc_persistant);
    %num_thrown_cor_associated = highindexb + lowindexb - numel(cor_assoc_persistant);
    
    %num_thrown_ref_associated = highindexa + lowindexa - numel(ref_assoc_persistant_tracks);
   % num_thrown_cor_associated = highindexb + lowindexb - numel(cor_unassoc_persistant_tracks);
    
    % complete_associated_nb is in a loop and variable is overwritten and thus output wrong in CME analysis (so far cannot identify under which loop)
    %complete_associated_nb = numel(track_stat_ref.ref_tracks_associated) - num_thrown_ref_associated
    
    
    %fprintf('In fraction:  %.3f persistant reference tracks, %.3f timepoint-touching reference tracks,and %.3f complete reference tracks\n', ...
        %ref_unassoc_persistant, num_thrown_ref / numel(track_stat_ref.ref_tracks_unassociated), ...
        %(numel(track_stat_ref.ref_tracks_unassociated) - num_thrown_ref) / numel(track_stat_ref.ref_tracks_unassociated));
    %fprintf('In fraction:  %.3f persistant correlated tracks, %.3f  timepoint-touching correlated tracks,and %.3f complete correlated tracks\n', ...
        %cor_unassoc_persistant, num_thrown_cor / numel(track_stat_cor.cor_tracks_unassociated), ...
        %(numel(track_stat_cor.cor_tracks_unassociated) - num_thrown_cor) / numel(track_stat_cor.cor_tracks_unassociated));
    %fprintf('In fraction:  %.3f  persistant associated reference tracks, %.3f  timepoint-touching associated reference tracks,and %.3f complete associated reference tracks\n', ...
        %ref_assoc_persistant, num_thrown_ref_associated / numel(track_stat_ref.ref_tracks_associated), ...
        %(numel(track_stat_ref.ref_tracks_associated) - num_thrown_ref_associated) / numel(track_stat_ref.ref_tracks_associated));
    %fprintf('In fraction:  %.3f persistant associated correlated tracks, %.3f  timepoint-touching associated correlated tracks,and %.3f complete associated correlated tracks\n', ...
        %cor_assoc_persistant, num_thrown_cor_associated / numel(track_stat_ref.ref_tracks_associated), ...
        %(numel(track_stat_ref.ref_tracks_associated) - num_thrown_cor_associated) / numel(track_stat_ref.ref_tracks_associated));
    
    
    %these criteria:
    % I is > 0 if there was a "low timepoint violator" in the associated reference
    % tracks
    
    % J is > 0 if there was a "high timepoint violator" in the assocaited
    % reference tracks
    
    if or(I>0, J>0)  % now delete appropriate vector entry to keep lists pure
        times_flagged_ref=times_flagged_ref+1;
        
        % diagnositc
%         i
%         I
%         J
%         num_thrown_tsrrta
%         num_thrown_ref
        
        % add a conditional here to prevent a crash when
        % num_thrown_tsrrta>0 (when there are tracks in track stat ref that
        % violoate timepoints but none happen to be associated)
        
%         This finds the first element in the list (except for and deletes it.
            
            track_stat_ref.ref_tracks_associated(i-num_thrown_tsrrta,:)=[];
            num_thrown_tsrrta=num_thrown_tsrrta+1;
        
        %         for j=1:numel(ROI_track_ref)
        %             numkdel=0;
        %             for k=1:size(ROI_track_ref(j).associated_ref_track,2)
        %                 if (ROI_track_ref(j).associated_ref_track(k-numkdel)==aa)
        %                     ROI_track_ref(j).associated_ref_track(k-numkdel)=[];
        %                     numkdel=numkdel+1;
        %                     num_thrown_roiart=num_thrown_roiart+1;
        %                 end
        %             end
        %         end
        %
        %
        %         numdel=0;
        %         for j=1:numel(track_stat_ref.ref_primary_track-num_thrown_tsrrpt)
        %             if (track_stat_ref.ref_primary_track(j-numdel)==aa)
        %                 track_stat_ref.ref_primary_track(j-numdel)= [];
        %                 numdel=numdel+1;
        %                 num_thrown_tsrrpt=num_thrown_tsrrpt+1;
        %             end
        %         end
        
    end
    
    %     if or(K>0, L>0)
    %         times_flagged_cor=times_flagged_cor+1;
    %         for j=1:numel(ROI_track_ref)
    %             for l=1:size(bb,2)
    %                 numkdel=0;
    %                 for k=1:size(ROI_track_ref(j).associated_cor_track,2)
    %                     if (ROI_track_ref(j).associated_cor_track(k-numkdel)==bb(l))
    %                         ROI_track_ref(j).associated_cor_track(k-numkdel)=[];
    %                         numkdel=numkdel+1;
    %                         num_thrown_cor=num_thrown_cor+1;
    %                     end
    %                 end
    %             end
    %         end
    %     end
end
%% save number of persistent tracks
 ref_unassoc_persistant_nb = numel(ref_unassoc_persistant_tracks);
    ref_assoc_persistant_nb = numel(ref_assoc_persistant_tracks);
    cor_unassoc_persistant_nb = numel(cor_unassoc_persistant_tracks);
    cor_assoc_persistant_nb = numel(cor_assoc_persistant_tracks);
    
    ref_unassoc_persistant = numel(ref_unassoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_unassociated);
    ref_assoc_persistant = numel(ref_assoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_associated);
    cor_unassoc_persistant = numel(cor_unassoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_unassociated);
    cor_assoc_persistant = numel(cor_assoc_persistant_tracks) / numel(track_stat_ref.ref_tracks_unassociated);
    
    % change by CK 180802, reference to the wrong variable
    %num_thrown_ref_associated = highindexa + lowindexa - numel(ref_assoc_persistant);
    %num_thrown_cor_associated = highindexb + lowindexb - numel(cor_assoc_persistant);
    
    num_thrown_ref_associated = highindexa + lowindexa - numel(ref_assoc_persistant_tracks);
    num_thrown_cor_associated = highindexb + lowindexb - numel(cor_unassoc_persistant_tracks);
    
    % complete_associated_nb is in a loop and variable is overwritten and thus output wrong in CME analysis (so far cannot identify under which loop)
    complete_associated_nb = numel(track_stat_ref.ref_tracks_associated) - num_thrown_ref_associated




%% save montages of the corresponding tracks that are not associated with ref.
unassociated_tracks = struct('ref_movie', {},'cor_movie', {});

if ~isempty(undeleted_tscctu)

for cor_track = undeleted_tscctu
    
    ROI_x = ROI_track_cor(cor_track,1).ROI_x;
    ROI_y = ROI_track_cor(cor_track,1).ROI_y;
    
    ROI_delx = ROI_x(1,2) - ROI_x(1,1) + 1;
    ROI_dely = ROI_y(1,2) - ROI_y(1,1) + 1;
    
    tp_prev = ROI_track_cor(cor_track, 1).tp_prev;
    tp_post = ROI_track_cor(cor_track, 1).tp_post;
    
    
    ROI_x_range = ROI_x(1):ROI_x(2);
    ROI_y_range = ROI_y(1):ROI_y(2);
    
    unassociated_tracks(cor_track, 1).cor_movie = reshape(cor_stack_3D(ROI_y_range, ROI_x_range, tp_prev:tp_post), [ROI_dely, ROI_delx*(tp_post - tp_prev+ 1), 1]);
    unassociated_tracks(cor_track, 1).ref_movie = reshape(ref_stack_3D(ROI_y_range - shift_y, ROI_x_range - shift_x,...
        tp_prev:tp_post), [ROI_dely, ROI_delx*(tp_post - tp_prev+ 1), 1]);
    
    unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected = NaN*ones(2, (tp_post - tp_prev+ 1));
    unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected = NaN*ones(2, (tp_post - tp_prev+ 1));
    
    
    unassociated_tracks(cor_track, 1).cor_montage = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
    unassociated_tracks(cor_track, 1).ref_montage = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
    unassociated_tracks(cor_track, 1).cor_intensity = NaN*ones(2, (tp_post - tp_prev+ 1));
    unassociated_tracks(cor_track, 1).ref_intensity = NaN*ones(2, (tp_post - tp_prev+ 1));
    unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected = NaN*ones(2, (tp_post - tp_prev+ 1));
    unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected = NaN*ones(2, (tp_post - tp_prev+ 1));
    
    unassociated_tracks(cor_track, 1).cor_montage_bg = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
    unassociated_tracks(cor_track, 1).ref_montage_bg = NaN*ones((2*montage_size+1),(2*montage_size+1)*(tp_post - tp_prev+ 1));
    
    
    
    % Determine the integrated intensity
    if ~any(any(isnan(ROI_track_cor(cor_track).cor_coord_extrapol)))
        for tp = ROI_track_cor(cor_track,1).tp_prev:ROI_track_cor(cor_track,1).tp_post
            rel_tp = (tp - ROI_track_cor(cor_track,1).tp_prev + 1);
            
            % save corresponding montage
            peak_x = round(ROI_track_cor(cor_track).cor_coord_extrapol(rel_tp, 1)/pixel_size);
            peak_y = round(ROI_track_cor(cor_track).cor_coord_extrapol(rel_tp, 2)/pixel_size);
            
            peak_region_image = cor_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp).*peak_region_mask;
            unassociated_tracks(cor_track, 1).cor_intensity(1, rel_tp) = sum(sum(peak_region_image));
            unassociated_tracks(cor_track, 1).cor_intensity(2, rel_tp) = std2(peak_region_image);
            unassociated_tracks(cor_track, 1).cor_montage(:,(rel_tp*7-6):(rel_tp*7)) = peak_region_image;
            
            bg_region_image = cor_stack_3D(peak_y-bg_size:peak_y+bg_size, peak_x-bg_size:peak_x+bg_size, tp).*bg_region_mask;
            bg_value = median(bg_region_image(bg_region_mask));
            peak_region_image_bg_corrected = (cor_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp)...
                - bg_value).*peak_region_mask;
            unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, rel_tp) = sum(sum(peak_region_image_bg_corrected));
            unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(2, rel_tp) = std2(peak_region_image_bg_corrected);
            unassociated_tracks(cor_track, 1).cor_montage_bg(:,(rel_tp*7-6):(rel_tp*7)) = peak_region_image;
            
            
            
            if ~any(isnan(ROI_track_cor(cor_track).ref_coord_extrapol))
                % save reference montage
                peak_x = round(ROI_track_cor(cor_track).ref_coord_extrapol(rel_tp, 1)/pixel_size);
                peak_y = round(ROI_track_cor(cor_track).ref_coord_extrapol(rel_tp, 2)/pixel_size);
                
                peak_region_image = ref_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp).*peak_region_mask;
                
                unassociated_tracks(cor_track, 1).ref_intensity(1, rel_tp) = sum(sum(peak_region_image));
                unassociated_tracks(cor_track, 1).ref_intensity(2, rel_tp) = std2(peak_region_image);
                unassociated_tracks(cor_track, 1).ref_montage(:,(rel_tp*7-6):(rel_tp*7)) = peak_region_image;
                
                bg_region_image = ref_stack_3D(peak_y-bg_size:peak_y+bg_size, peak_x-bg_size:peak_x+bg_size, tp).*bg_region_mask;
                bg_value = median(bg_region_image(bg_region_mask));
                peak_region_image_bg_corrected = (ref_stack_3D(peak_y-montage_size:peak_y+montage_size, peak_x-montage_size:peak_x+montage_size, tp)...
                    - bg_value).*peak_region_mask;
                unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, rel_tp) = sum(sum(peak_region_image_bg_corrected));
                unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(2, rel_tp) = std2(peak_region_image_bg_corrected);
            end
            
        end
    end
end
end
%% CD save old files that won't be modified in clean tracks

for i=1:numel(unassociated_tracks)
    if (isempty(unassociated_tracks(i).ref_intensity)==1)
        track_stat_ref.old_unassoc_ref_intensity(i).a=0;
    else
        track_stat_ref.old_unassoc_ref_intensity(i).a=unassociated_tracks(i).ref_intensity;
    end
    if (isempty(unassociated_tracks(i).cor_intensity)==1)
        track_stat_ref.old_unassoc_cor_intensity(i).a=0;
    else
        track_stat_ref.old_unassoc_cor_intensity(i).a=unassociated_tracks(i).cor_intensity;
    end
end

for i=1:numel(associated_tracks)
    if (isempty(associated_tracks(i).ref_intensity)==1)
        track_stat_ref.old_assoc_ref_intensity(i).a=0;
    else
        track_stat_ref.old_assoc_ref_intensity(i).a=associated_tracks(i).ref_intensity;
    end
    
    if (isempty(associated_tracks(i).cor_intensity)==1)
        track_stat_ref.old_assoc_cor_intensity(i).a=0;
    else
        track_stat_ref.old_assoc_cor_intensity(i).a=associated_tracks(i).cor_intensity;
    end
end



%% save file names
analysis_info.ref_csv_name = fullfile(ref_folder_name, ref_file);
analysis_info.cor_csv_name = fullfile(cor_folder_name, cor_file);
analysis_info.ref_image_name = fullfile(ref_img_folder_name, ref_img_file);
analysis_info.cor_image_name = fullfile(cor_img_folder_name, cor_img_file);
analysis_info.tp_max = tp_max;
analysis_info.pixel_size = pixel_size;
analysis_info.shift_y = shift_y;
analysis_info.shift_x = shift_x;
analysis_info.shift_y_nm = shift_y_nm;
analysis_info.shift_x_nm = shift_x_nm;
analysis_info.extra_tp = extra_tp;
%%
parent_folder = cor_img_folder_name(1:end - 4);
save(fullfile(parent_folder, 'associated_tracks_ref_cor.mat'),...
    'position_data_ref_original', 'position_data_cor_original', ...
    'associated_tracks', 'unassociated_tracks', 'track_stat_ref', 'track_stat_cor', 'ROI_track_ref', 'ROI_track_cor', 'analysis_info');

display(['The association results are saved as -', fullfile(cor_folder_name, 'associated_tracks_ref_cor.mat')]);
%%
% [y, Fs] = audioread('notify.wav');
% sound(y, Fs);
% 
% 
% 
display('debug nb of complete_associated_nb ');
complete_associated_nb
% display('The number of unassociated reference tracks within timepoint boundaries is');
% numel(track_stat_ref.ref_tracks_unassociated)
% display('The number of associated corresponding tracks within timepoint boundaries is');
% numel(unique([ROI_track_ref(track_stat_ref.ref_tracks_associated).associated_cor_track]))
% display('The number of unassociated corresponding tracks within timepoint boundaries is');
% numel(track_stat_cor.cor_tracks_unassociated)
% 
% display('The number of associated reference tracks above timepoint boundaries is');
% numel(track_stat_ref.ref_tracks_associated_high_breach)
% display('The number of unassociated reference tracks above timepoint boundaries is');
% numel(track_stat_ref.ref_tracks_unassociated_high_breach)
% display('The number of associated corresponding tracks above timepoint boundaries is');
% %numel(track_stat_cor.cor_tracks_associated_high_breach)
% display('The number of unassociated corresponding tracks above timepoint boundaries is');
% numel(track_stat_cor.cor_tracks_unassociated_high_breach)
% 
% display('The number of associated reference tracks below timepoint boundaries is');
% numel(track_stat_ref.ref_tracks_associated_low_breach)
% display('The number of unassociated reference tracks below timepoint boundaries is');
% numel(track_stat_ref.ref_tracks_unassociated_low_breach)
% display('The number of associated corresponding tracks below timepoint boundaries is');
% numel(track_stat_cor.cor_tracks_associated_low_breach)
% display('The number of unassociated corresponding tracks below timepoint boundaries is');
% numel(track_stat_cor.cor_tracks_unassociated_low_breach)