%% This program reads in a mat file of associated tracks determined with Associate_tracks.m
%% and manual_pick file from Manually_pick_associtated_tracks.m
%
%%

% Log
% 3/10/3015
% Asks user which of the following she wants to use for analysis:
% manually_cleaned_track_list.mat or cleaned_track_list.mat

% 4/30/2015
% 1) Reduced the number of plots.
% 2) Changed the size and the positioning of the plots.
% 3) Defined drawAveragePlots to simplify code.
% 4) Fixed the error in 'lifetime_max_intensity_assoc_cor
%         from
%         lifetime_max_intensity_assoc_cor(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).tp_last - ROI_track_ref(ref_track, 1).tp_first + 1)*interval, ...
%                         max(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :))];
%         to
%         lifetime_max_intensity_assoc_cor(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval, ...
%                         max(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :))];
% 5) Added in 'track_stat_asso,' a matrix for statistical information of
% the associated tracks. 'track_stat_asso' is saved in
% 'track_statistics.txt' file. along with 'track_stat_ref_only' and
% 'track_stat_cor_only.'
% 6) Output file name 'track_stat_associated_only.txt' was changed to 'track_statistics.txt'

% 10/2/2015
% MA
% Changed 'suplabel' on graphs so that title appears (added 't')

% 1/12/16
% MA
% Accommodates 0,1,2,3 option for track label, from "manually pick tracks"
% program and plots as individual subplots

% 5/31/16
% CD
% Outputs the lifetime_associated data and the field_positions to identify
% fields within the lifetime_associated data. This information can then be
% used to compare stats with an Anderson Darling test using the
% "AnDarksamtest.m" function


function [field_positions, lifetime_associated]= Plot_stats20151006()
%% Mark 'true' if you want to save all the variables

saveVariables = true;

%% Mark 'true' if you want to plot the "unusual" and "special case" tracks (#2 and 3)

plot_2and3_tracks = false;
plot_4hilo = false;         % MA this seems broken.

analyzeUnassociated = 0; % decide whether to plot and analyze unassociated tracks

%% Ask users for interval
interval = input('What was the interval between time points?')

%% Initialize variables for all files
% filenamestp_
foldername_all = struct;


% all associated, unassociated tracks

jj = 1;
kk = 1;
ll = 1;
jjj = 1;

allAssociatedTracks      = [];
allUnassociatedTracksRef = [];
allUnassociatedTracksCor = [];

% lifetime matrices
lifetime_ref_only  = [];
lifetime_cor_only  = [];
field_positions = [];
lifetime_associated= [];

lifetime_ref_only_2  = [];
lifetime_cor_only_2  = [];
lifetime_associated_2= [];

lifetime_ref_only_3   = [];
lifetime_cor_only_3   = [];
lifetime_associated_3 = [];

lifetime_ref_only_4_high   = [];
lifetime_cor_only_4_high   = [];
lifetime_associated_4_high = [];

lifetime_ref_only_4_low   = [];
lifetime_cor_only_4_low   = [];
lifetime_associated_4_low = [];

% profile matrices
%1
cor_aligned2cormax = [];
ref_aligned2cormax = [];

ref_only_ref_aligned2cormax = [];
ref_only_cor_aligned2cormax = [];

cor_only_ref_aligned2cormax = [];
cor_only_cor_aligned2cormax = [];
%2
cor_aligned2refmax = [];
ref_aligned2refmax = [];

ref_only_ref_aligned2refmax = [];
ref_only_cor_aligned2refmax = [];

cor_only_ref_aligned2refmax = [];
cor_only_cor_aligned2refmax = [];

%3
cor_aligned2ref_begin = [];
ref_aligned2ref_begin = [];

ref_only_ref_aligned2ref_begin = [];
ref_only_cor_aligned2ref_begin = [];

ref_x1_aligned2ref_begin = [];
ref_x2_aligned2ref_begin = [];

%4
cor_aligned2ref_end = [];
ref_aligned2ref_end = [];

ref_only_ref_aligned2ref_end = [];
ref_only_cor_aligned2ref_end = [];

%5
cor_aligned2cor_begin = [];
ref_aligned2cor_begin = [];

cor_only_ref_aligned2cor_begin = [];
cor_only_cor_aligned2cor_begin = [];

cor_x1_aligned2cor_begin = [];
cor_x2_aligned2cor_begin = [];


%6
cor_aligned2cor_end = [];
ref_aligned2cor_end = [];

cor_only_ref_aligned2cor_end = [];
cor_only_cor_aligned2cor_end = [];

%7
ref_only_x1_aligned2begin = [];
ref_only_x2_aligned2begin = [];

cor_only_x1_aligned2begin = [];
cor_only_x2_aligned2begin = [];


% Correlation of lifetime and maximum intensity.
lifetime_max_intensity_assoc_ref= [];
lifetime_max_intensity_ref_only = [];
lifetime_mean_intensity_assoc_ref= [];
lifetime_mean_intensity_ref_only = [];

lifetime_max_intensity_assoc_cor= [];
lifetime_max_intensity_cor_only = [];
lifetime_mean_intensity_assoc_cor= [];
lifetime_mean_intensity_cor_only = [];

% Stats for the ref-only and cor-only tracks
track_stat_ref_only = [];
track_stat_cor_only = [];
track_stat_asso = [];


%% Input folder for all data.
load('previous_folder.mat')
data_folder = uigetdir(previous_folder, 'Pick folder of all analysis.');
if isdir(data_folder)
    previous_folder = fullfile(data_folder, '');
    save('previous_folder', 'previous_folder')
else
    display('You did not pick a folder.')
    return
end


%% Asks the user which of the following she wants to use for analysis:
% manually_cleaned_track_list.mat or cleaned_track_list.mat

clean_list_choice = input(['Which of the following two do you want to use?.\n',...
    '1) manually_cleaned_track_list.mat \n2) cleaned_track_list.mat (automatic cleanup)\n',...
    'Type 1 or 2.\n']);
if clean_list_choice == 1
    clean_list_name = 'manually_cleaned_track_list.txt';
elseif clean_list_choice == 2
    clean_list_name = 'cleaned_track_list.mat';%'cleaned_track_list.mat';
    
else
    error('The answer needs to be either 1 or 2.')
    return
end

%% Name of the results files
lifetime_intensity_file = fullfile(data_folder, 'lifetime_intensity.txt');
track_stat_file_name = fullfile(data_folder, 'track_statistics.txt');
histogram_file_name = fullfile(data_folder, 'histogram_pool_data.txt');

%% Go through folders to pool data.
folder_ID1 = 0;
list1 = dir(data_folder);
for foldernum1 = 3:numel(list1)
    if list1(foldernum1).isdir
        list2 = dir(fullfile(data_folder, list1(foldernum1).name));
        folder_ID1 = folder_ID1 + 1;
        folder_ID2 = 0;
        for foldernum2 = 3:numel(list2)
            if list2(foldernum2).isdir
                
                % load tracking data and manual tag data.
                track_folder_name = fullfile(data_folder, list1(foldernum1).name, list2(foldernum2).name);
                if exist(fullfile(track_folder_name, 'associated_tracks_ref_cor.mat'), 'file');
                    load(fullfile(track_folder_name, 'associated_tracks_ref_cor.mat'));
                    pixel_size = analysis_info.pixel_size; % in nm
                    if exist(fullfile(track_folder_name, clean_list_name), 'file');
                        machine_tag_folder = track_folder_name;
                        if strcmp(clean_list_name((end-2):end), 'mat')
                            load(fullfile(track_folder_name, clean_list_name));
                        elseif strcmp(clean_list_name((end-2):end), 'txt')
                            track_flag_list = dlmread(fullfile(track_folder_name, clean_list_name), '\t', 10, 0);
                        end
                    else
                        display([fullfile(track_folder_name, clean_list_name), ' is not found.']);
                    end
                else
                    display([fullfile(track_folder_name, 'associated_tracks_ref_cor.mat'), ' is not found.']);
                    continue
                end
                folder_ID2 = folder_ID2 + 1;
                
                %% Record file names
                foldername_all(folder_ID1, folder_ID2).associated_track_folder = track_folder_name;
                foldername_all(folder_ID1, folder_ID2).associated_track_file = 'associated_tracks_ref_cor.mat';
                foldername_all(folder_ID1, folder_ID2).machine_tag_folder = machine_tag_folder;
                foldername_all(folder_ID1, folder_ID2).machine_tag_file = clean_list_name;
                
                %% Retrieve good tracks.
                ref_cor_tracks  = track_flag_list(and(track_flag_list(:,1) == 1, track_flag_list(:,3) == 1), 2);
                ref_only_tracks = track_flag_list(and(track_flag_list(:,1) == 2, track_flag_list(:,3) == 1),2);
                cor_only_tracks = track_flag_list(and(track_flag_list(:,1) == 3, track_flag_list(:,3) == 1),2);
                
                cor_only_num = numel(cor_only_tracks);
                ref_only_num = numel(ref_only_tracks);
                ref_cor_num = numel(ref_cor_tracks);
                
                %% Retrieve "unusual" tracks (track ID = 2)
                %% casey copy here for _4
                ref_cor_tracks_2  = track_flag_list(and(track_flag_list(:,1) == 1, track_flag_list(:,3) == 2), 2);
                ref_only_tracks_2  = track_flag_list(and(track_flag_list(:,1) == 2, track_flag_list(:,3) == 2), 2);
                cor_only_tracks_2  = track_flag_list(and(track_flag_list(:,1) == 3, track_flag_list(:,3) == 2), 2);
                
                ref_cor_2_num     = numel(ref_cor_tracks_2);
                ref_only_2_num    = numel(ref_only_tracks_2);
                cor_only_2_num    = numel(cor_only_tracks_2);
                
                %% Retrieve "special case" tracks (track ID = 3)
                
                ref_cor_tracks_3  = track_flag_list(and(track_flag_list(:,1) == 1, track_flag_list(:,3) == 3), 2);
                ref_only_tracks_3 = track_flag_list(and(track_flag_list(:,1) == 2, track_flag_list(:,3) == 3),2);
                cor_only_tracks_3 = track_flag_list(and(track_flag_list(:,1) == 3, track_flag_list(:,3) == 3),2);
                
                ref_cor_3_num     = numel(ref_cor_tracks_3);
                ref_only_3_num    = numel(ref_only_tracks_3);
                cor_only_3_num    = numel(cor_only_tracks_3);
                
                 %% Retrieve "timepoint violators" tracks
                % There are 2 subcategories: early (low) and late (high) tracks
                
                %MA: if the field doesn't exist (i.e. no tracks came from
                %this categroy) then make an empty vector.
                
                if isfield(track_stat_ref, 'ref_tracks_associated_low_breach')
                    ref_cor_tracks_4_low = [track_stat_ref.ref_tracks_associated_low_breach]';
                else
                    ref_cor_tracks_4_low = [];
                end
                
                if isfield(track_stat_ref, 'ref_tracks_unassociated_low_breach')
                    ref_only_tracks_4_low = [track_stat_ref.ref_tracks_unassociated_low_breach]';
                else
                    ref_only_tracks_4_low = [];
                end
                
                if isfield(track_stat_cor, 'cor_tracks_unassociated_low_breach')
                    cor_only_tracks_4_low = [track_stat_cor.cor_tracks_unassociated_low_breach]';
                else
                    cor_only_tracks_4_low = [];
                end
                
                ref_cor_tracks_4_low_num = numel(ref_cor_tracks_4_low);
                ref_only_tracks_4_low_num = numel(ref_only_tracks_4_low);
                cor_only_tracks_4_low_num = numel(cor_only_tracks_4_low);
                
                
                if isfield(track_stat_ref, 'ref_tracks_associated_high_breach')
                    ref_cor_tracks_4_high = [track_stat_ref.ref_tracks_associated_high_breach]';
                else
                    ref_cor_tracks_4_high = [];
                end
                
                if isfield(track_stat_ref, 'ref_tracks_unassociated_high_breach')
                    ref_only_tracks_4_high = [track_stat_ref.ref_tracks_unassociated_high_breach]';
                else
                    ref_only_tracks_4_high = [];
                end
                
                if isfield(track_stat_cor, 'cor_tracks_unassociated_high_breach')
                    cor_only_tracks_4_high = [track_stat_cor.cor_tracks_unassociated_high_breach]';
                else
                    cor_only_tracks_4_high = [];
                end
                
                ref_cor_tracks_4_high_num = numel(ref_cor_tracks_4_high);
                ref_only_tracks_4_high_num = numel(ref_only_tracks_4_high);
                cor_only_tracks_4_high_num = numel(cor_only_tracks_4_high);
                
                
                %% Retrive variables.
                tp_min = 1;
                tp_max = analysis_info.tp_max;
                
                
                %% Compile ref_cor and ref_ indecies, intensities
                for i = 1:ref_cor_num
                    cur_track = ref_cor_tracks(i);
                    
                    curCorTrack = (associated_tracks(cur_track, 1).cor_intensity_bg_corrected(1, :));
                    curRefTrack = (associated_tracks(cur_track, 1).ref_intensity_bg_corrected(1, :));
                    
                    allAssociatedTracks(jj).Cor = curCorTrack;
                    allAssociatedTracks(jj).Ref = curRefTrack;
                    jj = jj+1;
                end
                                
                
                curTrack = [];
                curCorTrack = [];
                curRefTrack = [];
                
                if analyzeUnassociated
                
                if ~isempty(unassociated_tracks)&&~isempty(ref_only_tracks)
                    
                    for i = 1:ref_only_num
                        cur_track = ref_only_tracks(i);
                        
                        curCorTrack = unassociated_tracks(cur_track,1).cor_intensity_bg_corrected(1,:);
                        curRefTrack = unassociated_tracks(cur_track,1).ref_intensity_bg_corrected(1,:);
                        
                        allUnassociatedTracksRef(kk).Cor = curCorTrack;
                        allUnassociatedTracksRef(kk).Ref = curRefTrack;
                        kk = kk+1;
                    end
                    
                else 
                    disp('Did not find any unassociaed tracks.')
                end
                end
                curTrack = [];
                curCorTrack = [];
                curRefTrack = [];
                
                if ~isempty(unassociated_tracks)&&~isempty(cor_only_tracks)
                    
                    for i = 1:cor_only_num
                        cur_track = cor_only_tracks(i);
                        
                        curCorTrack = unassociated_tracks(cur_track,1).cor_intensity_bg_corrected(1,:);
                        curRefTrack = unassociated_tracks(cur_track,1).ref_intensity_bg_corrected(1,:);
                        
                        allUnassociatedTracksCor(ll).Cor = curCorTrack;
                        allUnassociatedTracksCor(ll).Ref = curRefTrack;
                        ll = ll+1;
                    end
                end
                %% Determine lifetimes add it to the matrices
                previous_index = size(lifetime_associated, 1);
                lifetime_associated= [lifetime_associated; NaN*ones(ref_cor_num, 2)];
                field_positions=[field_positions; ref_cor_num];
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    cor_track = ROI_track_ref(ref_track, 1).associated_cor_track;
                    %%
                    lifetime_associated(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).ref_tp_last - ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval, ...
                        %%
                        (ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval];
                end
                
                previous_index = size(lifetime_ref_only, 1);
                lifetime_ref_only = [lifetime_ref_only; NaN*ones(ref_only_num, 1)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    lifetime_ref_only(previous_index + i,1) = (ROI_track_ref(ref_track, 1).tp_last - ROI_track_ref(ref_track, 1).tp_first + 1)*interval;
                end
                
                previous_index = size(lifetime_cor_only, 1);
                lifetime_cor_only = [lifetime_cor_only; NaN*ones(cor_only_num, 1)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    lifetime_cor_only(previous_index + i,1) =(ROI_track_cor(cor_track, 1).tp_last - ROI_track_cor(cor_track, 1).tp_first + 1)*interval;
                end
                
                %% Determine lifetimes of unusual tracks
                %% casey copy here for _4
                
                previous_index = size(lifetime_associated_2, 1);
                lifetime_associated_2= [lifetime_associated_2; NaN*ones(ref_cor_2_num, 2)];
                for i = 1:ref_cor_2_num
                    ref_track_2 = ref_cor_tracks_2(i);
                    cor_track_2 = ROI_track_ref(ref_track_2, 1).associated_cor_track;
                    lifetime_associated_2(previous_index + i,:) = [(ROI_track_ref(ref_track_2, 1).ref_tp_last - ROI_track_ref(ref_track_2, 1).ref_tp_first + 1)*interval, ...
                        (ROI_track_ref(ref_track_2, 1).cor_tp_last - ROI_track_ref(ref_track_2, 1).cor_tp_first + 1)*interval];
                end
                
                previous_index = size(lifetime_ref_only_2, 1);
                lifetime_ref_only_2 = [lifetime_ref_only_2; NaN*ones(ref_only_2_num, 1)];
                for i = 1:ref_only_2_num
                    ref_track_2 = ref_only_tracks_2(i);
                    lifetime_ref_only_2(previous_index + i,1) = (ROI_track_ref(ref_track_2, 1).tp_last - ROI_track_ref(ref_track_2, 1).tp_first + 1)*interval;
                end
                
                previous_index = size(lifetime_cor_only_2, 1);
                lifetime_cor_only_2 = [lifetime_cor_only_2; NaN*ones(cor_only_2_num, 1)];
                for i = 1:cor_only_2_num
                    cor_track_2 = cor_only_tracks_2(i);
                    lifetime_cor_only_2(previous_index + i,1) =(ROI_track_cor(cor_track_2, 1).tp_last - ROI_track_cor(cor_track_2, 1).tp_first + 1)*interval;
                end
                
                
                %% Determine lifetimes of special case tracks
                previous_index = size(lifetime_associated_3, 1);
                lifetime_associated_3= [lifetime_associated_3; NaN*ones(ref_cor_3_num, 2)];
                for i = 1:ref_cor_3_num
                    ref_track_3 = ref_cor_tracks_3(i);
                    cor_track_3 = ROI_track_ref(ref_track_3, 1).associated_cor_track;
                    lifetime_associated_3(previous_index + i,:) = [(ROI_track_ref(ref_track_3, 1).ref_tp_last - ROI_track_ref(ref_track_3, 1).ref_tp_first + 1)*interval, ...
                        (ROI_track_ref(ref_track_3, 1).cor_tp_last - ROI_track_ref(ref_track_3, 1).cor_tp_first + 1)*interval];
                end
                
                previous_index = size(lifetime_ref_only_3, 1);
                lifetime_ref_only_3 = [lifetime_ref_only_3; NaN*ones(ref_only_3_num, 1)];
                for i = 1:ref_only_3_num
                    ref_track_3 = ref_only_tracks_3(i);
                    lifetime_ref_only_3(previous_index + i,1) = (ROI_track_ref(ref_track_3, 1).tp_last - ROI_track_ref(ref_track_3, 1).tp_first + 1)*interval;
                end
                
                previous_index = size(lifetime_cor_only_3, 1);
                lifetime_cor_only_3 = [lifetime_cor_only_3; NaN*ones(cor_only_3_num, 1)];
                for i = 1:cor_only_3_num
                    cor_track_3 = cor_only_tracks_3(i);
                    lifetime_cor_only_3(previous_index + i,1) =(ROI_track_cor(cor_track_3, 1).tp_last - ROI_track_cor(cor_track_3, 1).tp_first + 1)*interval;
                end
                
                %% Determine lifetimes of timepoint violators for high timepoints
                
                previous_index = size(lifetime_associated_4_high, 1);
                lifetime_associated_4_high= [lifetime_associated_4_high; NaN*ones(ref_cor_tracks_4_high_num, 2)];
                for i = 1:ref_cor_tracks_4_high_num
                    ref_track_4_high = ref_cor_tracks_4_high(i);
                    cor_track_4_high = ROI_track_ref(ref_track_4_high, 1).associated_cor_track;
                    lifetime_associated_4_high(previous_index + i,:) = [(ROI_track_ref(ref_track_4_high, 1).ref_tp_last - ROI_track_ref(ref_track_4_high, 1).ref_tp_first + 1)*interval, ...
                        (ROI_track_ref(ref_track_4_high, 1).cor_tp_last - ROI_track_ref(ref_track_4_high, 1).cor_tp_first + 1)*interval];
                end
                
                previous_index = size(lifetime_ref_only_4_high, 1);
                lifetime_ref_only_4_high = [lifetime_ref_only_4_high; NaN*ones(ref_only_tracks_4_high_num, 1)];
                for i = 1:ref_only_tracks_4_high_num
                    ref_track_4_high = ref_only_tracks_4_high(i);
                    lifetime_ref_only_4_high(previous_index + i,1) = (ROI_track_ref(ref_track_4_high, 1).tp_last - ROI_track_ref(ref_track_4_high, 1).tp_first + 1)*interval;
                end
                
                previous_index = size(lifetime_cor_only_4_high, 1);
                lifetime_cor_only_4_high = [lifetime_cor_only_4_high; NaN*ones(cor_only_tracks_4_high_num, 1)];
                for i = 1:cor_only_tracks_4_high_num
                    cor_track_4_high = cor_only_tracks_4_high(i);
                    lifetime_cor_only_4_high(previous_index + i,1) =(ROI_track_cor(cor_track_4_high, 1).tp_last - ROI_track_cor(cor_track_4_high, 1).tp_first + 1)*interval;
                end
                
                %% Determine lifetimes of timepoint violators for low timepoints
                
                previous_index = size(lifetime_associated_4_low, 1);
                lifetime_associated_4_low= [lifetime_associated_4_low; NaN*ones(ref_cor_tracks_4_low_num, 2)];
                for i = 1:ref_cor_tracks_4_low_num
                    ref_track_4_low = ref_cor_tracks_4_low(i);
                    cor_track_4_low = ROI_track_ref(ref_track_4_low, 1).associated_cor_track;
                    lifetime_associated_4_low(previous_index + i,:) = [(ROI_track_ref(ref_track_4_low, 1).ref_tp_last - ROI_track_ref(ref_track_4_low, 1).ref_tp_first + 1)*interval, ...
                        (ROI_track_ref(ref_track_4_low, 1).cor_tp_last - ROI_track_ref(ref_track_4_low, 1).cor_tp_first + 1)*interval];
                end
                
                previous_index = size(lifetime_ref_only_4_low, 1);
                lifetime_ref_only_4_low = [lifetime_ref_only_4_low; NaN*ones(ref_only_tracks_4_low_num, 1)];
                for i = 1:ref_only_tracks_4_low_num
                    ref_track_4_low = ref_only_tracks_4_low(i);
                    lifetime_ref_only_4_low(previous_index + i,1) = (ROI_track_ref(ref_track_4_low, 1).tp_last - ROI_track_ref(ref_track_4_low, 1).tp_first + 1)*interval;
                end
                
                previous_index = size(lifetime_cor_only_4_low, 1);
                lifetime_cor_only_4_low = [lifetime_cor_only_4_low; NaN*ones(cor_only_tracks_4_low_num, 1)];
                for i = 1:cor_only_tracks_4_low_num
                    cor_track_4_low = cor_only_tracks_4_low(i);
                    lifetime_cor_only_4_low(previous_index + i,1) =(ROI_track_cor(cor_track_4_low, 1).tp_last - ROI_track_cor(cor_track_4_low, 1).tp_first + 1)*interval;
                end
                
                
                %% Save montage as tif file.
                montage_folder_name = fullfile(track_folder_name, 'montage');
                if ~exist(montage_folder_name, 'dir')
                    mkdir(montage_folder_name);
                else
                    rmdir(montage_folder_name,'s');
                    mkdir(montage_folder_name);
                end
                
                movie_folder_name = fullfile(track_folder_name, 'movie');
                if ~exist(movie_folder_name, 'dir')
                    mkdir(movie_folder_name);
                else
                    rmdir(movie_folder_name,'s');
                    mkdir(movie_folder_name);
                end
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated', num2str(ref_track, '%04i'),'.tif']),'tif');
                end
                
                % as PNG too
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_', num2str(ref_track, '%04i'),'.png']),'png');
                end
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_ref_only', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_ref_only', num2str(ref_track, '%04i'),'.tif']),'tif');
                end
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    imwrite([fit_contrast_ignore_zero(unassociated_tracks(cor_track, 1).ref_montage); fit_contrast_ignore_zero(unassociated_tracks(cor_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_cor_only', num2str(cor_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(unassociated_tracks(cor_track, 1).ref_movie); fit_contrast_ignore_zero(unassociated_tracks(cor_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_cor_only', num2str(cor_track, '%04i'),'.tif']),'tif');
                end
                
                % save montages of "unusual" (#2) tracks %MA
                
                for i = 1:ref_cor_2_num
                    
                    ref_track = ref_cor_tracks_2(i);
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_unusual', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_unusual', num2str(ref_track, '%04i'),'.tif']),'tif');
                    
                    % as PNG too
                    
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_unusual_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_unusual_', num2str(ref_track, '%04i'),'.png']),'png');
                end
                
                % save montages of "special case" (#3) tracks %MA
                
                for i = 1:ref_cor_3_num
                    
                    ref_track = ref_cor_tracks_3(i);
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_specialcase', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_specialcase', num2str(ref_track, '%04i'),'.tif']),'tif');
                    
                    % as PNG too
                    
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_specialcase_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_specialcase_', num2str(ref_track, '%04i'),'.png']),'png');
                end
                
                % save montages of high tracks %CD
                
                for i = 1:ref_cor_tracks_4_high_num
                    
                    ref_track = ref_cor_tracks_4_high(i);
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_high', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_high', num2str(ref_track, '%04i'),'.tif']),'tif');
                    
                    % as PNG too
                    
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_high_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_high_', num2str(ref_track, '%04i'),'.png']),'png');
                end
                
                % save montages of low tracks %CD
                
                for i = 1:ref_cor_tracks_4_low_num
                    
                    ref_track = ref_cor_tracks_4_low(i);
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_low', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_low', num2str(ref_track, '%04i'),'.tif']),'tif');
                    
                    % as PNG too
                    
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_low_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_low_', num2str(ref_track, '%04i'),'.png']),'png');
                end
                
                
                %1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to maximum of cor
                previous_index = size(cor_aligned2cormax, 1);
                cor_aligned2cormax = [cor_aligned2cormax; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2cormax = [ref_aligned2cormax; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    %smoothened_curve = medfilt1(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :), 5);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    %         align the curve so that maximum of the curve is coincident
                    %         [max_value, max_index] = max(smoothened_curve);
                    [max_value, max_index] = max(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_index = median(max_index);
                    cor_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2cormax(previous_index +i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                %% Align reference only tracks to maximum of cor
                previous_index = size(ref_only_ref_aligned2cormax, 1);
                ref_only_ref_aligned2cormax = [ref_only_ref_aligned2cormax; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2cormax = [ref_only_cor_aligned2cormax; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    smoothened_ref= medfilt1(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :), 5);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_ref);
                    ref_index = median(max_index);
                    ref_only_cor_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %% Align corresponding only tracks to maximum of cor
                previous_index = size(cor_only_ref_aligned2cormax, 1);
                cor_only_ref_aligned2cormax = [cor_only_ref_aligned2cormax; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2cormax = [cor_only_cor_aligned2cormax; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    smoothened_curve = medfilt1(unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :), 5);
                    profile_length = ROI_track_cor(cor_track, 1).tp_post - ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_curve);
                    ref_index = median(max_index);
                    cor_only_cor_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to maximum of ref
                previous_index = size(cor_aligned2refmax, 1);
                cor_aligned2refmax = [cor_aligned2refmax; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2refmax = [ref_aligned2refmax; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    %         smoothened_curve = medfilt1(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    %         align the curve so that maximum of the curve is coincident
                    %         [max_value, max_index] = max(smoothened_curve);
                    [max_value, max_index] = max(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                    ref_index = median(max_index);
                    cor_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2refmax(previous_index +i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                
                %% Align reference only tracks to maximum of ref
                previous_index = size(ref_only_ref_aligned2refmax, 1);
                ref_only_ref_aligned2refmax = [ref_only_ref_aligned2refmax; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2refmax = [ref_only_cor_aligned2refmax; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    smoothened_ref= medfilt1(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_ref);
                    ref_index = median(max_index);
                    ref_only_cor_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %% Align corresponding only tracks to maximum of ref
                previous_index = size(cor_only_ref_aligned2refmax, 1);
                cor_only_ref_aligned2refmax = [cor_only_ref_aligned2refmax; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2refmax = [cor_only_cor_aligned2refmax; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    smoothened_curve = medfilt1(unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = ROI_track_cor(cor_track, 1).tp_post - ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_curve);
                    ref_index = median(max_index);
                    cor_only_cor_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the beginning of ref
                previous_index = size(cor_aligned2ref_begin, 1);
                cor_aligned2ref_begin = [cor_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2ref_begin = [ref_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align curve to the beginning of the reference track.
                    ref_index = ROI_track_ref(ref_track, 1).ref_tp_first;
                    first_index = ROI_track_ref(ref_track, 1).ref_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    cor_aligned2ref_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2ref_begin(previous_index +i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Coordinates: Align reference tracks to the beginning of reference track
                previous_index = size(ref_x1_aligned2ref_begin, 1);
                ref_x1_aligned2ref_begin = [ref_x1_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max)];
                ref_x2_aligned2ref_begin = [ref_x2_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).ref_tp_last - ROI_track_ref(ref_track, 1).ref_tp_first + 1;
                    first_index = ROI_track_ref(ref_track, 1).ref_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    last_index = ROI_track_ref(ref_track, 1).ref_tp_last - ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    ref_x1_aligned2ref_begin(previous_index + i,1:profile_length) = ROI_track_ref(ref_track,1).ref_coord_extrapol(first_index:last_index,1)';
                    ref_x2_aligned2ref_begin(previous_index + i,1:profile_length) = ROI_track_ref(ref_track,1).ref_coord_extrapol(first_index:last_index,2)';
                end
                
                
                %% Align reference only tracks to beginning of ref
                previous_index = size(ref_only_ref_aligned2ref_begin, 1);
                ref_only_ref_aligned2ref_begin = [ref_only_ref_aligned2ref_begin; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2ref_begin = [ref_only_cor_aligned2ref_begin; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    smoothened_ref= medfilt1(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align curve to the beginning of the reference track.
                    ref_index = ROI_track_ref(ref_track, 1).ref_tp_first;
                    first_index = ROI_track_ref(ref_track, 1).ref_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    ref_only_cor_aligned2ref_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2ref_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the end of reference track
                previous_index = size(cor_aligned2ref_end, 1);
                cor_aligned2ref_end = [cor_aligned2ref_end; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2ref_end = [ref_aligned2ref_end; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index = ROI_track_ref(ref_track, 1).ref_tp_last;
                    last_index_from_end = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).ref_tp_last;
                    cor_aligned2ref_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2ref_end(previous_index +i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                
                %% Align reference only tracks to the end of reference track
                previous_index = size(ref_only_ref_aligned2ref_end, 1);
                ref_only_ref_aligned2ref_end = [ref_only_ref_aligned2ref_end; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2ref_end = [ref_only_cor_aligned2ref_end; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index = ROI_track_ref(ref_track, 1).ref_tp_last;
                    last_index_from_end = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).ref_tp_last;
                    ref_only_cor_aligned2ref_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2ref_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the beginning of corresponding track
                previous_index = size(cor_aligned2cor_begin, 1);
                cor_aligned2cor_begin = [cor_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2cor_begin = [ref_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the beginning of the curve is coincident
                    ref_index = ROI_track_ref(ref_track, 1).cor_tp_first;
                    first_index = ROI_track_ref(ref_track, 1).cor_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    cor_aligned2cor_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2cor_begin(previous_index +i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Coordinates: Align corresponding tracks to the beginning of corresponding track
                previous_index = size(cor_x1_aligned2cor_begin, 1);
                cor_x1_aligned2cor_begin = [cor_x1_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max)];
                cor_x2_aligned2cor_begin = [cor_x2_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).cor_tp_first + 1;
                    first_index = ROI_track_ref(ref_track, 1).cor_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    last_index = ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    cor_x1_aligned2cor_begin(previous_index + i,1:profile_length) = ROI_track_ref(ref_track,1).cor_coord_extrapol(first_index:last_index,1)';
                    cor_x2_aligned2cor_begin(previous_index + i,1:profile_length) = ROI_track_ref(ref_track,1).cor_coord_extrapol(first_index:last_index,2)';
                end
                
                
                %% Align corresponding only tracks to the beginning of corresponding track
                previous_index = size(cor_only_ref_aligned2cor_begin, 1);
                cor_only_ref_aligned2cor_begin = [cor_only_ref_aligned2cor_begin; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2cor_begin = [cor_only_cor_aligned2cor_begin; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    profile_length = ROI_track_cor(cor_track, 1).tp_post - ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that the beginning of the curve is coincident
                    ref_index = ROI_track_cor(cor_track, 1).tp_first;
                    first_index = ROI_track_cor(cor_track, 1).tp_first - ROI_track_cor(cor_track, 1).tp_prev + 1;
                    cor_only_cor_aligned2cor_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2cor_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                    %                     ref_index = ROI_track_ref(ref_track, 1).cor_tp_first;
                    %                     first_index = ROI_track_ref(ref_track, 1).cor_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                end
                
                
                %6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the disappearance of corresponding track
                previous_index = size(cor_aligned2cor_end, 1);
                cor_aligned2cor_end = [cor_aligned2cor_end; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2cor_end = [ref_aligned2cor_end; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index = ROI_track_ref(ref_track, 1).cor_tp_last;
                    last_index_from_end = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).cor_tp_last;
                    cor_aligned2cor_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2cor_end(previous_index +i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                
                %% Align corresponding only tracks to the disappearance of corresponding track
                previous_index = size(cor_only_ref_aligned2cor_end, 1);
                cor_only_ref_aligned2cor_end = [cor_only_ref_aligned2cor_end; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2cor_end = [cor_only_cor_aligned2cor_end; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    profile_length = ROI_track_cor(cor_track, 1).tp_post - ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index =  ROI_track_cor(cor_track, 1).tp_last;
                    last_index_from_end =  ROI_track_cor(cor_track, 1).tp_post -  ROI_track_cor(cor_track, 1).tp_last;
                    cor_only_cor_aligned2cor_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2cor_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Coordinates: Align reference only tracks to the beginning of reference track
                previous_index = size(ref_only_x1_aligned2begin, 1);
                ref_only_x1_aligned2begin = [ref_only_x1_aligned2begin; NaN*ones(ref_only_num, tp_max)];
                ref_only_x2_aligned2begin = [ref_only_x2_aligned2begin; NaN*ones(ref_only_num, tp_max)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    profile_length = ROI_track_ref(ref_track, 1).tp_last - ROI_track_ref(ref_track, 1).tp_first + 1;
                    ref_only_x1_aligned2begin(previous_index + i,1:profile_length) = ROI_track_ref(ref_track,1).center_coord(:,1)';
                    ref_only_x2_aligned2begin(previous_index + i,1:profile_length) = ROI_track_ref(ref_track,1).center_coord(:,2)';
                end
                
                %% Coordinates: Align corresponding only tracks to the beginning of corresponding track
                previous_index = size(cor_only_x1_aligned2begin, 1);
                cor_only_x1_aligned2begin = [cor_only_x1_aligned2begin; NaN*ones(cor_only_num, tp_max)];
                cor_only_x2_aligned2begin = [cor_only_x2_aligned2begin; NaN*ones(cor_only_num, tp_max)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    profile_length = ROI_track_cor(cor_track, 1).tp_last - ROI_track_cor(cor_track, 1).tp_first + 1;
                    cor_only_x1_aligned2begin(previous_index + i,1:profile_length) = ROI_track_cor(cor_track,1).center_coord(:,1)';
                    cor_only_x2_aligned2begin(previous_index + i,1:profile_length) = ROI_track_cor(cor_track,1).center_coord(:,2)';
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Get statistics of reference tracks associated with corresponding tracks and reference only tracks: Lifetime and maximum intensity
                previous_index = size(lifetime_max_intensity_assoc_ref, 1);
                lifetime_max_intensity_assoc_ref= [lifetime_max_intensity_assoc_ref; NaN*ones(ref_cor_num, 2)];
                lifetime_mean_intensity_assoc_ref= [lifetime_mean_intensity_assoc_ref; NaN*ones(ref_cor_num, 2)];
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    lifetime_max_intensity_assoc_ref(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).ref_tp_last - ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval, ...
                        max(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_assoc_ref(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).ref_tp_last - ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval, ...
                        mean(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
                end
                
                %% Statistics for associated tracks (folder num 1,folder num 2, track ID, ref track lifetime, cor track lifetime
                %% profile length, ref track beginninng/ending time, cor track beginninng/ending time and profile beginning/ending time,
                %% max intensity, mean intenskty, mean change of intensity and mean change of position of reference track
                %% max intensity, mean intenskty, mean change of intensity and mean change of position of corresponding track
                previous_index = size(track_stat_asso, 1);
                track_stat_asso = [track_stat_asso; NaN*ones(ref_cor_num, 20)];
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    
                    track_stat_asso(previous_index+i,1) = folder_ID1;
                    track_stat_asso(previous_index+i,2) = folder_ID2;
                    track_stat_asso(previous_index+i,3) = ref_track; %track ID
                    track_stat_asso(previous_index+i,4) = (ROI_track_ref(ref_track, 1).ref_tp_last - ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval; % ref track lifetime
                    track_stat_asso(previous_index+i,5) = (ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval; % cor track lifetime
                    track_stat_asso(previous_index+i,6) = (ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all)*interval; % profile length
                    track_stat_asso(previous_index+i,7:8) = [ROI_track_ref(ref_track, 1).ref_tp_first, ROI_track_ref(ref_track, 1).ref_tp_last]; % ref track beginninng/ending time
                    track_stat_asso(previous_index+i,9:10) = [ROI_track_ref(ref_track, 1).cor_tp_first, ROI_track_ref(ref_track, 1).cor_tp_last]; % cor track beginninng/ending time
                    track_stat_asso(previous_index+i,11:12) = [ROI_track_ref(ref_track, 1).tp_prev_all, ROI_track_ref(ref_track, 1).tp_post_all]; % profile beginninng/ending time
                    track_stat_asso(previous_index+i,13) = max(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % ref max intensity
                    track_stat_asso(previous_index+i,14) = mean(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % ref mean intensity
                    track_stat_asso(previous_index+i,17) = max(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :)); % cor max intensity
                    track_stat_asso(previous_index+i,18) = mean(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :)); % cor mean intensity
                    % determine change in intensity of the reference track
                    track_index = (ROI_track_ref(ref_track, 1).ref_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all + 1):...
                        (ROI_track_ref(ref_track, 1).ref_tp_last - ROI_track_ref(ref_track, 1).tp_prev_all + 1);
                    track_stat_asso(previous_index+i,15) = mean((associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(2:end)) -...
                        associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_asso(previous_index+i,16) = mean(sum((ROI_track_ref(ref_track, 1).center_coord(2:end,:) -...
                        ROI_track_ref(ref_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    % determine change in intensity of the corresponding track
                    track_index = (ROI_track_ref(ref_track, 1).cor_tp_first - ROI_track_ref(ref_track, 1).tp_prev_all + 1):...
                        (ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).tp_prev_all + 1);
                    track_stat_asso(previous_index+i,19) = mean((associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, track_index(2:end)) -...
                        associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_asso(previous_index+i,20) = mean(sum((ROI_track_ref(ref_track, 1).center_coord(2:end,:) -...
                        ROI_track_ref(ref_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    
                end
                
                
                %% Statistics for ref-only tracks (folder num 1,folder num 2, track ID, track lifetime,
                %% profile length, track beginninng/ending time and profile beginning/ending time,
                %% max intensity, mean intenskty, mean change of intensity and mean change of position
                previous_index = size(lifetime_max_intensity_ref_only, 1);
                track_stat_ref_only = [track_stat_ref_only; NaN*ones(ref_only_num, 13)];
                lifetime_max_intensity_ref_only = [lifetime_max_intensity_ref_only; NaN*ones(ref_only_num, 2)];
                lifetime_mean_intensity_ref_only = [lifetime_mean_intensity_ref_only; NaN*ones(ref_only_num, 2)];
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    lifetime_max_intensity_ref_only(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).tp_last - ROI_track_ref(ref_track, 1).tp_first + 1)*interval, ...
                        max(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_ref_only(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).tp_last - ROI_track_ref(ref_track, 1).tp_first + 1)*interval, ...
                        mean(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
                    
                    track_stat_ref_only(previous_index+i,1) = folder_ID1;
                    track_stat_ref_only(previous_index+i,2) = folder_ID2;
                    track_stat_ref_only(previous_index+i,3) = ref_track; %track ID
                    track_stat_ref_only(previous_index+i,4) = (ROI_track_ref(ref_track, 1).tp_last - ROI_track_ref(ref_track, 1).tp_first + 1)*interval; % track lifetime
                    track_stat_ref_only(previous_index+i,5) = (ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all)*interval; % profile length
                    track_stat_ref_only(previous_index+i,6:7) = [ROI_track_ref(ref_track, 1).tp_first, ROI_track_ref(ref_track, 1).tp_last]; % track beginninng/ending time
                    track_stat_ref_only(previous_index+i,8:9) = [ROI_track_ref(ref_track, 1).tp_prev_all, ROI_track_ref(ref_track, 1).tp_post_all]; % profile beginninng/ending time
                    track_stat_ref_only(previous_index+i,10) = max(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % max intensity
                    track_stat_ref_only(previous_index+i,11) = mean(associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % mean intensity
                    % determine change in intensity
                    track_index = (ROI_track_ref(ref_track, 1).tp_first - ROI_track_ref(ref_track, 1).tp_prev_all + 1):...
                        (ROI_track_ref(ref_track, 1).tp_last - ROI_track_ref(ref_track, 1).tp_prev_all + 1);
                    track_stat_ref_only(previous_index+i,12) = mean((associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(2:end)) -...
                        associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_ref_only(previous_index+i,13) = mean(sum((ROI_track_ref(ref_track, 1).center_coord(2:end,:) -...
                        ROI_track_ref(ref_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    
                end
                
                
                %% Get statistics of corresponding tracks associated with reference tracks: Lifetime and maximum intensity
                previous_index = size(lifetime_max_intensity_assoc_cor, 1);
                lifetime_max_intensity_assoc_cor= [lifetime_mean_intensity_assoc_cor; NaN*ones(ref_cor_num, 2)];
                lifetime_mean_intensity_assoc_cor= [lifetime_mean_intensity_assoc_cor; NaN*ones(ref_cor_num, 2)];
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    lifetime_max_intensity_assoc_cor(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval, ...
                        max(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_assoc_cor(previous_index + i,:) = [(ROI_track_ref(ref_track, 1).cor_tp_last - ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval, ...
                        mean(associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :))];
                end
                
                %% Statistics for cor-only tracks (folder num 1,folder num 2, track ID, track lifetime,
                %% profile length, track beginninng/ending time and profile beginning/ending time,
                %% max intensity, mean intenskty, mean change of intensity and mean change of position
                previous_index = size(lifetime_max_intensity_cor_only, 1);
                track_stat_cor_only = [track_stat_cor_only; NaN*ones(cor_only_num, 13)];
                lifetime_mean_intensity_cor_only = [lifetime_mean_intensity_cor_only; NaN*ones(cor_only_num, 2)];
                lifetime_mean_intensity_cor_only = [lifetime_mean_intensity_cor_only; NaN*ones(cor_only_num, 2)];
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    lifetime_max_intensity_cor_only(previous_index + i,:) = [(ROI_track_cor(cor_track, 1).tp_last - ROI_track_cor(cor_track, 1).tp_first + 1)*interval, ...
                        max(unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_cor_only(previous_index + i,:) = [(ROI_track_cor(cor_track, 1).tp_last - ROI_track_cor(cor_track, 1).tp_first + 1)*interval, ...
                        mean(unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :))];
                    
                    track_stat_cor_only(previous_index+i,1) = folder_ID1;
                    track_stat_cor_only(previous_index+i,2) = folder_ID2;
                    track_stat_cor_only(previous_index+i,3) = cor_track; %track ID
                    track_stat_cor_only(previous_index+i,4) = (ROI_track_cor(cor_track, 1).tp_last - ROI_track_cor(cor_track, 1).tp_first + 1)*interval; % track lifetime
                    track_stat_cor_only(previous_index+i,5) = (ROI_track_cor(cor_track, 1).tp_post - ROI_track_cor(cor_track, 1).tp_prev)*interval; % profile length
                    track_stat_cor_only(previous_index+i,6:7) = [ROI_track_cor(cor_track, 1).tp_first, ROI_track_cor(cor_track, 1).tp_last]; % track beginninng/ending time
                    track_stat_cor_only(previous_index+i,8:9) = [ROI_track_cor(cor_track, 1).tp_prev, ROI_track_cor(cor_track, 1).tp_post]; % profile beginninng/ending time
                    track_stat_cor_only(previous_index+i,10) = max(unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :)); % max intensity
                    track_stat_cor_only(previous_index+i,11) = mean(unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :)); % mean intensity
                    % determine change in intensity
                    track_index = (ROI_track_cor(cor_track, 1).tp_first - ROI_track_cor(cor_track, 1).tp_prev + 1):...
                        (ROI_track_cor(cor_track, 1).tp_last - ROI_track_cor(cor_track, 1).tp_prev + 1);
                    track_stat_cor_only(previous_index+i,12) = mean((unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, track_index(2:end)) -...
                        unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_cor_only(previous_index+i,13) = mean(sum((ROI_track_cor(cor_track, 1).center_coord(2:end,:) -...
                        ROI_track_cor(cor_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    
                end
                % inventory of the indecies of "allAssociatedTracks" that
                % separate each dataset
                betweenDatasets(jjj) = jj-1;
                jjj = jjj+1;
            end
        end
    end
end
%% Determine ref/corcolor
if strcmp(analysis_info.ref_channel, 'RFP')
    analysis_info.ref_color = [1, 0, 0.8];
    analysis_info.cor_color = [0, 0.8, 0.1];
else
    analysis_info.ref_color = [0, 0.8, 0.1];
    analysis_info.cor_color = [1, 0, 0.8];
end
analysis_info.interval = interval;
%% Get screen size for figures.
scrsz = get(groot,'ScreenSize');
fig_width = 400;
fig_height = 250;


%% Save variables

if saveVariables
    
    fileLabel = input('What is the name of this experiment? ', 's');
    save(fileLabel)
    
    
end

%% Figure: Histogram of lifetime.

bins = linspace((tp_min-1)*interval, (tp_max-1)*interval, 13);
binedges = linspace((tp_min-1)*interval, (tp_max-1)*interval, 14);
bin_size = bins(2) - bins(1);
bin_center = (bins(1:end-1) + bins(2:end))/2;
bin_center = [bin_center, bin_center(end) + bin_size];

lifetime_hist_ref_only = histc(lifetime_ref_only, bins);
lifetime_hist_cor_only = histc(lifetime_cor_only, bins);
lifetime_hist_ref= histc(lifetime_associated(:,1), bins);
lifetime_hist_cor= histc(lifetime_associated(:,2), bins);

% transpose any data that's a row (if there are 1 or 0 values) so you can combine

if isrow(lifetime_hist_ref_only)
    lifetime_hist_ref_only = lifetime_hist_ref_only';
end

if isrow(lifetime_hist_cor_only)
    lifetime_hist_cor_only = lifetime_hist_cor_only';
end

if isrow(lifetime_hist_ref)
    lifetime_hist_ref = lifetime_hist_ref';
end

if isrow(lifetime_hist_cor)
    lifetime_hist_cor = lifetime_hist_cor';
end

lifetime_hist_data_all = [lifetime_hist_ref_only, lifetime_hist_cor_only,...
    lifetime_hist_ref, lifetime_hist_cor];

%figure('Position',[1 scrsz(4)-fig_height fig_width fig_height])
figure('OuterPosition',[1, scrsz(4)-fig_height, fig_width, fig_height])
bar(bin_center,lifetime_hist_data_all);
axis([-5, tp_max*interval, 0, max(max(lifetime_hist_data_all))+2]);
legend('refonly', 'coronly', 'ref', 'cor', 'Location', 'SouthEast')
set(gca,'XTick', bins)
bins_label = cell(1, numel(bins)+1);
for i = 1:(numel(bins))
    bins_label{i} = num2str(bins(i), '%4.0f');
end
bins_label{i+1} =  '+';
set(gca,'XTickLabel',bins_label)
text(tp_max*interval/2, max(max(lifetime_hist_data_all))*0.9,['ref lifetime=', num2str(nanmean(lifetime_associated(:,1)), '%4.2f'), '+-', num2str(nanstd(lifetime_associated(:,1)),'%4.2f')]);
text(tp_max*interval/2, max(max(lifetime_hist_data_all))*0.8,['cor lifetime=', num2str(nanmean(lifetime_associated(:,2)), '%4.2f'), '+-', num2str(nanstd(lifetime_associated(:,2)),'%4.2f')]);
text(tp_max*interval/2, max(max(lifetime_hist_data_all))*0.7,['ref only lifetime=', num2str(nanmean(lifetime_ref_only), '%4.2f'), '+-', num2str(nanstd(lifetime_ref_only),'%4.2f')]);
text(tp_max*interval/2, max(max(lifetime_hist_data_all))*0.6,['cor only lifetime=', num2str(nanmean(lifetime_cor_only), '%4.2f'), '+-', num2str(nanstd(lifetime_cor_only),'%4.2f')]);
xlabel('Lifetime (sec)')
ylabel('Number of incidents');
set(gcf,'color','w');
title('Histogram of lifetimes')

%% Histogram of lifetime, show only the associated tracks only.

figure('Position', [76   558   560   248]);
bar(bin_center,[lifetime_hist_ref,lifetime_hist_cor]);
axis([-5, tp_max*interval, 0, max(max([lifetime_hist_ref, lifetime_hist_cor]))+2]);
legend('ref','cor', 'Location','SouthEast')
set(gca,'XTick', bins)
bins_label = cell(1, numel(bins)+1);
for i = 1:(numel(bins))
    bins_label{i} = num2str(bins(i), '%4.0f');
end
bins_label{i+1} =  '+';
set(gca,'XTickLabel',bins_label)
text(tp_max*interval/2, max(max([lifetime_hist_cor, lifetime_hist_ref]))*0.9,['ref lifetime=', num2str(nanmean(lifetime_associated(:,1)), '%4.2f'), '+-', num2str(nanstd(lifetime_associated(:,1)),'%4.2f')]);
text(tp_max*interval/2, max(max([lifetime_hist_cor, lifetime_hist_ref]))*0.8,['cor lifetime=', num2str(nanmean(lifetime_associated(:,2)), '%4.2f'), '+-', num2str(nanstd(lifetime_associated(:,2)),'%4.2f')]);
xlabel('Lifetime (sec)')
ylabel('Number of incidents');
% title(strrep(folder_name, '_', '-'));
% suplabel(strrep(csv_file, '_', '-'));
set(gcf,'color','w');

%% Save histogram data.
%% Record the name of files
fid = fopen(histogram_file_name, 'w');
for foldernum1 = 1:folder_ID1
    for foldernum2 = 1:folder_ID2
        if ~isempty(foldername_all(foldernum1, foldernum2).associated_track_file)
            fprintf(fid, '%s\n', ['Track file(', num2str(foldernum1),',',num2str(foldernum2),'):',...
                fullfile(foldername_all(foldernum1, foldernum2).associated_track_folder, foldername_all(foldernum1, foldernum2).associated_track_file)]);
            fprintf(fid, '%s\n', ['Cleanup tag file(', num2str(foldernum1),',',num2str(foldernum2),'):',...
                fullfile(foldername_all(foldernum1, foldernum2).machine_tag_folder, foldername_all(foldernum1, foldernum2).machine_tag_file)]);
        end
    end
end
fprintf(fid, '%s\n', ['Interval btw time points:', num2str(interval), ' sec.']);
fprintf(fid, '%s\n\n', ['Number of time points in the time series:', num2str(tp_max)]);

% save reference lifetime
fprintf(fid, '%s\n', 'Lifetime Distribution');
fprintf(fid, strcat(repmat('%s\t',[1,6]),'%s\n'), 'Bin (min)', 'Bin (max)', 'Bin center', 'ref (associated)','cor (associated)',....
    'ref only', 'cor only');
fclose(fid);
dlmwrite(histogram_file_name, [bins', (bins+bin_size)', bin_center', lifetime_hist_ref, lifetime_hist_cor, lifetime_hist_ref_only,lifetime_hist_cor_only], '-append', 'delimiter', '\t');

%% Figure: Fraction of associated and unasociated tracks.
num_ref_only = sum(lifetime_hist_ref_only);
num_associated= sum(lifetime_hist_cor);
num_cor_only = sum(lifetime_hist_cor_only);

figure('OuterPosition', [fig_width+50, scrsz(4)-fig_height, fig_width,   fig_height]);
subplot(1,2,1);
bar([num_ref_only, num_associated]/(num_ref_only + num_associated))
set(gca,'XTickLabel',{'Unassociated','Associated'},'Fontsize', 8);
title(['reftrack, n=', num2str(num_ref_only + num_associated)])
subplot(1,2,2);
bar([num_cor_only, num_associated]/(num_cor_only + num_associated))
set(gca,'XTickLabel',{'Unassociated','Associated'}, 'Fontsize', 8);
title(['cortrack, n=', num2str(num_cor_only + num_associated)])
%% Plot to show the number of tracks shorter or longer than 20 sec
%{
    % reference tracks
    ref_w_cor_less20 = sum(lifetime_associated(:,1) < 20);
    ref_w_cor_larger20 = sum(lifetime_associated(:,1) >= 20);
    ref_only_less20 = sum(lifetime_ref_only < 20);
    ref_only_larger20 = sum(lifetime_ref_only >= 20);

    figure('Position', [69   584   560   202]);
    subplot(1,2,1);
    less_than_20 = ref_w_cor_less20/(ref_only_less20 + ref_w_cor_less20)*100;
    more_than_20 = ref_w_cor_larger20/(ref_only_larger20 + ref_w_cor_larger20)*100;
    bar([less_than_20,more_than_20]);
    text(0.8, 1,[num2str(less_than_20, '%2.1f'), '%'], 'Color', 'w');
    text(1.8, 1,[num2str(more_than_20, '%2.1f'), '%'], 'Color', 'w');
    set(gca,'XTickLabel',{'<20','>= 20'});
    title('reftracks with the corresponding tracks')
    subplot(1,2,2);
    less_than_20 = (ref_w_cor_less20 + ref_only_less20)/...
    (ref_w_cor_less20 + ref_only_less20 +ref_w_cor_larger20 + ref_only_larger20)*100;
    more_than_20 = (ref_w_cor_larger20 + ref_only_larger20)/...
    (ref_w_cor_less20 + ref_only_less20 +ref_w_cor_larger20 + ref_only_larger20)*100;
    bar([less_than_20,more_than_20]);
    text(0.8, 1,[num2str(less_than_20, '%2.1f'), '%'], 'Color', 'w');
    text(1.8, 1,[num2str(more_than_20, '%2.1f'), '%'], 'Color', 'w');
    set(gca,'XTickLabel',{'<20','>= 20'});
    title('% of reference tracks')

    % corresponding tracks
    cor_w_ref_less20 = sum(lifetime_associated(:,2) < 20);
    cor_w_ref_larger20 = sum(lifetime_associated(:,2) >= 20);
    cor_only_less20 = sum(lifetime_cor_only < 20);
    cor_only_larger20 = sum(lifetime_cor_only >= 20);

    figure('Position', [69   584   560   202])
    subplot(1,2,1);
    less_than_20 = cor_w_ref_less20/(cor_only_less20 + cor_w_ref_less20)*100;
    more_than_20 = cor_w_ref_larger20/(cor_only_larger20 + cor_w_ref_larger20)*100;
    bar([less_than_20,more_than_20]);
    text(0.8, 1,[num2str(less_than_20, '%2.1f'), '%'], 'Color', 'w');
    text(1.8, 1,[num2str(more_than_20, '%2.1f'), '%'], 'Color', 'w');
    set(gca,'XTickLabel',{'<20','>= 20'});
    title('cortracks with the reference tracks')
    subplot(1,2,2);
    less_than_20 = (cor_w_ref_less20 + cor_only_less20)/...
    (cor_w_ref_less20 + cor_only_less20 +cor_w_ref_larger20 + cor_only_larger20)*100;
    more_than_20 = (cor_w_ref_larger20 + cor_only_larger20)/...
    (cor_w_ref_less20 + cor_only_less20 +cor_w_ref_larger20 + cor_only_larger20)*100;
    bar([less_than_20,more_than_20]);
    text(0.8, 1,[num2str(less_than_20, '%2.1f'), '%'], 'Color', 'w');
    text(1.8, 1,[num2str(more_than_20, '%2.1f'), '%'], 'Color', 'w');
    set(gca,'XTickLabel',{'<20','>= 20'});
    title('% of corresponding tracks')
%}

%% plot all profiles
figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]);
for i=1:min(80, size(cor_aligned2cormax, 1))
    subplot(8,10,i);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_aligned2cormax(i, :),'Color', analysis_info.cor_color); hold on;
    subplot(8,10,i);plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_aligned2cormax(i, :),'Color',  analysis_info.ref_color);
    first_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'first');
    last_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'last');
    axis([(first_tp-tp_max)*2, (last_tp-tp_max)*2, min(min(ref_aligned2cormax)), max(max(ref_aligned2cormax))]);
end

% figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]);
% for i=1:min(80, size(cor_aligned2cormax, 1))
%     %scale cor I so its max is similar to max of ref
%     subplot(1,7,i);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_aligned2cormax(i, :)/nanmax(cor_aligned2cormax(:))*nanmax(ref_aligned2cormax(:)),'Color', analysis_info.cor_color); hold on;
%     subplot(1,7,i);plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_aligned2cormax(i, :),'Color',  analysis_info.ref_color);
%     first_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'first');
%     last_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'last');
%     axis([(first_tp-tp_max)*2, (last_tp-tp_max)*2, min(min(ref_aligned2cormax)), max(max(ref_aligned2cormax))]);
% end


%% plot all profiles (moving average with 3 tp)
%{
    figure;
    for i=1:min(80, size(cor_aligned2cormax, 1))
    first_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'first');
    last_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'last');
    subplot(8,10,i);plot((first_tp - tp_max)*interval:interval:(last_tp - tp_max)*interval, medfilt1(cor_aligned2cormax(i, first_tp:last_tp),3),'Color',  analysis_info.cor_color); hold on;
    subplot(8,10,i);plot((first_tp - tp_max)*interval:interval:(last_tp - tp_max)*interval, medfilt1(ref_aligned2cormax(i, first_tp:last_tp),3),'Color',  analysis_info.ref_color);
    axis([(first_tp-tp_max)*2, (last_tp-tp_max)*2, min(min(ref_aligned2cormax)), max(max(ref_aligned2cormax))]);
    end

    figure;
    subplot(1,2,1);plot(-(tp_max-1)*interval:interval:tp_max*interval, nanmean(cor_aligned2cormax));
    title('cor')
    axis([-tp_max*interval, tp_max*interval, min(min(cor_aligned2cormax)), max(max(cor_aligned2cormax))]);
    subplot(1,2,2);plot(-(tp_max-1)*interval:interval:tp_max*interval, nanmean(ref_aligned2cormax));
    title('ref')
    axis([-tp_max*interval, tp_max*interval, min(min(ref_aligned2cormax)), max(max(ref_aligned2cormax))]);
    suplabel('ref-cortracks. Aligned to maximum');
%}
if num_associated>0
    %% 1. Plot tracks of reference and corresponding tracks aligne to the max of corresponding track.
    % substitue NaNs with zero
    h = figure('OuterPosition', [1, scrsz(4)-fig_height*2, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_aligned2cormax, cor_aligned2cormax, analysis_info, 'associated', lifetime_associated, h);
    figure(h)
    suplabel(['Aligned to the max of cor. track, n=', num2str(N_sample)],'t');
    
    
    %% 2. Plot tracks of reference and corresponding tracks aligne to the max of reference track.
    h = figure('OuterPosition', [fig_width, scrsz(4)-fig_height*2, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_aligned2refmax, cor_aligned2refmax, analysis_info, 'associated', lifetime_associated, h);
    figure(h)
    suplabel(['Aligned to the max of ref. track, n=', num2str(N_sample)],'t');
    
    %% 3. Plot tracks of reference and corresponding tracks aligne to the beinning of reference track.
    h = figure('OuterPosition', [fig_width*2, scrsz(4)-fig_height*2, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_aligned2ref_begin, cor_aligned2ref_begin, analysis_info, 'associated', lifetime_associated, h);
    figure(h)
    suplabel(['Aligned to the beginning of ref. track, n=', num2str(N_sample)],'t');
    
    %% 3-2. Plot tracks of reference and corresponding tracks aligne to the beinning of reference track.
    %{
    figure('Position', [1140, 566, 560, 420]);
    subplot(3,2,1);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_aligned2ref_begin);
    axis([-tp_max*interval, tp_max*interval, min(min(cor_aligned2ref_begin)), max(max(cor_aligned2ref_begin))]);
    title('cor')
    subplot(3,2,2);plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_aligned2ref_begin);
    axis([-tp_max*interval, tp_max*interval, min(min(ref_aligned2ref_begin)), max(max(ref_aligned2ref_begin))]);
    title('ref')
    suplabel('ref-cortracks. Aligned to the beginning of reference track');

    % substitue NaNs with zero
    subplot(3,2,3:6)
    plot(-(tp_max*interval-interval):interval:tp_max*interval,...
        mean(cor_aligned2ref_begin_zero)/max(mean(cor_aligned2ref_begin_zero)), 'color',analysis_info.cor_color); hold on
    plot(-(tp_max*interval-interval):interval:tp_max*interval,...
        mean(ref_aligned2ref_begin_zero)/max(mean(ref_aligned2ref_begin_zero)), 'color', analysis_info.ref_color);

    axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval, -0.1, 1.1]);
    legend('cor', 'ref')
    xlabel('Time (sec)')
    ylabel('Normalized Average Intensity')
    suplabel('ref-cortracks. Aligned to the beginning of reference track');
    %}
    %% 4. Plot tracks of reference and corresponding tracks aligned to the ending of reference track.
    h = figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*2, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_aligned2ref_end, cor_aligned2ref_end, analysis_info, 'associated', lifetime_associated, h);
    figure(h)
    suplabel(['Aligned to the ending of ref. track, n=', num2str(N_sample)],'t');
    %% 4-2. Plot tracks of reference and corresponding tracks aligned to the ending of reference track.
    %{
    figure('Position', [14    71   560   420]);
    subplot(3,2,1);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_aligned2ref_end);
    axis([-tp_max*interval, tp_max*interval, min(min(cor_aligned2ref_end)), max(max(cor_aligned2ref_end))]);
    title('cor')
    subplot(3,2,2);plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_aligned2ref_end);
    axis([-tp_max*interval, tp_max*interval, min(min(ref_aligned2ref_end)), max(max(ref_aligned2ref_end))]);
    title('ref')
    suplabel('ref-cortracks. Aligned to to the ending of reference track');

    %
    % substitue NaNs with zero
    subplot(3,2,3:6)
    plot(-(tp_max*interval-interval):interval:tp_max*interval,...
        mean(cor_aligned2ref_end_zero)/max(mean(cor_aligned2ref_end_zero)), 'color',analysis_info.cor_color); hold on
    plot(-(tp_max*interval-interval):interval:tp_max*interval,...
        mean(ref_aligned2ref_end_zero)/max(mean(ref_aligned2ref_end_zero)), 'color', analysis_info.ref_color);

    axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval, -0.1, 1.1]);
    legend('cor', 'ref')
    xlabel('Time (sec)')
    ylabel('Normalized Average Intensity')
    suplabel(['ref-cortracks. Aligned to to the ending of reference track, n=', num2str(size(cor_aligned2ref_end,1))]);
    %}
    
    %% 5. Plot tracks of reference and corresponding tracks aligned to the beinning of corresponding track.
    h = figure('OuterPosition', [1, scrsz(4)-fig_height*3, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_aligned2cor_begin, cor_aligned2cor_begin, analysis_info, 'associated', lifetime_associated, h);
    figure(h)
    suplabel(['Aligned to the beinning of cor. track, n=', num2str(N_sample)],'t');
    
    %% 6. Plot tracks of reference and corresponding tracks aligned to the disappearance of corresponding track.
    h = figure('OuterPosition', [fig_width, scrsz(4)-fig_height*3, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_aligned2cor_end, cor_aligned2cor_end, analysis_info, 'associated', lifetime_associated, h);
    figure(h)
    suplabel(['Aligned to the disappearance of cor. track, n=', num2str(N_sample)],'t');
    
    %% 6-2. Plot tracks of reference and corresponding tracks aligned to the disappearance of corresponding track.
    %{
    figure('Position', [1136,89, 560, 420]);
    subplot(3,2,1);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_aligned2cor_end);
    axis([-tp_max*interval, tp_max*interval, min(min(cor_aligned2cor_end)), max(max(cor_aligned2cor_end))]);
    title('cor')
    subplot(3,2,2);plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_aligned2cor_end);
    axis([-tp_max*interval, tp_max*interval, min(min(ref_aligned2cor_end)), max(max(ref_aligned2cor_end))]);
    title('ref')
    suplabel('ref-cortracks. Aligned to the disappearance of corresponding track');

    % substitue NaNs with zero
    subplot(3,2,3:6)
    plot(-(tp_max*interval-interval):interval:tp_max*interval,...
        mean(cor_aligned2cor_end_zero)/max(mean(cor_aligned2cor_end_zero)), 'color',analysis_info.cor_color); hold on
    plot(-(tp_max*interval-interval):interval:tp_max*interval,...
        mean(ref_aligned2cor_end_zero)/max(mean(ref_aligned2cor_end_zero)), 'color', analysis_info.ref_color);


    axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval, -0.1, 1.1]);
    legend('cor', 'ref')
    xlabel('Time (sec)')
    ylabel('Normalized Average Intensity')
    suplabel(['ref-cortracks. Aligned to the disappearance of corresponding track, n=', num2str(size(cor_aligned2cor_end,1))]);
    %}
end
%% Plot unassociated reference tracks
if num_ref_only>0
    h = figure('OuterPosition', [fig_width*2, scrsz(4)-fig_height*3, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_only_ref_aligned2ref_begin, ref_only_cor_aligned2ref_begin, analysis_info, 'ref', lifetime_ref_only, h);
    figure(h)
    suplabel(['Ref only. Aligned to the appearance of ref. track, n=', num2str(N_sample)],'t');
    
    h = figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*3, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_only_ref_aligned2ref_end, ref_only_cor_aligned2ref_end, analysis_info, 'ref', lifetime_ref_only, h);
    figure(h)
    suplabel(['Ref only. Aligned to the disappearance of ref. track, n=', num2str(N_sample)],'t');
    
    h = figure('OuterPosition', [1, scrsz(4)-fig_height*4, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(ref_only_ref_aligned2refmax, ref_only_cor_aligned2refmax, analysis_info, 'ref', lifetime_ref_only, h);
    figure(h)
    suplabel(['Ref only. Aligned to the maximum of ref. track, n=', num2str(N_sample)],'t');
end
%% Plot unassociated corresponding tracks
if num_cor_only>0
    h = figure('OuterPosition', [fig_width, scrsz(4)-fig_height*4, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(cor_only_ref_aligned2cor_begin, cor_only_cor_aligned2cor_begin, analysis_info, 'cor', lifetime_cor_only, h);
    figure(h)
    suplabel(['Cor only. Aligned to the appearance of cor. track, n=', num2str(N_sample)],'t');
    
    h = figure('OuterPosition', [fig_width*2, scrsz(4)-fig_height*4, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(cor_only_ref_aligned2cor_end, cor_only_cor_aligned2cor_end, analysis_info, 'cor', lifetime_cor_only, h);
    figure(h)
    suplabel(['Cor only. Aligned to the disappearance of cor. track, n=', num2str(N_sample)],'t');
    
    h = figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*4, fig_width,   fig_height]);
    [h, N_sample] = drawAveragePlots(cor_only_ref_aligned2cormax, cor_only_cor_aligned2cormax, analysis_info, 'cor', lifetime_cor_only, h);
    figure(h)
    suplabel(['Cor only. Aligned to the maximum of cor. track, n=', num2str(N_sample)],'t');
end
%% Plot statistics of reference tracks associated with corresponding tracks and reference only tracks.
if or(size(lifetime_max_intensity_assoc_ref,1)>0,size(lifetime_max_intensity_ref_only, 1)>0)
    figure('OuterPosition', [fig_width*2, scrsz(4)-fig_height*3, fig_height,   fig_height]);
    if size(lifetime_max_intensity_assoc_ref,1)>0
        plot(lifetime_max_intensity_assoc_ref(:,1), lifetime_max_intensity_assoc_ref(:,2), 'bo'); hold on
    end
    if size(lifetime_max_intensity_ref_only, 1)>0
        plot(lifetime_max_intensity_ref_only(:,1), lifetime_max_intensity_ref_only(:,2), 'ro');
    end
    xlabel('Track lifetime (Sec)')
    ylabel('Max intensity')
    title('Reference track stat. Blue: reference associated with cor. Red: reference only.')
end
%{
    figure('Position', [99   168   333   260]);
    plot(lifetime_mean_intensity_assoc_ref(:,1), lifetime_mean_intensity_assoc_ref(:,2), 'bo'); hold on
    plot(lifetime_mean_intensity_ref_only(:,1), lifetime_mean_intensity_ref_only(:,2), 'ro');
    xlabel('Track lifetime (Sec)')
    ylabel('Mean intensity')
    title('refstat. Blue: reference associated with cor. Red: reference only.')
<<<<<<< HEAD
%}
%% Plot statistics of corresponding tracks associated with reference tracks and corresponding only tracks.
if or(size(lifetime_max_intensity_assoc_cor,1)>0, size(lifetime_max_intensity_cor_only, 1)>0)
    figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*3, fig_height,   fig_height]);
    if size(lifetime_max_intensity_assoc_cor,1)>0
        plot(lifetime_max_intensity_assoc_cor(:,1), lifetime_max_intensity_assoc_cor(:,2), 'bo'); hold on
    end
    if size(lifetime_max_intensity_cor_only, 1)>0
        plot(lifetime_max_intensity_cor_only(:,1), lifetime_max_intensity_cor_only(:,2), 'ro');
    end
    xlabel('Track lifetime (Sec)')
    ylabel('Max intensity')
    title('Corresponding track stat. Blue: corresponding associated with ref. Red: corresponding only.')
end

%% Plot individual "unusual" tracks
%% CASEY look here! copy and paste for ref_cor_4_num, etc

if plot_2and3_tracks
    
    
    ii = 1;
    
    for i = 1:ref_cor_2_num
        ref_track = ref_cor_tracks_2(i);
        
        curCorTrack = (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
        curRefTrack = (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
        
        unusualAssocTracks(ii).Cor = curCorTrack;
        unusualAssocTracks(ii).Ref = curRefTrack;
        
        sizes_ref_cor_tracks_2(ii) = length(curCorTrack);
        
        ii = ii + 1;
    end
    
    % Ref only
    ii = 1;
    
    for i = 1:ref_only_2_num
        ref_track = ref_only_tracks_2(i);
        
        curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
        curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
        
        
        unusualRefOnlyTracks(ii).Ref = curRefTrack;
        
        % Note that "RefOnlyTracks.Cor" is the intensity trace in the Cor channel
        % that corresponds to the track in the Ref channel, even though
        % this track is unassociated. So in normal cases, when the ref
        % track really is unassociated,
        % "RefOnlyTracks.Cor" should have little fluorescence
        % intensity.
        
        unusualRefOnlyTracks(ii).Cor = curCorTrack;
        sizes_ref_only_tracks_2(ii) = length(curCorTrack);
        
        
        ii = ii + 1;
    end
    
    
    % Cor only
    ii = 1;
    
    for i = 1:cor_only_2_num
        ref_track = cor_only_tracks_2(i);
        
        curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
        curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
        
        unusualCorOnlyTracks(ii).Cor = curCorTrack;
        unusualCorOnlyTracks(ii).Ref = curRefTrack;
        sizes_cor_only_tracks_2(ii) = length(curCorTrack);
        
        ii = ii + 1;
        %}
        %% Plot statistics of corresponding tracks associated with reference tracks and corresponding only tracks.
        if or(size(lifetime_max_intensity_assoc_cor,1)>0, size(lifetime_max_intensity_cor_only, 1)>0)
            figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*3, fig_height,   fig_height]);
            if size(lifetime_max_intensity_assoc_cor,1)>0
                plot(lifetime_max_intensity_assoc_cor(:,1), lifetime_max_intensity_assoc_cor(:,2), 'bo'); hold on
            end
            if size(lifetime_max_intensity_cor_only, 1)>0
                plot(lifetime_max_intensity_cor_only(:,1), lifetime_max_intensity_cor_only(:,2), 'ro');
            end
            xlabel('Track lifetime (Sec)')
            ylabel('Max intensity')
            title('Corresponding track stat. Blue: corresponding associated with ref. Red: corresponding only.')
        end
        
    end
        %% Plot individual "unusual" tracks
        
        if plot_2and3_tracks
            
            
            ii = 1;
            
            for i = 1:ref_cor_2_num
                ref_track = ref_cor_tracks_2(i);
                
                curCorTrack = (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                curRefTrack = (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                
                unusualAssocTracks(ii).Cor = curCorTrack;
                unusualAssocTracks(ii).Ref = curRefTrack;
                
                sizes_ref_cor_tracks_2(ii) = length(curCorTrack);
                
                ii = ii + 1;
            end
            
            % Ref only
            ii = 1;
            
            for i = 1:ref_only_2_num
                ref_track = ref_only_tracks_2(i);
                
                curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                
                
                unusualRefOnlyTracks(ii).Ref = curRefTrack;
                
                % Note that "RefOnlyTracks.Cor" is the intensity trace in the Cor channel
                % that corresponds to the track in the Ref channel, even though
                % this track is unassociated. So in normal cases, when the ref
                % track really is unassociated,
                % "RefOnlyTracks.Cor" should have little fluorescence
                % intensity.
                
                unusualRefOnlyTracks(ii).Cor = curCorTrack;
                sizes_ref_only_tracks_2(ii) = length(curCorTrack);
                
                
                ii = ii + 1;
            end
            
            
            % Cor only
            ii = 1;
            
            for i = 1:cor_only_2_num
                ref_track = cor_only_tracks_2(i);
                
                curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                
                unusualCorOnlyTracks(ii).Cor = curCorTrack;
                unusualCorOnlyTracks(ii).Ref = curRefTrack;
                sizes_cor_only_tracks_2(ii) = length(curCorTrack);
                
                ii = ii + 1;
            end
            
            
            % Plot
            if ref_cor_2_num>0
                
                figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
                
                nbSubplotRows = ceil(ref_cor_2_num/8);
                
                subplot(ceil(ref_cor_2_num/nbSubplotRows), nbSubplotRows, 1);
                title('Unusual tracks associated'); hold on;
                
                
                for ii = 1:ref_cor_2_num
                    
                    
                    subplot(ceil(ref_cor_2_num/nbSubplotRows), nbSubplotRows, ii);
                    
                    plot(interval:interval:length(unusualAssocTracks(ii).Cor)*interval, unusualAssocTracks(ii).Cor,  'g'); hold on;
                    plot(interval:interval:length(unusualAssocTracks(ii).Ref)*interval, unusualAssocTracks(ii).Ref,  'm'); hold on;
                    
                    xlim([0 max(sizes_ref_cor_tracks_2)*interval+interval]);
                    
                end
                xlabel('Time (s)');
            end
            
            
            if ref_only_2_num>0
                
                figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
                
                nbSubplotRows = ceil(ref_only_2_num/8);
                
                subplot(ceil(ref_only_2_num/nbSubplotRows), nbSubplotRows, 1);
                
                title('Unusual tracks ref only'); hold on;
                
                for ii = 1:ref_only_2_num
                    
                    subplot(ceil(ref_only_2_num/nbSubplotRows), nbSubplotRows, ii);
                    
                    plot(interval:interval:length(unusualRefOnlyTracks(ii).Cor)*interval, unusualRefOnlyTracks(ii).Cor, 'g'); hold on;
                    plot(interval:interval:length(unusualRefOnlyTracks(ii).Cor)*interval, unusualRefOnlyTracks(ii).Ref, 'm'); hold on;
                    xlim([0 max(sizes_ref_only_tracks_2)*interval+interval]);
                    
                end
                xlabel('Time (s)');
                
            end
            
            if cor_only_2_num>0
                
                figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
                
                nbSubplotRows = ceil(cor_only_2_num/8);
                
                subplot(ceil(cor_only_2_num/nbSubplotRows), nbSubplotRows, 1);
                title('Unusual tracks cor only'); hold on;
                
                for ii = 1:cor_only_2_num
                    
                    
                    subplot(ceil(cor_only_2_num/nbSubplotRows), nbSubplotRows, ii);
                    
                    plot(interval:interval:length(unusualCorOnlyTracks(ii).Cor)*interval, unusualCorOnlyTracks(ii).Cor, 'g'); hold on;
                    plot(interval:interval:length(unusualCorOnlyTracks(ii).Cor)*interval, unusualCorOnlyTracks(ii).Ref, 'm'); hold on;
                    xlim([0 max(sizes_cor_only_tracks_2)*interval+interval]);
                    
                end
                
                xlabel('Time (s)')
            end
        end
        
        
        %% Plot individual "special case" tracks
        
        if plot_2and3_tracks
            
            
            
            ii = 1;
            
            for i = 1:ref_cor_3_num
                ref_track = ref_cor_tracks_3(i);
                
                curCorTrack = (associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                curRefTrack = (associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                
                specialCaseAssocTracks(ii).Cor = curCorTrack;
                specialCaseAssocTracks(ii).Ref = curRefTrack;
                
                sizes_ref_cor_tracks_3(ii) = length(curCorTrack);
                
                ii = ii + 1;
            end
            
            % Ref only
            ii = 1;
            
            for i = 1:ref_only_3_num
                ref_track = ref_only_tracks_3(i);
                
                curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                
                
                specialCaseRefOnlyTracks(ii).Ref = curRefTrack;
                
                % Note that "RefOnlyTracks.Cor" is the intensity trace in the Cor channel
                % that corresponds to the track in the Ref channel, even though
                % this track is unassociated. So in normal cases, when the ref
                % track really is unassociated,
                % "RefOnlyTracks.Cor" should have little fluorescence
                % intensity.
                
                specialCaseRefOnlyTracks(ii).Cor = curCorTrack;
                sizes_ref_only_tracks_3(ii) = length(curCorTrack);
                
                
                ii = ii + 1;
            end
            
            
            % Cor only
            ii = 1;
            
            for i = 1:cor_only_3_num
                ref_track = cor_only_tracks_3(i);
                
                curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                
                specialCaseCorOnlyTracks(ii).Cor = curCorTrack;
                specialCaseCorOnlyTracks(ii).Ref = curRefTrack;
                
                sizes_cor_only_tracks_3(ii) = length(curCorTrack);
                
                ii = ii + 1;
            end
            
            
            
            
            
            % Plot
            
            if ref_cor_3_num>0
                figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
                
                
                nbSubplotRows = ceil(ref_cor_3_num/8);
                
                subplot(ceil(ref_cor_3_num/nbSubplotRows), nbSubplotRows, 1);
                
                title('Special case tracks associated'); hold on;
                
                for ii = 1:ref_cor_3_num
                    
                    subplot(ceil(ref_cor_3_num/nbSubplotRows), nbSubplotRows, ii);
                    
                    plot(interval:interval:length(specialCaseAssocTracks(ii).Cor)*interval, specialCaseAssocTracks(ii).Cor, 'g'); hold on;
                    plot(interval:interval:length(specialCaseAssocTracks(ii).Cor)*interval, specialCaseAssocTracks(ii).Ref, 'm'); hold on;
                    xlim([0 max(sizes_ref_cor_tracks_3)*interval+interval]);
                    
                end
                xlabel('Time (s)');
            end
            
            if ref_only_3_num>0
                figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
                
                
                nbSubplotRows = ceil(ref_only_3_num/8);
                
                subplot(ceil(ref_only_3_num/nbSubplotRows), nbSubplotRows, 1);
                
                title('Special case tracks ref only'); hold on;
                
                for ii = 1:ref_only_3_num
                    
                    
                    subplot(ceil(ref_only_3_num/nbSubplotRows), nbSubplotRows, ii);
                    
                    plot(interval:interval:length(specialCaseRefOnlyTracks(ii).Cor)*interval, specialCaseRefOnlyTracks(ii).Cor, 'g'); hold on;
                    plot(interval:interval:length(specialCaseRefOnlyTracks(ii).Cor)*interval, specialCaseRefOnlyTracks(ii).Ref, 'm'); hold on;
                    xlim([0 max(sizes_ref_only_tracks_3)*interval+interval]);
                    
                end
                xlabel('Time (s)');
                
            end
            
            if cor_only_3_num>0
                figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
                
                nbSubplotRows = ceil(cor_only_3_num/8);
                subplot(ceil(cor_only_3_num/nbSubplotRows), nbSubplotRows, 1);
                title('Special case tracks cor only'); hold on;
                
                for ii = 1:cor_only_3_num
                    
                    
                    subplot(ceil(cor_only_3_num/nbSubplotRows), nbSubplotRows, ii);
                    
                    plot(interval:interval:length(specialCaseCorOnlyTracks(ii).Cor)*interval, specialCaseCorOnlyTracks(ii).Cor, 'g'); hold on;
                    plot(interval:interval:length(specialCaseCorOnlyTracks(ii).Cor)*interval, specialCaseCorOnlyTracks(ii).Ref, 'm'); hold on;
                    xlim([0 max(sizes_cor_only_tracks_3)*interval+interval]);
                    
                end
                xlabel('Time (s)');
                
            end
        end
    end
    
    %% Plot Timepoint violators High
    
    if plot_4hilo
        
        % MA this seems broken.
        ii = 1;
        
        for i = 1:ref_cor_tracks_4_high_num
            ref_track = ref_cor_tracks_4_high(i);
            
            curCorTrack = (associated_tracks(ref_track, 1).cor_intensity(1, :));
            curRefTrack = (associated_tracks(ref_track, 1).ref_intensity(1, :));
            
            highAssocTracks(ii).Cor = curCorTrack;
            highAssocTracks(ii).Ref = curRefTrack;
            
            sizes_ref_cor_tracks_4_high(ii) = length(curCorTrack);
            
            ii = ii + 1;
        end
        
        % Ref only
        ii = 1;
        
        for i = 1:ref_only_tracks_4_high_num
            ref_track = ref_only_tracks_4_high(i);
            
            curCorTrack = track_stat_ref.old_unassoc_ref_intensity(ref_track).a(1, :);
            curRefTrack = track_stat_ref.old_unassoc_cor_intensity(ref_track).a(1, :);
            
            
            highRefOnlyTracks(ii).Ref = curRefTrack;
            
            % Note that "RefOnlyTracks.Cor" is the intensity trace in the Cor channel
            % that corresponds to the track in the Ref channel, even though
            % this track is unassociated. So in normal cases, when the ref
            % track really is unassociated,
            % "RefOnlyTracks.Cor" should have little fluorescence
            % intensity.
            
            highRefOnlyTracks(ii).Cor = curCorTrack;
            sizes_ref_only_tracks_4_high(ii) = length(curCorTrack);
            
            
            ii = ii + 1;
        end
        
        
        % Cor only
        ii = 1;
        
        for i = 1:cor_only_tracks_4_high_num
            ref_track = cor_only_tracks_4_high(i);
            
            curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity(1, :));
            curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity(1, :));
            
            highCorOnlyTracks(ii).Cor = curCorTrack;
            highCorOnlyTracks(ii).Ref = curRefTrack;
            sizes_cor_only_tracks_4_high(ii) = length(curCorTrack);
            
            ii = ii + 1;
        end
        
        
        % Plot
        if ref_cor_tracks_4_high_num>0
            
            figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
            
            nbSubplotRows = ceil(ref_cor_tracks_4_high_num/8);
            
            subplot(ceil(ref_cor_tracks_4_high_num/nbSubplotRows), nbSubplotRows, 1);
            title('High tracks associated'); hold on;
            
            
            for ii = 1:ref_cor_tracks_4_high_num
                
                
                subplot(ceil(ref_cor_tracks_4_high_num/nbSubplotRows), nbSubplotRows, ii);
                
                plot(interval:interval:length(highAssocTracks(ii).Cor)*interval, highAssocTracks(ii).Cor,  'g'); hold on;
                plot(interval:interval:length(highAssocTracks(ii).Ref)*interval, highAssocTracks(ii).Ref,  'm'); hold on;
                
                xlim([0 max(sizes_ref_cor_tracks_4_high)*interval+interval]);
                
            end
            xlabel('Time (s)');
        end
        
        
        if ref_only_tracks_4_high_num>0
            
            figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
            
            nbSubplotRows = ceil(ref_only_tracks_4_high_num/8);
            
            subplot(ceil(ref_only_tracks_4_high_num/nbSubplotRows), nbSubplotRows, 1);
            
            title('High tracks ref only'); hold on;
            
            for ii = 1:ref_only_tracks_4_high_num
                
                subplot(ceil(ref_only_tracks_4_high_num/nbSubplotRows), nbSubplotRows, ii);
                
                plot(interval:interval:length(highRefOnlyTracks(ii).Cor)*interval, highRefOnlyTracks(ii).Cor, 'g'); hold on;
                plot(interval:interval:length(highRefOnlyTracks(ii).Cor)*interval, highRefOnlyTracks(ii).Ref, 'm'); hold on;
                xlim([0 max(sizes_ref_only_tracks_4_high)*interval+interval]);
                
            end
            xlabel('Time (s)');
            
        end
        
        if cor_only_tracks_4_high_num>0
            
            figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
            
            nbSubplotRows = ceil(cor_only_tracks_4_high_num/8);
            
            subplot(ceil(cor_only_tracks_4_high_num/nbSubplotRows), nbSubplotRows, 1);
            title('High tracks cor only'); hold on;
            
            for ii = 1:cor_only_tracks_4_high_num
                
                
                subplot(ceil(cor_only_tracks_4_high_num/nbSubplotRows), nbSubplotRows, ii);
                
                plot(interval:interval:length(highCorOnlyTracks(ii).Cor)*interval, highCorOnlyTracks(ii).Cor, 'g'); hold on;
                plot(interval:interval:length(highCorOnlyTracks(ii).Cor)*interval, highCorOnlyTracks(ii).Ref, 'm'); hold on;
                xlim([0 max(sizes_cor_only_tracks_4_high)*interval+interval]);
                
            end
            
            xlabel('Time (s)')
        end
    end
    
    
    %% Plot Timepoint Violators Low
    
    
    if plot_4hilo
        
        
        ii = 1;
        
        for i = 1:ref_cor_tracks_4_low_num
            ref_track = ref_cor_tracks_4_low(i);
            
            curCorTrack = (associated_tracks(ref_track, 1).cor_intensity(1, :));
            curRefTrack = (associated_tracks(ref_track, 1).ref_intensity(1, :));
            
            lowAssocTracks(ii).Cor = curCorTrack;
            lowAssocTracks(ii).Ref = curRefTrack;
            
            sizes_ref_cor_tracks_4_low(ii) = length(curCorTrack);
            
            ii = ii + 1;
        end
        
        % Ref only
        ii = 1;
        
        for i = 1:ref_only_tracks_4_low_num
            ref_track = ref_only_tracks_4_low(i);
            
            curCorTrack = track_stat_ref.old_unassoc_ref_intensity(ref_track).a(1, :);
            curRefTrack = track_stat_ref.old_unassoc_cor_intensity(ref_track).a(1, :);
            
            
            lowRefOnlyTracks(ii).Ref = curRefTrack;
            
            % Note that "RefOnlyTracks.Cor" is the intensity trace in the Cor channel
            % that corresponds to the track in the Ref channel, even though
            % this track is unassociated. So in normal cases, when the ref
            % track really is unassociated,
            % "RefOnlyTracks.Cor" should have little fluorescence
            % intensity.
            
            lowRefOnlyTracks(ii).Cor = curCorTrack;
            sizes_ref_only_tracks_4_low(ii) = length(curCorTrack);
            
            
            ii = ii + 1;
        end
        
        
        % Cor only
        ii = 1;
        
        for i = 1:cor_only_tracks_4_low_num
            ref_track = cor_only_tracks_4_low(i);
            
            curCorTrack = (unassociated_tracks(ref_track, 1).cor_intensity(1, :));
            curRefTrack = (unassociated_tracks(ref_track, 1).ref_intensity(1, :));
            
            lowCorOnlyTracks(ii).Cor = curCorTrack;
            lowCorOnlyTracks(ii).Ref = curRefTrack;
            sizes_cor_only_tracks_4_low(ii) = length(curCorTrack);
            
            ii = ii + 1;
        end
        
        
        % Plot
        if ref_cor_tracks_4_low_num>0
            
            figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
            
            nbSubplotRows = ceil(ref_cor_tracks_4_low_num/8);
            
            subplot(ceil(ref_cor_tracks_4_low_num/nbSubplotRows), nbSubplotRows, 1);
            title('low tracks associated'); hold on;
            
            
            for ii = 1:ref_cor_tracks_4_low_num
                
                
                subplot(ceil(ref_cor_tracks_4_low_num/nbSubplotRows), nbSubplotRows, ii);
                
                plot(interval:interval:length(lowAssocTracks(ii).Cor)*interval, lowAssocTracks(ii).Cor,  'g'); hold on;
                plot(interval:interval:length(lowAssocTracks(ii).Ref)*interval, lowAssocTracks(ii).Ref,  'm'); hold on;
                
                xlim([0 max(sizes_ref_cor_tracks_4_low)*interval+interval]);
                
            end
            xlabel('Time (s)');
        end
        
        
        if ref_only_tracks_4_low_num>0
            
            figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
            
            nbSubplotRows = ceil(ref_only_tracks_4_low_num/8);
            
            subplot(ceil(ref_only_tracks_4_low_num/nbSubplotRows), nbSubplotRows, 1);
            
            title('low tracks ref only'); hold on;
            
            for ii = 1:ref_only_tracks_4_low_num
                
                subplot(ceil(ref_only_tracks_4_low_num/nbSubplotRows), nbSubplotRows, ii);
                
                plot(interval:interval:length(lowRefOnlyTracks(ii).Cor)*interval, lowRefOnlyTracks(ii).Cor, 'g'); hold on;
                plot(interval:interval:length(lowRefOnlyTracks(ii).Cor)*interval, lowRefOnlyTracks(ii).Ref, 'm'); hold on;
                xlim([0 max(sizes_ref_only_tracks_4_low)*interval+interval]);
                
            end
            xlabel('Time (s)');
            
        end
        
        if cor_only_tracks_4_low_num>0
            
            figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]); hold on;
            
            nbSubplotRows = ceil(cor_only_tracks_4_low_num/8);
            
            subplot(ceil(cor_only_tracks_4_low_num/nbSubplotRows), nbSubplotRows, 1);
            title('low tracks cor only'); hold on;
            
            for ii = 1:cor_only_tracks_4_low_num
                
                
                subplot(ceil(cor_only_tracks_4_low_num/nbSubplotRows), nbSubplotRows, ii);
                
                plot(interval:interval:length(lowCorOnlyTracks(ii).Cor)*interval, lowCorOnlyTracks(ii).Cor, 'g'); hold on;
                plot(interval:interval:length(lowCorOnlyTracks(ii).Cor)*interval, lowCorOnlyTracks(ii).Ref, 'm'); hold on;
                xlim([0 max(sizes_cor_only_tracks_4_low)*interval+interval]);
                
            end
            
            xlabel('Time (s)')
        end
    end
    
    
    %%
    %{    
    
    %{
>>>>>>> Casey_stats_branch
    figure('Position', [440   168   333   260])
    plot(lifetime_mean_intensity_assoc_cor(:,1), lifetime_mean_intensity_assoc_cor(:,2), 'bo'); hold on
    plot(lifetime_mean_intensity_cor_only(:,1), lifetime_mean_intensity_cor_only(:,2), 'ro');
    xlabel('Track lifetime (Sec)')
    ylabel('Mean intensity')
    title('corstat. Blue: corresponding associated with ref. Red: corresponding only.')
    %}
    %% Plot statistics of corresponding only and reference only tracks
    %{
    figure('Position', [840   168   333   260]);
    plot(track_stat_cor_only(:,4), track_stat_cor_only(:,12), 'o', 'color', analysis_info.cor_color); hold on
    plot(track_stat_ref_only(:,4), track_stat_ref_only(:,12), 'o', 'color', analysis_info.ref_color);
    xlabel('Track lifetime (Sec)')
    ylabel('mean change of intensity')

    figure('Position', [840   540   333   260]);
    plot(track_stat_cor_only(:,4), track_stat_cor_only(:,13), 'o', 'color', analysis_info.cor_color); hold on
    plot(track_stat_ref_only(:,4), track_stat_ref_only(:,13), 'o', 'color', analysis_info.ref_color);
    xlabel('Track lifetime (Sec)')
    ylabel('mean change of position')
    %}
    %% Plot mean squared displacement and distance traveled for corresponding and reference tracks
    %{
    ref_only_MSD = (ref_only_x1_aligned2begin(:, 2:end) - repmat(ref_only_x1_aligned2begin(:,1), [1, tp_max - 1])).^2 +...
    (ref_only_x2_aligned2begin(:, 2:end) - repmat(ref_only_x2_aligned2begin(:,1), [1, tp_max - 1])).^2 ;
    cor_only_MSD = (cor_only_x1_aligned2begin(:, 2:end) - repmat(cor_only_x1_aligned2begin(:,1), [1, tp_max - 1])).^2 +...
    (cor_only_x2_aligned2begin(:, 2:end) - repmat(cor_only_x2_aligned2begin(:,1), [1, tp_max - 1])).^2 ;

    h = figure('Position', [1240   540   333   260]);
    plot([1:(tp_max-1)]*interval, nanmean(ref_only_MSD, 1), 'color', analysis_info.ref_color); hold on
    plot([1:(tp_max-1)]*interval, nanmean(cor_only_MSD, 1), 'color', analysis_info.cor_color); hold on
    title('MSD')
    %
    ref_only_dist = cumsum(sqrt((ref_only_x1_aligned2begin(:, 2:end) - ref_only_x1_aligned2begin(:,1:(end-1))).^2 +...
    (ref_only_x2_aligned2begin(:, 2:end) - ref_only_x2_aligned2begin(:,1:(end-1))).^2), 2);
    cor_only_dist = cumsum(sqrt((cor_only_x1_aligned2begin(:, 2:end) - cor_only_x1_aligned2begin(:,1:(end-1))).^2 +...
    (cor_only_x2_aligned2begin(:, 2:end) - cor_only_x2_aligned2begin(:,1:(end-1))).^2), 2);
    figure('Position', [1240   168   333   260]);
    plot([1:(tp_max-1)]*interval, nanmean(ref_only_dist, 1), 'color', analysis_info.ref_color); hold on
    plot([1:(tp_max-1)]*interval, nanmean(cor_only_dist, 1), 'color', analysis_info.cor_color); hold on
    title('distance traveled')
    %}
    %% Plot mean squared displacement and distance traveled for associated corresponding and reference tracks
    %{
    ref_MSD = (ref_x1_aligned2ref_begin(:, 2:end) - repmat(ref_x1_aligned2ref_begin(:,1), [1, tp_max - 1])).^2 +...
    (ref_x2_aligned2ref_begin(:, 2:end) - repmat(ref_x2_aligned2ref_begin(:,1), [1, tp_max - 1])).^2 ;
    cor_MSD = (cor_x1_aligned2cor_begin(:, 2:end) - repmat(cor_x1_aligned2cor_begin(:,1), [1, tp_max - 1])).^2 +...
    (cor_x2_aligned2cor_begin(:, 2:end) - repmat(cor_x2_aligned2cor_begin(:,1), [1, tp_max - 1])).^2 ;

    figure(h)
    plot([1:(tp_max-1)]*interval, nanmean(ref_MSD, 1), 'color', analysis_info.ref_color*.75); hold on
    plot([1:(tp_max-1)]*interval, nanmean(cor_MSD, 1), 'color', analysis_info.cor_color*.75); hold on
    title('MSD for associated tracks')
    %
    ref_dist = cumsum(sqrt((ref_x1_aligned2ref_begin(:, 2:end) - ref_x1_aligned2ref_begin(:,1:(end-1))).^2 +...
    (ref_x2_aligned2ref_begin(:, 2:end) - ref_x2_aligned2ref_begin(:,1:(end-1))).^2), 2);
    cor_dist = cumsum(sqrt((cor_x1_aligned2cor_begin(:, 2:end) - cor_x1_aligned2cor_begin(:,1:(end-1))).^2 +...
    (cor_x2_aligned2cor_begin(:, 2:end) - cor_x2_aligned2cor_begin(:,1:(end-1))).^2), 2);
    figure('Position', [1240   168   333   260])
    plot([1:(tp_max-1)]*interval, nanmean(ref_dist, 1), 'color', analysis_info.ref_color); hold on
    plot([1:(tp_max-1)]*interval, nanmean(cor_dist, 1), 'color', analysis_info.cor_color); hold on
    title('distance traveled')}
    %}
    
    %% save intensity and lifetime data
    fid = fopen(lifetime_intensity_file, 'w');
    % Record the name of files
    for foldernum1 = 1:folder_ID1
        for foldernum2 = 1:folder_ID2
            if ~isempty(foldername_all(foldernum1, foldernum2).associated_track_file)
                fprintf(fid, '%s\n', ['Track file(', num2str(foldernum1),',',num2str(foldernum2),'):',...
                    fullfile(foldername_all(foldernum1, foldernum2).associated_track_folder, foldername_all(foldernum1, foldernum2).associated_track_file)]);
                fprintf(fid, '%s\n', ['Manual tag file(', num2str(foldernum1),',',num2str(foldernum2),'):',...
                    fullfile(foldername_all(foldernum1, foldernum2).machine_tag_folder, foldername_all(foldernum1, foldernum2).machine_tag_file)]);
            end
        end
    end
    % save reference lifetime
    fprintf(fid, '%s\n', 'ref lifetime (cor-refassociated)');
    fprintf(fid, '%s\t%s\n', 'lifetime', 'max intensity');
    fclose(fid);
    dlmwrite(lifetime_intensity_file, lifetime_max_intensity_assoc_ref, '-append', 'delimiter', '\t');
    fid = fopen(lifetime_intensity_file, 'a');
    fprintf(fid, '%s\n', 'corlifetime (cor-refassociated)');
    fprintf(fid, '%s\t%s\n', 'lifetime', 'max intensity');
    fclose(fid);
    dlmwrite(lifetime_intensity_file, lifetime_max_intensity_assoc_cor, '-append', 'delimiter', '\t');
    fid = fopen(lifetime_intensity_file, 'a');
    fprintf(fid, '%s\n', 'refonly lifetime');
    fprintf(fid, '%s\t%s\n', 'lifetime', 'max intensity');
    fclose(fid);
    dlmwrite(lifetime_intensity_file, lifetime_max_intensity_ref_only, '-append', 'delimiter', '\t');
    fid = fopen(lifetime_intensity_file, 'a');
    fprintf(fid, '%s\n', 'coronly lifetime');
    fprintf(fid, '%s\t%s\n', 'lifetime', 'max intensity');
    fclose(fid);
    dlmwrite(lifetime_intensity_file, lifetime_max_intensity_cor_only, '-append', 'delimiter', '\t');
    
    %% Save statistics for ref-only and cor-only tracks.
    fid_stat = fopen(track_stat_file_name, 'w');
    %% Record the name of files
    for foldernum1 = 1:folder_ID1
        for foldernum2 = 1:folder_ID2
            if ~isempty(foldername_all(foldernum1, foldernum2).associated_track_file)
                fprintf(fid_stat, '%s\n', ['Track file(', num2str(foldernum1),',',num2str(foldernum2),'):',...
                    fullfile(foldername_all(foldernum1, foldernum2).associated_track_folder, foldername_all(foldernum1, foldernum2).associated_track_file)]);
                fprintf(fid_stat, '%s\n', ['Manual tag file(', num2str(foldernum1),',',num2str(foldernum2),'):',...
                    fullfile(foldername_all(foldernum1, foldernum2).machine_tag_folder, foldername_all(foldernum1, foldernum2).machine_tag_file)]);
            end
        end
    end
    fprintf(fid_stat, '%s\n', ['Interval btw time points:', num2str(interval), ' sec.']);
    fprintf(fid_stat, '%s\n\n', ['Number of time points in the time series:', num2str(tp_max)]);
    
    fprintf(fid_stat, '%s\n\n', 'Assicated tracks ');
    fprintf(fid_stat, [repmat('%s\t', [1,19]),'%s\n\n'], ...
        'folder 1', 'folder 2', 'track ID', ...
        'ref track lifetime (sec)','cor track lifetime (sec)', 'profile length (sec)', ...
        'ref track beginninng time point', 'ref track disappearance time point', ...
        'cor track beginninng time point', 'cor track disappearance time point', ...
        'profile beginning time point','profile disappearance time point', ...
        'ref max intensity','ref mean intensity','mean change of intensity of ref track', 'mean change of position of ref track',...
        'cor max intensity','cor mean intensity','mean change of intensity of cor track', 'mean change of position of cor track');
    fclose(fid_stat);
    dlmwrite(track_stat_file_name, track_stat_asso, '-append', 'delimiter', '\t');
    fid_stat = fopen(track_stat_file_name, 'a');
    fprintf(fid_stat, '%s\n\n', 'reference track only');
    fprintf(fid_stat, [repmat('%s\t', [1,12]),'%s\n\n'], 'folder 1', 'folder 2', 'track ID', ...
        'track lifetime (sec)', 'profile length (sec)', 'track beginninng time point', ...
        'track disappearance time point', 'profile beginning time point','profile disappearance time point', ...
        'max intensity','mean intensity','mean change of intensity', 'mean change of position');
    fclose(fid_stat);
    dlmwrite(track_stat_file_name, track_stat_ref_only, '-append', 'delimiter', '\t');
    fid_stat = fopen(track_stat_file_name, 'a');
    fprintf(fid_stat, '\n%s\n\n', 'corresponding track only');
    fprintf(fid_stat, [repmat('%s\t', [1,12]),'%s\n\n'], 'folder 1', 'folder 2', 'track ID', ...
        'track lifetime (sec)', 'profile length (sec)', 'track beginninng time point', ...
        'track disappearance time point', 'profile beginning time point','profile disappearance time point', ...
        'max intensity','mean intensity','mean change of intensity', 'mean change of position');
    fclose(fid_stat);
    dlmwrite(track_stat_file_name, track_stat_cor_only, '-append', 'delimiter', '\t');
    % %%
    % [y, Fs] = audioread('notify.wav');
    % sound(y, Fs);
    
end
%}
function [h, N_sample] = drawAveragePlots(ref_profile, cor_profile, analysis_info, category, track_lifetimes, h)
        if or(numel(ref_profile)== 0, numel(cor_profile)== 0 )
            error('No tracks are found', ref_profile, cor_profile)
        end
        % h is a figure handler
        % cor_profile and ref_profile are coordinates of reference tracks and
        % corresponding tracks.
        tp_max = analysis_info.tp_max;
        interval = analysis_info.interval;
        
        % substitue NaNs with zero
        cor_profile_zero = cor_profile;
        cor_profile_zero(isnan(cor_profile_zero)) = 0;
        ref_profile_zero = ref_profile;
        ref_profile_zero(isnan(ref_profile_zero)) = 0;
        
        axis_min = min(find(~isnan(nanmean(cor_profile)), 1,'first'),find(~isnan(nanmean(ref_profile)), 1,'first'));
        axis_max = max(find(~isnan(nanmean(cor_profile)), 1,'last'),find(~isnan(nanmean(ref_profile)), 1,'last'));
        intensity_max = max(max(mean(cor_profile_zero)), max(mean(ref_profile_zero)));
        intensity_min = min(min(mean(cor_profile_zero)), min(mean(ref_profile_zero)));
        
        figure(h)
        subplot(5,2,[1,3]);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_profile);
        axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval,  min(min(cor_profile)), max(max(cor_profile))]);
        title('cor')
        subplot(5,2,[2,4]);plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_profile);
        axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval,  min(min(ref_profile)), max(max(ref_profile))]);
        title('ref')
        
        subplot(5,2,[7,9])
        plot(-(tp_max*interval-interval):interval:tp_max*interval, mean(cor_profile_zero), 'color',analysis_info.cor_color); hold on
        plot(-(tp_max*interval-interval):interval:tp_max*interval, mean(ref_profile_zero), 'color', analysis_info.ref_color);
        axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval, intensity_min, intensity_max]);
        % legend('cor', 'ref')
        xlabel('Time (sec)')
        ylabel('Average Intensity')
        
        % Divide samples to cohorts:
        % cohort colors
        if strcmp(analysis_info.ref_channel, 'RFP')
            % Each row is a different color.
            analysis_info.ref_cohort_color = [[240 51 163]/255;[185 45 153]/255;[146 39 143]/255];
            analysis_info.cor_cohort_color = [[0 187 246]/255;[0 167 157]/255;[0 148 68]/255];
        else
            analysis_info.ref_cohort_color = [[0 187 246]/255;[0 167 157]/255;[0 148 68]/255];
            analysis_info.cor_cohort_color = [[240 51 163]/255;[185 45 153]/255;[146 39 143]/255];
        end
        subplot(5,2,[8,10])
        % divide profiles so that the channel with longer lifetime is divided.
        if strcmp(category, 'associated')
            if mean(track_lifetimes(:,2))> mean(track_lifetimes(:,1))
                track_lifetimes=track_lifetimes(:,2);
                title_text='cor';
            else
                track_lifetimes=track_lifetimes(:,1);
                title_text='ref';
            end
        elseif strcmp(category, 'ref')
            title_text='ref';
        elseif strcmp(category, 'cor')
            title_text='cor';
        else
            display('Category should be one of the following three: associated, ref, cor');
            error
        end
        lifetime_range=linspace(min(track_lifetimes), max(track_lifetimes)+0.01, 4);
        intensity_max = -Inf;
        intensity_min = Inf;
        lifetime_text='';
        for i=1:3
            cohort_id = find(and(track_lifetimes>=lifetime_range(i),track_lifetimes<lifetime_range(i+1)));
            if numel(cohort_id)>=2
                plot(-(tp_max*interval-interval):interval:tp_max*interval, mean(cor_profile_zero(cohort_id,:)), 'color', analysis_info.cor_cohort_color(i,:)); hold on
                plot(-(tp_max*interval-interval):interval:tp_max*interval, mean(ref_profile_zero(cohort_id,:)), 'color', analysis_info.ref_cohort_color(i,:));
            elseif numel(cohort_id)==1
                plot(-(tp_max*interval-interval):interval:tp_max*interval, (cor_profile_zero(cohort_id,:)), 'color', analysis_info.cor_cohort_color(i,:)); hold on
                plot(-(tp_max*interval-interval):interval:tp_max*interval, (ref_profile_zero(cohort_id,:)), 'color', analysis_info.ref_cohort_color(i,:));
            end
            intensity_max=max(intensity_max,max(max(mean(cor_profile_zero(cohort_id,:))), max(mean(ref_profile_zero(cohort_id,:)))));
            intensity_min=min(intensity_min,min(min(mean(cor_profile_zero(cohort_id,:))), min(mean(ref_profile_zero(cohort_id,:)))));
            % Get info about the cohort
            lifetime_text=[lifetime_text, num2str(lifetime_range(i),'%.2f'),'-',num2str(lifetime_range(i+1),'%.2f'),' (n=', num2str(numel(cohort_id)),'), '];
        end
        
        for i=1:3
            
        end
        title([title_text,'-track lifetimes:', lifetime_text(1:(end-1))])
        axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval, intensity_min, intensity_max]);
        
        N_sample = size(cor_profile,1);
  