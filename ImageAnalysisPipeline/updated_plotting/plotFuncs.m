%% This program reads in a mat file of associated tracks determined with Associate_tracks.m
%% and manual_pick file from Manually_pick_associtated_tracks.m

%% csv files are saved with tracks horizontal, ref/corr for intensity and mean/std for averages

function [field_positions, lifetime_associated]= plotFuncs(S)
% list of all the plots that this file does
%% Mark 'true' if you want to save all the variables

funcs = {@plotAssociatedCorr, @plotAssociatedRef,...
    @plotUnassociatedCorr, @plotUnassociatedRef, @plotAlignCorrDisappear, @plotAlignCorrStart, ...
    @plotAlignRefEnd, @plotAlignRefStart, @plotAlignRefMax, @plotAlignCorrMax, @plotProfiles, ...
    @lifetimeHistogram, @AssociatedLifetimeHistogram, @plotHalfDropoff, @compare_stats};

toPlot = S.CNT;
saveVariables = S.save;
analyzeUnassociated = 0; % decide whether to plot and analyze unassociated tracks
    
%% Ask users for interval
interval = str2double(S.freq);
%input('What was the interval between time points?')

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

cor_aligned_half = [];
ref_aligned_half = [];


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
%load('previous_folder.mat')
data_folder = S.folder;
%uigetdir(previous_folder, 'Pick folder of all analysis.');
if isdir(data_folder)
    previous_folder = fullfile(data_folder, '');
    save('previous_folder', 'previous_folder')
else
    disp('You did not pick a folder.')
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

%% which channel to align to, for 50 percent dropoff alignment

if S.CNT(14)
    align_ref = input('Use reference or correlated channel to align from 50 percent dropoff? 1 for ref, 0 for corr: ');
end


%% Name of the results files
histogram_file_name = fullfile(data_folder, 'histogram_pool_data.txt');

%% Go through folders to pool data.
folder_ID1 = 0;
found = false;
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
                if exist(fullfile(track_folder_name, 'associated_tracks_ref_cor.mat'), 'file')
                    if (~S.plotAll)
                        prompt = ['Plot data at ', track_folder_name, ' ? 1 for yes, 0 for no'];
                        analyze = input(prompt);
                        if analyze == 0
                            continue;
                        end
                    end
                    found = true;
                    track_data = load(fullfile(track_folder_name, 'associated_tracks_ref_cor.mat'));
                    pixel_size = track_data.analysis_info.pixel_size; % in nm
                    if exist(fullfile(track_folder_name, clean_list_name), 'file')
                        machine_tag_folder = track_folder_name;
                        if strcmp(clean_list_name((end-2):end), 'mat')
                            load(fullfile(track_folder_name, clean_list_name));
                        elseif strcmp(clean_list_name((end-2):end), 'txt')
                            track_flag_list = dlmread(fullfile(track_folder_name, clean_list_name), '\t', 10, 0);
                        end
                    else
                        disp([fullfile(track_folder_name, clean_list_name), ' is not found.']);
                    end
                    disp(['using', fullfile(track_folder_name)]);
                else
                    disp([fullfile(track_folder_name, 'associated_tracks_ref_cor.mat'), ' is not found.']);
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
                
                %% Retrieve "timepoint violators" tracks
                % There are 2 subcategories: early (low) and late (high) tracks
                
                %MA: if the field doesn't exist (i.e. no tracks came from
                %this categroy) then make an empty vector.
                
                if isfield(track_data.track_stat_ref, 'ref_tracks_associated_low_breach')
                    ref_cor_tracks_4_low = [track_data.track_stat_ref.ref_tracks_associated_low_breach]';
                else
                    ref_cor_tracks_4_low = [];
                end
                
                if isfield(track_data.track_stat_ref, 'ref_tracks_unassociated_low_breach')
                    ref_only_tracks_4_low = [track_data.track_stat_ref.ref_tracks_unassociated_low_breach]';
                else
                    ref_only_tracks_4_low = [];
                end
                
                if isfield(track_data.track_stat_cor, 'cor_tracks_unassociated_low_breach')
                    cor_only_tracks_4_low = [track_data.track_stat_cor.cor_tracks_unassociated_low_breach]';
                else
                    cor_only_tracks_4_low = [];
                end
                
                ref_cor_tracks_4_low_num = numel(ref_cor_tracks_4_low);
                ref_only_tracks_4_low_num = numel(ref_only_tracks_4_low);
                cor_only_tracks_4_low_num = numel(cor_only_tracks_4_low);
                
                
                if isfield(track_data.track_stat_ref, 'ref_tracks_associated_high_breach')
                    ref_cor_tracks_4_high = [track_data.track_stat_ref.ref_tracks_associated_high_breach]';
                else
                    ref_cor_tracks_4_high = [];
                end
                
                if isfield(track_data.track_stat_ref, 'ref_tracks_unassociated_high_breach')
                    ref_only_tracks_4_high = [track_data.track_stat_ref.ref_tracks_unassociated_high_breach]';
                else
                    ref_only_tracks_4_high = [];
                end
                
                if isfield(track_data.track_stat_cor, 'cor_tracks_unassociated_high_breach')
                    cor_only_tracks_4_high = [track_data.track_stat_cor.cor_tracks_unassociated_high_breach]';
                else
                    cor_only_tracks_4_high = [];
                end
                
                ref_cor_tracks_4_high_num = numel(ref_cor_tracks_4_high);
                ref_only_tracks_4_high_num = numel(ref_only_tracks_4_high);
                cor_only_tracks_4_high_num = numel(cor_only_tracks_4_high);
                
                
                %% Retrive variables.
                tp_min = 1;
                tp_max = track_data.analysis_info.tp_max;
                
                
                %% Compile ref_cor and ref_ indecies, intensities
                for i = 1:ref_cor_num
                    cur_track = ref_cor_tracks(i);
                    
                    curCorTrack = (track_data.associated_tracks(cur_track, 1).cor_intensity_bg_corrected(1, :));
                    curRefTrack = (track_data.associated_tracks(cur_track, 1).ref_intensity_bg_corrected(1, :));
                    
                    allAssociatedTracks(jj).Cor = curCorTrack;
                    allAssociatedTracks(jj).Ref = curRefTrack;
                    jj = jj+1;
                end
                     
                curTrack = [];
                curCorTrack = [];
                curRefTrack = [];
                
                if analyzeUnassociated
                    
                    if ~isempty(track_data.unassociated_tracks)&&~isempty(ref_only_tracks)
                        
                        for i = 1:ref_only_num
                            cur_track = ref_only_tracks(i);
                            
                            curCorTrack = track_data.unassociated_tracks(cur_track,1).cor_intensity_bg_corrected(1,:);
                            curRefTrack = track_data.unassociated_tracks(cur_track,1).ref_intensity_bg_corrected(1,:);
                            
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
                
                if ~isempty(track_data.unassociated_tracks)&&~isempty(cor_only_tracks)
                    
                    for i = 1:cor_only_num
                        cur_track = cor_only_tracks(i);
                        
                        curCorTrack = track_data.unassociated_tracks(cur_track,1).cor_intensity_bg_corrected(1,:);
                        curRefTrack = track_data.unassociated_tracks(cur_track,1).ref_intensity_bg_corrected(1,:);
                        
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
                    cor_track = track_data.ROI_track_ref(ref_track, 1).associated_cor_track;
                    %%
                    lifetime_associated(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track, 1).ref_tp_last - track_data.ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval, ...
                        %%
                        (track_data.ROI_track_ref(ref_track, 1).cor_tp_last - track_data.ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval];
                end
                
                previous_index = size(lifetime_ref_only, 1);
                lifetime_ref_only = [lifetime_ref_only; NaN*ones(ref_only_num, 1)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    lifetime_ref_only(previous_index + i,1) = (track_data.ROI_track_ref(ref_track, 1).tp_last - track_data.ROI_track_ref(ref_track, 1).tp_first + 1)*interval;
                end
                
                previous_index = size(lifetime_cor_only, 1);
                lifetime_cor_only = [lifetime_cor_only; NaN*ones(cor_only_num, 1)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    lifetime_cor_only(previous_index + i,1) =(track_data.ROI_track_cor(cor_track, 1).tp_last - track_data.ROI_track_cor(cor_track, 1).tp_first + 1)*interval;
                end
               
                
                %% Determine lifetimes of timepoint violators for low timepoints
                
                previous_index = size(lifetime_associated_4_low, 1);
                lifetime_associated_4_low= [lifetime_associated_4_low; NaN*ones(ref_cor_tracks_4_low_num, 2)];
                for i = 1:ref_cor_tracks_4_low_num
                    ref_track_4_low = ref_cor_tracks_4_low(i);
                    cor_track_4_low = track_data.ROI_track_ref(ref_track_4_low, 1).associated_cor_track;
                    lifetime_associated_4_low(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track_4_low, 1).ref_tp_last - track_data.ROI_track_ref(ref_track_4_low, 1).ref_tp_first + 1)*interval, ...
                        (track_data.ROI_track_ref(ref_track_4_low, 1).cor_tp_last - track_data.ROI_track_ref(ref_track_4_low, 1).cor_tp_first + 1)*interval];
                end
                
                previous_index = size(lifetime_ref_only_4_low, 1);
                lifetime_ref_only_4_low = [lifetime_ref_only_4_low; NaN*ones(ref_only_tracks_4_low_num, 1)];
                for i = 1:ref_only_tracks_4_low_num
                    ref_track_4_low = ref_only_tracks_4_low(i);
                    lifetime_ref_only_4_low(previous_index + i,1) = (track_data.ROI_track_ref(ref_track_4_low, 1).tp_last - track_data.ROI_track_ref(ref_track_4_low, 1).tp_first + 1)*interval;
                end
                
                previous_index = size(lifetime_cor_only_4_low, 1);
                lifetime_cor_only_4_low = [lifetime_cor_only_4_low; NaN*ones(cor_only_tracks_4_low_num, 1)];
                for i = 1:cor_only_tracks_4_low_num
                    cor_track_4_low = cor_only_tracks_4_low(i);
                    lifetime_cor_only_4_low(previous_index + i,1) =(track_data.ROI_track_cor(cor_track_4_low, 1).tp_last - track_data.ROI_track_cor(cor_track_4_low, 1).tp_first + 1)*interval;
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
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated', num2str(ref_track, '%04i'),'.tif']),'tif');
                end
                
                % as PNG too
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_', num2str(ref_track, '%04i'),'.png']),'png');
                end
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_ref_only', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_ref_only', num2str(ref_track, '%04i'),'.tif']),'tif');
                end
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    imwrite([fit_contrast_ignore_zero(track_data.unassociated_tracks(cor_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.unassociated_tracks(cor_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_cor_only', num2str(cor_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(track_data.unassociated_tracks(cor_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.unassociated_tracks(cor_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_cor_only', num2str(cor_track, '%04i'),'.tif']),'tif');
                end
               
                
                % save montages of high tracks %CD
                
                for i = 1:ref_cor_tracks_4_high_num
                    
                    ref_track = ref_cor_tracks_4_high(i);
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_high', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_high', num2str(ref_track, '%04i'),'.tif']),'tif');
                    
                    % as PNG too
                    
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_high_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_high_', num2str(ref_track, '%04i'),'.png']),'png');
                end
                
                % save montages of low tracks %CD
                
                for i = 1:ref_cor_tracks_4_low_num
                    
                    ref_track = ref_cor_tracks_4_low(i);
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_low', num2str(ref_track, '%04i'),'.tif']),'tif');
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_low', num2str(ref_track, '%04i'),'.tif']),'tif');
                    
                    % as PNG too
                    
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_montage); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_montage)],...
                        fullfile(montage_folder_name, ['montage_associated_low_', num2str(ref_track, '%04i'),'.png']),'png');
                    imwrite([fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).ref_movie); fit_contrast_ignore_zero(track_data.associated_tracks(ref_track, 1).cor_movie)],...
                        fullfile(movie_folder_name, ['movie_associated_low_', num2str(ref_track, '%04i'),'.png']),'png');
                end                
                
                %1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to maximum of cor
                previous_index = size(cor_aligned2cormax, 1);
                cor_aligned2cormax = [cor_aligned2cormax; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2cormax = [ref_aligned2cormax; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    %smoothened_curve = medfilt1(track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :), 5);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    %         align the curve so that maximum of the curve is coincident
                    %         [max_value, max_index] = max(smoothened_curve);
                    [max_value, max_index] = max(track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_index = median(max_index);
                    cor_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2cormax(previous_index +i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                %% Align reference only tracks to maximum of cor
                previous_index = size(ref_only_ref_aligned2cormax, 1);
                ref_only_ref_aligned2cormax = [ref_only_ref_aligned2cormax; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2cormax = [ref_only_cor_aligned2cormax; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    smoothened_ref= medfilt1(track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :), 5);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_ref);
                    ref_index = median(max_index);
                    ref_only_cor_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %% Align corresponding only tracks to maximum of cor
                previous_index = size(cor_only_ref_aligned2cormax, 1);
                cor_only_ref_aligned2cormax = [cor_only_ref_aligned2cormax; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2cormax = [cor_only_cor_aligned2cormax; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    smoothened_curve = medfilt1(track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :), 5);
                    profile_length = track_data.ROI_track_cor(cor_track, 1).tp_post - track_data.ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_curve);
                    ref_index = median(max_index);
                    cor_only_cor_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2cormax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        track_data.unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to maximum of ref
                previous_index = size(cor_aligned2refmax, 1);
                cor_aligned2refmax = [cor_aligned2refmax; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2refmax = [ref_aligned2refmax; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    %         smoothened_curve = medfilt1(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    %         align the curve so that maximum of the curve is coincident
                    %         [max_value, max_index] = max(smoothened_curve);
                    [max_value, max_index] = max(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                    ref_index = median(max_index);
                    cor_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2refmax(previous_index +i, (tp_max - ref_index):(tp_max + profile_length - ref_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                
                %% Align reference only tracks to maximum of ref
                previous_index = size(ref_only_ref_aligned2refmax, 1);
                ref_only_ref_aligned2refmax = [ref_only_ref_aligned2refmax; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2refmax = [ref_only_cor_aligned2refmax; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    smoothened_ref= medfilt1(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_ref);
                    ref_index = median(max_index);
                    ref_only_cor_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %% Align corresponding only tracks to maximum of ref
                previous_index = size(cor_only_ref_aligned2refmax, 1);
                cor_only_ref_aligned2refmax = [cor_only_ref_aligned2refmax; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2refmax = [cor_only_cor_aligned2refmax; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    smoothened_curve = medfilt1(track_data.unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = track_data.ROI_track_cor(cor_track, 1).tp_post - track_data.ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that maximum of the curve is coincident
                    [max_value, max_index] = max(smoothened_curve);
                    ref_index = median(max_index);
                    cor_only_cor_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2refmax(previous_index + i, (tp_max - ref_index):(tp_max + profile_length - ref_index - 1)) = ...
                        track_data.unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the beginning of ref
                previous_index = size(cor_aligned2ref_begin, 1);
                cor_aligned2ref_begin = [cor_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2ref_begin = [ref_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align curve to the beginning of the reference track.
                    ref_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_first;
                    first_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    cor_aligned2ref_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2ref_begin(previous_index +i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Coordinates: Align reference tracks to the beginning of reference track
                previous_index = size(ref_x1_aligned2ref_begin, 1);
                ref_x1_aligned2ref_begin = [ref_x1_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max)];
                ref_x2_aligned2ref_begin = [ref_x2_aligned2ref_begin; NaN*ones(ref_cor_num, tp_max)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).ref_tp_last - track_data.ROI_track_ref(ref_track, 1).ref_tp_first + 1;
                    first_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    last_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_last - track_data.ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    ref_x1_aligned2ref_begin(previous_index + i,1:profile_length) = track_data.ROI_track_ref(ref_track,1).ref_coord_extrapol(first_index:last_index,1)';
                    ref_x2_aligned2ref_begin(previous_index + i,1:profile_length) = track_data.ROI_track_ref(ref_track,1).ref_coord_extrapol(first_index:last_index,2)';
                end
                
                
                %% Align reference only tracks to beginning of ref
                previous_index = size(ref_only_ref_aligned2ref_begin, 1);
                ref_only_ref_aligned2ref_begin = [ref_only_ref_aligned2ref_begin; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2ref_begin = [ref_only_cor_aligned2ref_begin; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    smoothened_ref= medfilt1(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :), 5);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align curve to the beginning of the reference track.
                    ref_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_first;
                    first_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    ref_only_cor_aligned2ref_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2ref_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the end of reference track
                previous_index = size(cor_aligned2ref_end, 1);
                cor_aligned2ref_end = [cor_aligned2ref_end; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2ref_end = [ref_aligned2ref_end; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_last;
                    last_index_from_end = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).ref_tp_last;
                    cor_aligned2ref_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2ref_end(previous_index +i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                
                %% Align reference only tracks to the end of reference track
                previous_index = size(ref_only_ref_aligned2ref_end, 1);
                ref_only_ref_aligned2ref_end = [ref_only_ref_aligned2ref_end; NaN*ones(ref_only_num, tp_max*2)];
                ref_only_cor_aligned2ref_end = [ref_only_cor_aligned2ref_end; NaN*ones(ref_only_num, tp_max*2)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index = track_data.ROI_track_ref(ref_track, 1).ref_tp_last;
                    last_index_from_end = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).ref_tp_last;
                    ref_only_cor_aligned2ref_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                    ref_only_ref_aligned2ref_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the beginning of corresponding track
                previous_index = size(cor_aligned2cor_begin, 1);
                cor_aligned2cor_begin = [cor_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2cor_begin = [ref_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the beginning of the curve is coincident
                    ref_index = track_data.ROI_track_ref(ref_track, 1).cor_tp_first;
                    first_index = track_data.ROI_track_ref(ref_track, 1).cor_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    cor_aligned2cor_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2cor_begin(previous_index +i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        (track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Coordinates: Align corresponding tracks to the beginning of corresponding track
                previous_index = size(cor_x1_aligned2cor_begin, 1);
                cor_x1_aligned2cor_begin = [cor_x1_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max)];
                cor_x2_aligned2cor_begin = [cor_x2_aligned2cor_begin; NaN*ones(ref_cor_num, tp_max)];
                
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).cor_tp_last - track_data.ROI_track_ref(ref_track, 1).cor_tp_first + 1;
                    first_index = track_data.ROI_track_ref(ref_track, 1).cor_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    last_index = track_data.ROI_track_ref(ref_track, 1).cor_tp_last - track_data.ROI_track_ref(ref_track, 1).tp_prev_all +1;
                    cor_x1_aligned2cor_begin(previous_index + i,1:profile_length) = track_data.ROI_track_ref(ref_track,1).cor_coord_extrapol(first_index:last_index,1)';
                    cor_x2_aligned2cor_begin(previous_index + i,1:profile_length) = track_data.ROI_track_ref(ref_track,1).cor_coord_extrapol(first_index:last_index,2)';
                end
                
                
                %% Align corresponding only tracks to the beginning of corresponding track
                previous_index = size(cor_only_ref_aligned2cor_begin, 1);
                cor_only_ref_aligned2cor_begin = [cor_only_ref_aligned2cor_begin; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2cor_begin = [cor_only_cor_aligned2cor_begin; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    profile_length = track_data.ROI_track_cor(cor_track, 1).tp_post - track_data.ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that the beginning of the curve is coincident
                    ref_index = track_data.ROI_track_cor(cor_track, 1).tp_first;
                    first_index = track_data.ROI_track_cor(cor_track, 1).tp_first - track_data.ROI_track_cor(cor_track, 1).tp_prev + 1;
                    cor_only_cor_aligned2cor_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2cor_begin(previous_index + i, (tp_max - first_index):(tp_max + profile_length - first_index-1)) = ...
                        track_data.unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                    %                     ref_index = track_data.ROI_track_ref(ref_track, 1).cor_tp_first;
                    %                     first_index = track_data.ROI_track_ref(ref_track, 1).cor_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                end
                
                
                %6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Align tracks (associated corresponding and ref) to the disappearance of corresponding track
                previous_index = size(cor_aligned2cor_end, 1);
                cor_aligned2cor_end = [cor_aligned2cor_end; NaN*ones(ref_cor_num, tp_max*2)];
                ref_aligned2cor_end = [ref_aligned2cor_end; NaN*ones(ref_cor_num, tp_max*2)];
                
                for i = 1:ref_cor_num
                    %%
                    ref_track = ref_cor_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index = track_data.ROI_track_ref(ref_track, 1).cor_tp_last;
                    last_index_from_end = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).cor_tp_last;
                    cor_aligned2cor_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :));
                    ref_aligned2cor_end(previous_index +i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        (track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :));
                end
                
                
                %% Align corresponding only tracks to the disappearance of corresponding track
                previous_index = size(cor_only_ref_aligned2cor_end, 1);
                cor_only_ref_aligned2cor_end = [cor_only_ref_aligned2cor_end; NaN*ones(cor_only_num, tp_max*2)];
                cor_only_cor_aligned2cor_end = [cor_only_cor_aligned2cor_end; NaN*ones(cor_only_num, tp_max*2)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    profile_length = track_data.ROI_track_cor(cor_track, 1).tp_post - track_data.ROI_track_cor(cor_track, 1).tp_prev + 1;
                    % align the curve so that the ending of the curve is coincident
                    ref_index =  track_data.ROI_track_cor(cor_track, 1).tp_last;
                    last_index_from_end =  track_data.ROI_track_cor(cor_track, 1).tp_post -  track_data.ROI_track_cor(cor_track, 1).tp_last;
                    cor_only_cor_aligned2cor_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                    cor_only_ref_aligned2cor_end(previous_index + i, (tp_max - (profile_length - last_index_from_end) + 1):(tp_max + last_index_from_end)) = ...
                        track_data.unassociated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                end
                
                %% Align corresponding and reference tracks to 50% dropoff of the reference or correlated track 
                previous_index = size(ref_aligned_half, 1);
                ref_aligned_half = [ref_aligned_half; NaN*ones(ref_cor_num, tp_max*2)];
                cor_aligned_half = [cor_aligned_half; NaN*ones(ref_cor_num, tp_max*2)];
                
                if toPlot(14) == 1
                    if align_ref
                        for i = 1:ref_cor_num
                            ref_track = ref_cor_tracks(i);
                            smoothened_ref= medfilt1(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :), 5);
                            profile_length = track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1;
                            % align the curve so that maximum of the curve is coincident
                            [max_value, ~] = max(smoothened_ref);
                            max_value = median(max_value);
                            n = numel(smoothened_ref);
                            while(smoothened_ref(n) < (max_value / 2))
                                n = n - 1;
                            end
                            ref_aligned_half(previous_index + i, (tp_max - n):(tp_max + profile_length - n - 1)) = ...
                                track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :);
                            cor_aligned_half(previous_index + i, (tp_max - n):(tp_max + profile_length - n - 1)) = ...
                                track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :);
                        end
                    else
                        for i = 1:ref_cor_num
                            cor_track = ref_cor_tracks(i);
                            smoothened_cor= medfilt1(track_data.associated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :), 5);
                            profile_length = track_data.ROI_track_ref(cor_track, 1).tp_post_all - track_data.ROI_track_ref(cor_track, 1).tp_prev_all + 1;
                            % align the curve so that maximum of the curve is coincident
                            [max_value, ~] = max(smoothened_cor);
                            n = numel(smoothened_cor);
                            while(smoothened_cor(n) < (max_value / 2))
                                n = n - 1;
                            end
                            ref_aligned_half(previous_index + i, (tp_max - n):(tp_max + profile_length - n - 1)) = ...
                                track_data.associated_tracks(cor_track, 1).ref_intensity_bg_corrected(1, :);
                            cor_aligned_half(previous_index + i, (tp_max - n):(tp_max + profile_length - n - 1)) = ...
                                track_data.associated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :);
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Coordinates: Align reference only tracks to the beginning of reference track
                previous_index = size(ref_only_x1_aligned2begin, 1);
                ref_only_x1_aligned2begin = [ref_only_x1_aligned2begin; NaN*ones(ref_only_num, tp_max)];
                ref_only_x2_aligned2begin = [ref_only_x2_aligned2begin; NaN*ones(ref_only_num, tp_max)];
                
                for i = 1:ref_only_num
                    ref_track = ref_only_tracks(i);
                    profile_length = track_data.ROI_track_ref(ref_track, 1).tp_last - track_data.ROI_track_ref(ref_track, 1).tp_first + 1;
                    ref_only_x1_aligned2begin(previous_index + i,1:profile_length) = track_data.ROI_track_ref(ref_track,1).center_coord(:,1)';
                    ref_only_x2_aligned2begin(previous_index + i,1:profile_length) = track_data.ROI_track_ref(ref_track,1).center_coord(:,2)';
                end
                
                %% Coordinates: Align corresponding only tracks to the beginning of corresponding track
                previous_index = size(cor_only_x1_aligned2begin, 1);
                cor_only_x1_aligned2begin = [cor_only_x1_aligned2begin; NaN*ones(cor_only_num, tp_max)];
                cor_only_x2_aligned2begin = [cor_only_x2_aligned2begin; NaN*ones(cor_only_num, tp_max)];
                
                for i = 1:cor_only_num
                    cor_track = cor_only_tracks(i);
                    profile_length = track_data.ROI_track_cor(cor_track, 1).tp_last - track_data.ROI_track_cor(cor_track, 1).tp_first + 1;
                    cor_only_x1_aligned2begin(previous_index + i,1:profile_length) = track_data.ROI_track_cor(cor_track,1).center_coord(:,1)';
                    cor_only_x2_aligned2begin(previous_index + i,1:profile_length) = track_data.ROI_track_cor(cor_track,1).center_coord(:,2)';
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Get statistics of reference tracks associated with corresponding tracks and reference only tracks: Lifetime and maximum intensity
                previous_index = size(lifetime_max_intensity_assoc_ref, 1);
                lifetime_max_intensity_assoc_ref= [lifetime_max_intensity_assoc_ref; NaN*ones(ref_cor_num, 2)];
                lifetime_mean_intensity_assoc_ref= [lifetime_mean_intensity_assoc_ref; NaN*ones(ref_cor_num, 2)];
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    lifetime_max_intensity_assoc_ref(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track, 1).ref_tp_last - track_data.ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval, ...
                        max(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_assoc_ref(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track, 1).ref_tp_last - track_data.ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval, ...
                        mean(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
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
                    track_stat_asso(previous_index+i,4) = (track_data.ROI_track_ref(ref_track, 1).ref_tp_last - track_data.ROI_track_ref(ref_track, 1).ref_tp_first + 1)*interval; % ref track lifetime
                    track_stat_asso(previous_index+i,5) = (track_data.ROI_track_ref(ref_track, 1).cor_tp_last - track_data.ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval; % cor track lifetime
                    track_stat_asso(previous_index+i,6) = (track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all)*interval; % profile length
                    track_stat_asso(previous_index+i,7:8) = [track_data.ROI_track_ref(ref_track, 1).ref_tp_first, track_data.ROI_track_ref(ref_track, 1).ref_tp_last]; % ref track beginninng/ending time
                    track_stat_asso(previous_index+i,9:10) = [track_data.ROI_track_ref(ref_track, 1).cor_tp_first, track_data.ROI_track_ref(ref_track, 1).cor_tp_last]; % cor track beginninng/ending time
                    track_stat_asso(previous_index+i,11:12) = [track_data.ROI_track_ref(ref_track, 1).tp_prev_all, track_data.ROI_track_ref(ref_track, 1).tp_post_all]; % profile beginninng/ending time
                    track_stat_asso(previous_index+i,13) = max(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % ref max intensity
                    track_stat_asso(previous_index+i,14) = mean(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % ref mean intensity
                    track_stat_asso(previous_index+i,17) = max(track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :)); % cor max intensity
                    track_stat_asso(previous_index+i,18) = mean(track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :)); % cor mean intensity
                    % determine change in intensity of the reference track
                    track_index = (track_data.ROI_track_ref(ref_track, 1).ref_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1):...
                        (track_data.ROI_track_ref(ref_track, 1).ref_tp_last - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1);
                    track_stat_asso(previous_index+i,15) = mean((track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(2:end)) -...
                        track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_asso(previous_index+i,16) = mean(sum((track_data.ROI_track_ref(ref_track, 1).center_coord(2:end,:) -...
                        track_data.ROI_track_ref(ref_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    % determine change in intensity of the corresponding track
                    track_index = (track_data.ROI_track_ref(ref_track, 1).cor_tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1):...
                        (track_data.ROI_track_ref(ref_track, 1).cor_tp_last - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1);
                    track_stat_asso(previous_index+i,19) = mean((track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, track_index(2:end)) -...
                        track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_asso(previous_index+i,20) = mean(sum((track_data.ROI_track_ref(ref_track, 1).center_coord(2:end,:) -...
                        track_data.ROI_track_ref(ref_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    
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
                    lifetime_max_intensity_ref_only(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track, 1).tp_last - track_data.ROI_track_ref(ref_track, 1).tp_first + 1)*interval, ...
                        max(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_ref_only(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track, 1).tp_last - track_data.ROI_track_ref(ref_track, 1).tp_first + 1)*interval, ...
                        mean(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :))];
                    
                    track_stat_ref_only(previous_index+i,1) = folder_ID1;
                    track_stat_ref_only(previous_index+i,2) = folder_ID2;
                    track_stat_ref_only(previous_index+i,3) = ref_track; %track ID
                    track_stat_ref_only(previous_index+i,4) = (track_data.ROI_track_ref(ref_track, 1).tp_last - track_data.ROI_track_ref(ref_track, 1).tp_first + 1)*interval; % track lifetime
                    track_stat_ref_only(previous_index+i,5) = (track_data.ROI_track_ref(ref_track, 1).tp_post_all - track_data.ROI_track_ref(ref_track, 1).tp_prev_all)*interval; % profile length
                    track_stat_ref_only(previous_index+i,6:7) = [track_data.ROI_track_ref(ref_track, 1).tp_first, track_data.ROI_track_ref(ref_track, 1).tp_last]; % track beginninng/ending time
                    track_stat_ref_only(previous_index+i,8:9) = [track_data.ROI_track_ref(ref_track, 1).tp_prev_all, track_data.ROI_track_ref(ref_track, 1).tp_post_all]; % profile beginninng/ending time
                    track_stat_ref_only(previous_index+i,10) = max(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % max intensity
                    track_stat_ref_only(previous_index+i,11) = mean(track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, :)); % mean intensity
                    % determine change in intensity
                    track_index = (track_data.ROI_track_ref(ref_track, 1).tp_first - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1):...
                        (track_data.ROI_track_ref(ref_track, 1).tp_last - track_data.ROI_track_ref(ref_track, 1).tp_prev_all + 1);
                    track_stat_ref_only(previous_index+i,12) = mean((track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(2:end)) -...
                        track_data.associated_tracks(ref_track, 1).ref_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_ref_only(previous_index+i,13) = mean(sum((track_data.ROI_track_ref(ref_track, 1).center_coord(2:end,:) -...
                        track_data.ROI_track_ref(ref_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    
                end
                
                
                %% Get statistics of corresponding tracks associated with reference tracks: Lifetime and maximum intensity
                previous_index = size(lifetime_max_intensity_assoc_cor, 1);
                lifetime_max_intensity_assoc_cor= [lifetime_mean_intensity_assoc_cor; NaN*ones(ref_cor_num, 2)];
                lifetime_mean_intensity_assoc_cor= [lifetime_mean_intensity_assoc_cor; NaN*ones(ref_cor_num, 2)];
                for i = 1:ref_cor_num
                    ref_track = ref_cor_tracks(i);
                    lifetime_max_intensity_assoc_cor(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track, 1).cor_tp_last - track_data.ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval, ...
                        max(track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_assoc_cor(previous_index + i,:) = [(track_data.ROI_track_ref(ref_track, 1).cor_tp_last - track_data.ROI_track_ref(ref_track, 1).cor_tp_first + 1)*interval, ...
                        mean(track_data.associated_tracks(ref_track, 1).cor_intensity_bg_corrected(1, :))];
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
                    lifetime_max_intensity_cor_only(previous_index + i,:) = [(track_data.ROI_track_cor(cor_track, 1).tp_last - track_data.ROI_track_cor(cor_track, 1).tp_first + 1)*interval, ...
                        max(track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :))];
                    lifetime_mean_intensity_cor_only(previous_index + i,:) = [(track_data.ROI_track_cor(cor_track, 1).tp_last - track_data.ROI_track_cor(cor_track, 1).tp_first + 1)*interval, ...
                        mean(track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :))];
                    
                    track_stat_cor_only(previous_index+i,1) = folder_ID1;
                    track_stat_cor_only(previous_index+i,2) = folder_ID2;
                    track_stat_cor_only(previous_index+i,3) = cor_track; %track ID
                    track_stat_cor_only(previous_index+i,4) = (track_data.ROI_track_cor(cor_track, 1).tp_last - track_data.ROI_track_cor(cor_track, 1).tp_first + 1)*interval; % track lifetime
                    track_stat_cor_only(previous_index+i,5) = (track_data.ROI_track_cor(cor_track, 1).tp_post - track_data.ROI_track_cor(cor_track, 1).tp_prev)*interval; % profile length
                    track_stat_cor_only(previous_index+i,6:7) = [track_data.ROI_track_cor(cor_track, 1).tp_first, track_data.ROI_track_cor(cor_track, 1).tp_last]; % track beginninng/ending time
                    track_stat_cor_only(previous_index+i,8:9) = [track_data.ROI_track_cor(cor_track, 1).tp_prev, track_data.ROI_track_cor(cor_track, 1).tp_post]; % profile beginninng/ending time
                    track_stat_cor_only(previous_index+i,10) = max(track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :)); % max intensity
                    track_stat_cor_only(previous_index+i,11) = mean(track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, :)); % mean intensity
                    % determine change in intensity
                    track_index = (track_data.ROI_track_cor(cor_track, 1).tp_first - track_data.ROI_track_cor(cor_track, 1).tp_prev + 1):...
                        (track_data.ROI_track_cor(cor_track, 1).tp_last - track_data.ROI_track_cor(cor_track, 1).tp_prev + 1);
                    track_stat_cor_only(previous_index+i,12) = mean((track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, track_index(2:end)) -...
                        track_data.unassociated_tracks(cor_track, 1).cor_intensity_bg_corrected(1, track_index(1:(end-1)))).^2); % mean change of intensity
                    track_stat_cor_only(previous_index+i,13) = mean(sum((track_data.ROI_track_cor(cor_track, 1).center_coord(2:end,:) -...
                        track_data.ROI_track_cor(cor_track, 1).center_coord(1:(end-1),:)).^2, 2)); % mean change of position
                    
                end
                % inventory of the indecies of "allAssociatedTracks" that
                % separate each dataset
                betweenDatasets(jjj) = jj-1;
                jjj = jjj+1;
            end
        end
        if found == false
            disp('could not find any data');
        end
    end
end
%% Determine ref/corcolor
if strcmp(track_data.analysis_info.ref_channel, 'RFP')
    track_data.analysis_info.ref_color = [1, 0, 0.8];
    track_data.analysis_info.cor_color = [0, 0.8, 0.1];
else
    track_data.analysis_info.ref_color = [0, 0.8, 0.1];
    track_data.analysis_info.cor_color = [1, 0, 0.8];
end
track_data.analysis_info.interval = interval;
%% Get screen size for figures.
scrsz = get(groot,'ScreenSize');
fig_width = 400;
fig_height = 250;

bins = linspace((tp_min-1)*interval, (tp_max-1)*interval, 13);
bin_size = bins(2) - bins(1);
bin_center = (bins(1:end-1) + bins(2:end))/2;
bin_center = [bin_center, bin_center(end) + bin_size];

lifetime_hist_ref_only = histc(lifetime_ref_only, bins);
lifetime_hist_cor_only = histc(lifetime_cor_only, bins);
lifetime_hist_ref= histc(lifetime_associated(:,1), bins);
lifetime_hist_cor= histc(lifetime_associated(:,2), bins);

save(fullfile(data_folder, 'lifetimeData.mat'), 'lifetime_associated')

%% Record the name of files
fid = fopen(histogram_file_name, 'w');
for foldernum1 = 1:folder_ID1
    for foldernum2  = 1:folder_ID2
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

% Some issue here, needed to transpose one array for some reason
end

size_x_cor_only = size(lifetime_hist_cor_only,1);
size_y_cor_only = size(lifetime_hist_cor_only,2);

if size_y_cor_only>size_x_cor_only
      lifetime_hist_cor_only=lifetime_hist_cor_only';
      disp('transposed lifetime hist cor only to save to text file')
end

% transpose if need to (not sure why it's sometimes not transposed)

size_x_ref_only = size(lifetime_hist_ref_only,1);
size_y_ref_only = size(lifetime_hist_ref_only,2);

if size_y_ref_only>size_x_ref_only
      lifetime_hist_ref_only=lifetime_hist_ref_only';
      disp('transposed lifetime hist ref only to save to text file')

dlmwrite(histogram_file_name, [bins', (bins+bin_size)', bin_center', lifetime_hist_ref, lifetime_hist_cor, lifetime_hist_ref_only,lifetime_hist_cor_only], '-append', 'delimiter', '\t');
% dlmwrite(histogram_file_name, [bins', (bins+bin_size)', bin_center', lifetime_hist_ref, lifetime_hist_cor, lifetime_hist_ref_only,lifetime_hist_cor_only'], '-append', 'delimiter', '\t');

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
%% Select which functions to plot
disp(toPlot);
for k = 1:length(toPlot)
    if (toPlot(k) == 1)
        funcs{k}();
    end
end

%% Write data out to CSV files
if saveVariables
    loc = strsplit(track_data.analysis_info.ref_csv_name, '/');
    assoc = strcat('/',fullfile(loc{1:numel(loc) - 3}, 'associated_intensity.csv'));
    ref = strcat('/',fullfile(loc{1:numel(loc) - 3}, 'ref_intensity.csv'));
    corr = strcat('/',fullfile(loc{1:numel(loc) - 3}, 'corr_intensity.csv'));
    assoc_mean = strcat('/',fullfile(loc{1:numel(loc) - 3}, 'associated_mean.csv'));
    ref_mean = strcat('/',fullfile(loc{1:numel(loc) - 3}, 'ref_mean.csv'));
    corr_mean = strcat('/',fullfile(loc{1:numel(loc) - 3}, 'corr_mean.csv'));
    
    index = 1;
    for arr = 1:numel(track_data.associated_tracks)
        if ~isempty(track_data.associated_tracks(arr).ref_intensity_bg_corrected)...
                && (track_data.associated_tracks(arr).ref_intensity_bg_corrected(1) ~= 1)
            dlmwrite(assoc,[track_data.associated_tracks(arr).ref_intensity_bg_corrected;...
                track_data.associated_tracks(arr).cor_intensity_bg_corrected], '-append');
            dlmwrite(assoc_mean, [mean(track_data.associated_tracks(arr).ref_intensity_bg_corrected);...
                std(track_data.associated_tracks(arr).ref_intensity_bg_corrected); ...
                mean(track_data.associated_tracks(arr).cor_intensity_bg_corrected);
                std(track_data.associated_tracks(arr).cor_intensity_bg_corrected)], '-append');
        end
    end
    
    for arr = 1:numel(track_data.unassociated_tracks)
        if ~isempty(track_data.unassociated_tracks(arr).ref_intensity_bg_corrected) ...
                && ~isnan(track_data.unassociated_tracks(arr).ref_intensity_bg_corrected(1) ~= 1)
            dlmwrite(ref,track_data.unassociated_tracks(arr).ref_intensity_bg_corrected, '-append');
            dlmwrite(ref_mean, [mean(track_data.unassociated_tracks(arr).ref_intensity_bg_corrected);...
                std(track_data.unassociated_tracks(arr).ref_intensity_bg_corrected)], '-append');
        end
    end
    for arr = 1:numel(track_data.unassociated_tracks)
        if ~isempty(track_data.unassociated_tracks(arr).cor_intensity_bg_corrected) ...
                && (track_data.unassociated_tracks(arr).cor_intensity_bg_corrected(1) ~= 1)
            dlmwrite(corr,track_data.unassociated_tracks(arr).cor_intensity_bg_corrected,'-append');
            dlmwrite(corr_mean, [mean(track_data.unassociated_tracks(arr).cor_intensity_bg_corrected);...
                std(track_data.unassociated_tracks(arr).cor_intensity_bg_corrected)], '-append');
        end
    end
end
%% Figure: Histogram of lifetime.
    function [] = lifetimeHistogram
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
    end

%% Histogram of lifetime, show only the associated tracks only.
    function [] = AssociatedLifetimeHistogram
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
    end

    %% plot all profiles
    function [] = plotProfiles
        figure('OuterPosition', [50+fig_width*4, scrsz(4)-fig_height*2, fig_width*2,   fig_height*2]);
        for i=1:min(20, size(cor_aligned2cormax, 1))
            subplot(4,5,i);
            plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_aligned2cormax(i, :),'Color', track_data.analysis_info.cor_color); 
            hold on;
            subplot(4,5,i);
            plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_aligned2cormax(i, :),'Color',  track_data.analysis_info.ref_color);
            first_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'first');
            last_tp = find(~isnan(cor_aligned2cormax(i, :)), 1, 'last');
            axis([(first_tp-tp_max)*2, (last_tp-tp_max)*2, min(min(ref_aligned2cormax)), max(max(ref_aligned2cormax))]);
        end
    end

    %% Plot tracks of reference and corresponding tracks aligned to the max of corresponding track.
    function [] = plotAlignCorrMax
        if num_associated > 0
            % substitue NaNs with zero
            h = figure('OuterPosition', [1, scrsz(4)-fig_height*2, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_aligned2cormax, cor_aligned2cormax, track_data.analysis_info, 'associated', lifetime_associated, h);
            figure(h)
            suplabel(['Aligned to the max of cor. track, n=', num2str(N_sample)],'t');
        end
    end

    %% Plot tracks of reference and corresponding tracks aligned to the max of reference track.
    function [] = plotAlignRefMax
        if num_associated > 0
            h = figure('OuterPosition', [fig_width, scrsz(4)-fig_height*2, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_aligned2refmax, cor_aligned2refmax, track_data.analysis_info, 'associated', lifetime_associated, h);
            figure(h)
            suplabel(['Aligned to the max of ref. track, n=', num2str(N_sample)],'t');
        end
    end
    %% Plot tracks of reference and corresponding tracks aligned to half of reference or correlated
    function [] = plotHalfDropoff
        if num_associated > 0
            h = figure('OuterPosition', [fig_width, scrsz(4)-fig_height*2, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_aligned_half, cor_aligned_half, track_data.analysis_info, 'associated', lifetime_associated, h);
            figure(h)
            
            if align_ref
                suplabel(['Aligned to 50% dropoff of ref, n=', num2str(N_sample)],'t');
            else
                suplabel(['Aligned to 50% dropoff of cor, n=', num2str(N_sample)],'t');
            end
        end
    end

    %% Plot tracks of reference and corresponding tracks aligned to the beinning of reference track.
    function [] = plotAlignRefStart
        if num_associated > 0
            h = figure('OuterPosition', [fig_width*2, scrsz(4)-fig_height*2, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_aligned2ref_begin, cor_aligned2ref_begin, track_data.analysis_info, 'associated', lifetime_associated, h);
            figure(h)
            suplabel(['Aligned to the beginning of ref. track, n=', num2str(N_sample)],'t');
        end
    end
    %% Plot tracks of reference and corresponding tracks aligned to the ending of reference track.
    function [] = plotAlignRefEnd
        if num_associated > 0
            h = figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*2, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_aligned2ref_end, cor_aligned2ref_end, track_data.analysis_info, 'associated', lifetime_associated, h);
            figure(h)
            suplabel(['Aligned to the ending of ref. track, n=', num2str(N_sample)],'t');
        end
    end
    
    %% Plot tracks of reference and corresponding tracks aligned to the beinning of corresponding track.
    function [] = plotAlignCorrStart
        if num_associated > 0
            h = figure('OuterPosition', [1, scrsz(4)-fig_height*3, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_aligned2cor_begin, cor_aligned2cor_begin, track_data.analysis_info, 'associated', lifetime_associated, h);
            figure(h)
            suplabel(['Aligned to the beinning of cor. track, n=', num2str(N_sample)],'t');
        end
    end

    %% 6. Plot tracks of reference and corresponding tracks aligned to the disappearance of corresponding track.
    function [] = plotAlignCorrDisappear
        if num_associated > 0
            h = figure('OuterPosition', [fig_width, scrsz(4)-fig_height*3, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_aligned2cor_end, cor_aligned2cor_end, track_data.analysis_info, 'associated', lifetime_associated, h);
            figure(h)
            suplabel(['Aligned to the disappearance of cor. track, n=', num2str(N_sample)],'t');
        end
    end
%% Plot unassociated reference tracks
    function [] = plotUnassociatedRef
        if num_ref_only>0
            h = figure('OuterPosition', [fig_width*2, scrsz(4)-fig_height*3, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_only_ref_aligned2ref_begin, ref_only_cor_aligned2ref_begin, track_data.analysis_info, 'ref', lifetime_ref_only, h);
            figure(h)
            suplabel(['Ref only. Aligned to the appearance of ref. track, n=', num2str(N_sample)],'t');
            
            h = figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*3, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_only_ref_aligned2ref_end, ref_only_cor_aligned2ref_end, track_data.analysis_info, 'ref', lifetime_ref_only, h);
            figure(h)
            suplabel(['Ref only. Aligned to the disappearance of ref. track, n=', num2str(N_sample)],'t');
            
            h = figure('OuterPosition', [1, scrsz(4)-fig_height*4, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(ref_only_ref_aligned2refmax, ref_only_cor_aligned2refmax, track_data.analysis_info, 'ref', lifetime_ref_only, h);
            figure(h)
            suplabel(['Ref only. Aligned to the maximum of ref. track, n=', num2str(N_sample)],'t');
        end
    end
%% Plot unassociated corresponding tracks
    function [] = plotUnassociatedCorr
        if num_cor_only>0
            h = figure('OuterPosition', [fig_width, scrsz(4)-fig_height*4, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(cor_only_ref_aligned2cor_begin, cor_only_cor_aligned2cor_begin, track_data.analysis_info, 'cor', lifetime_cor_only, h);
            figure(h)
            suplabel(['Cor only. Aligned to the appearance of cor. track, n=', num2str(N_sample)],'t');
            
            h = figure('OuterPosition', [fig_width*2, scrsz(4)-fig_height*4, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(cor_only_ref_aligned2cor_end, cor_only_cor_aligned2cor_end, track_data.analysis_info, 'cor', lifetime_cor_only, h);
            figure(h)
            suplabel(['Cor only. Aligned to the disappearance of cor. track, n=', num2str(N_sample)],'t');
            
            h = figure('OuterPosition', [fig_width*3, scrsz(4)-fig_height*4, fig_width,   fig_height]);
            [h, N_sample] = drawAveragePlots(cor_only_ref_aligned2cormax, cor_only_cor_aligned2cormax, track_data.analysis_info, 'cor', lifetime_cor_only, h);
            figure(h)
            suplabel(['Cor only. Aligned to the maximum of cor. track, n=', num2str(N_sample)],'t');
        end
    end
%% Plot statistics of reference tracks associated with corresponding tracks and reference only tracks.
    function [] =  plotAssociatedRef()
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
    end

%% Plot statistics of corresponding tracks associated with reference tracks and corresponding only tracks.
    function [] = plotAssociatedCorr()
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

    %% Compare statistics across experiments, makes a GUI and calls another function
    function [] = compare_stats()
        num_compare = input("How many experiments would you like to compare?");
        experiments = cell(num_compare);
        for n = 1:num_compare
            experiments{n} = uigetdir('Select next experiment');
        end
        if (num_compare == 1)
        end
%         field_number = size(field_positions,1);
%         andarl_gfp = zeros(size(lifetime_associated));
%         andarl_rfp = zeros(size(lifetime_associated));
%         field_start = zeros(size(field_positions));%find starting positions of each field
%         index = 0;
%         for i = 1:field_number
%             field_start(i) = index + 1; %add starting position
%             for j = 1:field_positions(i, 1)
%                 index = index + 1;
%                 andarl_gfp(index, 1) = lifetime_associated(index, 1);
%                 andarl_gfp(index, 2) = i;
%                 andarl_rfp(index, 1) = lifetime_associated(index, 2);
%                 andarl_rfp(index, 2) = i;
%             end
%         end
%         alpha = 0.05;
%         
%         stringa = ['For all samples the Anderson Darling statistics, or chance',...
%             ' that these data \n are drawn from the same distribution is: '];
%         P_GFP = AnDarksamtest(andarl_gfp, alpha);
%         P_RFP = AnDarksamtest(andarl_rfp, alpha);
%         stringb = [stringa, sprintf('for GFP %.4f and for RFP %.4f', P_GFP, P_RFP)];
%         disp(stringb);
%         
%         a = 1:field_number; %generate unique pairs of fields to test
%         b = 1:field_number;
%         [p, q] = meshgrid(a,b);
%         mask = triu(ones(field_number), 1) > 0.5;
%         pairs = [p(mask) q(mask)];
%         
%         
%         for prz = 1:size(pairs,1)
%             current_pair = pairs(prz,:);
%             num1 = current_pair(1);
%             num2 = current_pair(2);
%             current_gfp = [andarl_gfp(field_start(num1):field_start(num1)+field_positions(num1)-1,1); ...
%                 andarl_gfp(field_start(num2):field_start(num2)+field_positions(num2)-1,1)];
%             current_gfp = [current_gfp, ones(size(current_gfp,1),1)];
%             current_gfp(field_positions(num1)+1 : field_positions(num2)...
%                 +field_positions(num1), 2) = 2;
%             
%             current_rfp = [andarl_rfp(field_start(num1):field_start(num1) + field_positions(num1) - 1, 1) ...
%                 ; andarl_rfp(field_start(num2):field_start(num2) + field_positions(num2) - 1, 1)];
%             current_rfp=[current_rfp, ones(size(current_rfp,1),1)];
%             current_rfp(field_positions(num1) + 1 : field_positions(num2)...
%                 + field_positions(num1), 2) = 2;
%             P_GFP = AnDarksamtest(current_gfp, alpha);
%             P_RFP = AnDarksamtest(current_rfp, alpha);
%             stringc = sprintf('For fields %d and %d the AD stat is %.4f for GFP and %.4f for RFP', num1, num2, P_GFP, P_RFP);
%             disp(stringc);
%         end
%         
%         figure
%         hindex = 0;
%         for i = 1:field_number %go between fields
%             h.count = 0; %gfp
%             h.range = [0, max(max(lifetime_associated)) + 1];
%             h.binwidth = 1; %change this to change resolution of histogram
%             h2.count = 0; %rfp
%             h2.range = [0, max(max(lifetime_associated)) + 1];
%             h2.binwidth = 1;
%             for k=1:field_positions(i,1) %advance within each field
%                 hindex = hindex + 1;
%                 h = histo(h, lifetime_associated(hindex, 1)); %bin in histo for gfp
%                 h2 = histo(h2, lifetime_associated(hindex, 2)); %bin in histo2 for rfp
%             end
%             string1=sprintf('Histogram of GFP Lifetimes in Field %d', i);
%             string2=sprintf('Histogram of RFP Lifetimes in Field %d', i);
%             
%             subplot(field_number, 1, i);
%             bar(h.vals, h.hist);
%             title(string1);
%             xlabel('Lifetime in seconds');
%             ylabel('Frequency');
%         end
%         
%         figure
%         hindex = 0;
%         for i = 1:field_number %go between fields
%             h.count = 0; %gfp
%             h.range = [0, max(max(lifetime_associated)) + 1];
%             h.binwidth = 1; %change this to change resolution of histogram
%             h2.count = 0; %rfp
%             h2.range = [0, max(max(lifetime_associated)) + 1];
%             h2.binwidth = 1;
%             for k = 1:field_positions(i, 1) %advance within each field
%                 hindex = hindex + 1;
%                 h = histo(h,lifetime_associated(hindex, 1)); %bin in histo for gfp
%                 h2 = histo(h2,lifetime_associated(hindex, 2)); %bin in histo2 for rfp
%             end
%             string1 = sprintf('Histogram of GFP Lifetimes in Field %d', i);
%             string2 = sprintf('Histogram of RFP Lifetimes in Field %d', i);
%             
%             subplot(field_number, 1, i);
%             bar(h2.vals, h2.hist);
%             title(string2);
%             xlabel('Lifetime in seconds');
%             ylabel('Frequency');
%         end
    end
end
 
 function [h, N_sample] = drawAveragePlots(ref_profile, cor_profile, analysis_info, category, track_lifetimes, h)
 if or(numel(ref_profile)== 0, numel(cor_profile)== 0 )
     error('No tracks are found', ref_profile, cor_profile)
 end
 % h is a figure handler
 % cor_profile and ref_profile are coordinates of reference tracks and corresponding tracks.
 tp_max = analysis_info.tp_max;
 interval = analysis_info.interval;
 
  
    normalize_intensities = 1;
    %  normalize intensity to max within each track, 
%     and subtract bg intensity (min intensity each track)
    reserve_ref_profile = ref_profile;
    reserve_cor_profile = cor_profile;
    if normalize_intensities
        for i = 1:size(ref_profile,1)

            ref_profile(i,:) = (ref_profile(i,:)-nanmin(ref_profile(i,:)))/(nanmax(ref_profile(i,:))-nanmin(ref_profile(i,:)));
        end

        for i = 1:size(cor_profile,1)

            cor_profile(i,:) = (cor_profile(i,:)-nanmin(cor_profile(i,:)))/(nanmax(cor_profile(i,:))-nanmin(cor_profile(i,:)));
        end
    end
 
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
%  subplot(5,2,[1,3]);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_profile(1:10,:)); % MA temporary!!
 subplot(5,2,[1,3]);plot(-(tp_max-1)*interval:interval:tp_max*interval, cor_profile);

 axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval,  min(min(cor_profile)), max(max(cor_profile))]);
 title('cor')
 subplot(5,2,[2,4]);plot(-(tp_max-1)*interval:interval:tp_max*interval, ref_profile);
 axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval,  min(min(ref_profile)), max(max(ref_profile))]);
 title('ref')
 
%  mean and std intensity over time, not normalized
 subplot(5,2,[7,9])

 x = -(tp_max*interval-interval):interval:tp_max*interval;
 y1 = mean(cor_profile_zero); %./2; % MA Temporary!!
 y2 = mean(ref_profile_zero);
 % need the error to be a 2 x 1 x 242
 b = zeros(size(y1, 2), 1, 2);
 b(:, 1, 1) = std(cor_profile_zero); %./2; % MA Temporary!!
 b(:, 1, 2) = std(ref_profile_zero);
 boundedline(x, [y1; y2], b, 'alpha');
 axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval, intensity_min, intensity_max]);
 legend('cor', 'ref')
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
     disp('Category should be one of the following three: associated, ref, cor');
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
 
 title([title_text,'-track lifetimes:', lifetime_text(1:(end-1))])
 axis([(axis_min-tp_max)*interval, (axis_max-tp_max)*interval, intensity_min, intensity_max]);
 
 N_sample = size(cor_profile,1);
 end
