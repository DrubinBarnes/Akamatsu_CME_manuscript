% Read in imaris result file exported as csv.
% Let users determine the shown tracks are good for further analysis.
% The results are saved as text file.
%
%
%      Sun Hae Hong
%      David Drubin lab
%      University of  California, Berkeley
%
%      Copyright 2015
%
% log
% Julian: try to stop from crashing when you don't input 0,1, or 2 
% 4/30/2015
% from        goodflag = input([num2str(i), '-th ref track from' num2str(asso_track_num) ' associated tracks.\nType 1 if for a good track and 0 for a bad track.'...
%             '\n Type 2 for unusual tracks.']);
% to        goodflag = input(['ref track #', num2str(i), '-from' num2str(asso_track_num) '.\nType 1 if for a good track and 0 for a bad track.'...
%             '\n Type 2 for unusual tracks.']);
% Matt Jan/12/16: add third option for special case classification
% Matt: add randomization of tracks to associate. so you don'thave to do
% all. 

% trackPlottingOrder = randperm(nbTracks, nbTracksToPlot);



%%
function Manually_pick_asso_ref_cor20150312(data)

%% choose associated_track file
if (nargin == 0)
    load('previous_file.mat')
    [track_file, track_folder_name] = uigetfile('*.mat','Pick mat file of associated tracks.',previous_file);
    if isdir(track_folder_name)
        previous_file = fullfile(track_folder_name, track_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
else
    track_file = 'associated_tracks_ref_cor.mat';
    track_folder_name = data.source(1:end - 4);
end

%% Load track results
load(fullfile(track_folder_name, track_file));

%% Parameter to use
pixel_size = analysis_info.pixel_size; % in nm
shift_y = analysis_info.shift_y;
shift_x = analysis_info.shift_x;
shift_y_nm = analysis_info.shift_y_nm;
shift_x_nm = analysis_info.shift_x_nm;
ref_channel = analysis_info.ref_channel;
cor_channel = analysis_info.cor_channel;
extra_tp = analysis_info.extra_tp;

%% Read in reference and corresponding raw image.
if exist(analysis_info.cor_image_name, 'file') && exist(analysis_info.ref_image_name, 'file')
    cor_stack_3D = openTiffStack(analysis_info.cor_image_name);
    ref_stack_3D = openTiffStack(analysis_info.ref_image_name);
else
    
    %% choose reference images
    load('previous_file.mat')
    [ref_img_file, ref_img_folder_name] = uigetfile('*.tif','Pick reference movie.',previous_file);
    if isdir(ref_img_folder_name)
        previous_file = fullfile(ref_img_folder_name, ref_img_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
    
    %% choose corresponding images
    load('previous_file.mat')
    [cor_img_file, cor_img_folder_name] = uigetfile('*.tif','Pick corresponding movie.',previous_file);
    if isdir(cor_img_folder_name)
        previous_file = fullfile(cor_img_folder_name, cor_img_file);
        save('previous_file', 'previous_file')
    else
        display('You did not pick a folder.')
        return
    end
    cor_stack_3D = openTiffStack(fullfile(cor_img_folder_name,cor_img_file));
    ref_stack_3D = openTiffStack(fullfile(ref_img_folder_name,ref_img_file));
end

%% Retrieve the number of tracks from the saved results.
asso_track_num = numel(track_stat_ref.ref_tracks_associated);
ref_track_num = numel(track_stat_ref.ref_tracks_unassociated);
cor_track_num = numel(track_stat_cor.cor_tracks_unassociated);



%% Retrieve variables from the saved results
ref_primary_track = track_stat_ref.ref_primary_track;
xcoord_ref = track_stat_ref.xcoord;
ycoord_ref = track_stat_ref.ycoord;
xcoord_cor = track_stat_cor.xcoord;
ycoord_cor = track_stat_cor.ycoord;
% xcoord_cor_aligned = track_stat_cor.xcoord;
% ycoord_cor_aligned = track_stat_cor.ycoord - shift_y_nm;
%% Make a file to save manual decisions
% This step will be skipped if a "manually_cleaned_track_list.txt" file
% already exists in the target folder.

manual_pick_file_name = 'manually_cleaned_track_list.txt';
initialize_flag = false;
if ~exist(fullfile(track_folder_name, manual_pick_file_name), 'file')
    initialize_flag = true;
    
    % Generate a random but defined order to classify tracks
    
    assocTrackClassificationOrder   = randperm(asso_track_num);
    refOnlyTrackClassificationOrder = randperm(ref_track_num);
    corOnlyTrackClassificationOrder = randperm(cor_track_num);
    
else
    %% open output file.
    fid = fopen(fullfile(track_folder_name, manual_pick_file_name), 'r');
    whole_text = textscan(fid,'%s','HeaderLines', 10);
    fclose(fid);
    if numel(whole_text{1}) == 0
        % if none of the previous manual tag was recorded,then initialize
        % the file.
        initialize_flag = true;
        clear whole_text
    else
        
        % Modified by MA so you can go back to classify any group of tracks
        % later.
        
        track_flag_list = dlmread(fullfile(track_folder_name, manual_pick_file_name), '\t', 10, 0);
        %         track_group = track_flag_list(end, 1);
        
        track_flag_list_probe_1_mask = track_flag_list(:,1)==1;
        track_flag_list_probe_2_mask = track_flag_list(:,1)==2;
        track_flag_list_probe_3_mask = track_flag_list(:,1)==3;
        
        track_flag_list_probe_1 = track_flag_list(track_flag_list_probe_1_mask,:);
        track_flag_list_probe_2 = track_flag_list(track_flag_list_probe_2_mask,:);
        track_flag_list_probe_3 = track_flag_list(track_flag_list_probe_3_mask,:);
        
        
        nbTracksClassified_1 = length(track_flag_list_probe_1);
        nbTracksClassified_2 = length(track_flag_list_probe_2);
        nbTracksClassified_3 = length(track_flag_list_probe_3);
        
%         asso_track_num = numel(track_stat_ref.ref_tracks_associated);
% ref_track_num = numel(track_stat_ref.ref_tracks_unassociated);
% cor_track_num = numel(track_stat_cor.cor_tracks_unassociated);
%         
        
        track_group = 1;
        
        % MA modify so that the length not the
        % precise value is compared. (when you're opening an existing
        % manually_cleaned_list.txt file
        
        %  MAKE SURE THIS WORKS
        
        if  (nbTracksClassified_1 < asso_track_num)
            first_track_index_1 = nbTracksClassified_1+1;
        elseif (nbTracksClassified_1 == asso_track_num && nbTracksClassified_2 == 0)
%             first_track = 1;
            track_group = 2;
        end
        if (nbTracksClassified_2 < ref_track_num)
            first_track_index_2 = nbTracksClassified_2+1;
        elseif (nbTracksClassified_2 == ref_track_num && nbTracksClassified_3 == 0)
%             first_track = 1;
            track_group = 3;
        end
        if (nbTracksClassified_3 < cor_track_num)
            first_track_index_3 = nbTracksClassified_3 + 1;
        elseif (nbTracksClassified_3 == cor_track_num)
            %first_track = find(track_stat_cor.cor_tracks_unassociated == track_flag_list_probe_3(end, 2)) + 1;
                         
        end
        
        if (nbTracksClassified_1 == asso_track_num && nbTracksClassified_2 == ref_track_num && nbTracksClassified_3 == cor_track_num)
            display('Manual selection is already done.')
                         return
        end
        
%         if (nbTracksClassified_2 < ref_track_num)
%             first_track_index_2 = nbTracksClassified_2+1;
%         end
%         if (nbTracksClassified_3 < cor_track_num)
%             first_track_index_3 = nbTracksClassified_3 + 1;
%         end
%         if  (track_flag_list_probe_1(end, 2) < track_stat_ref.ref_tracks_associated(end))
%             first_track = find(track_stat_ref.ref_tracks_associated == track_flag_list_probe_1(end, 2)) + 1;
%         elseif (track_flag_list_probe_1(end, 2) == track_stat_ref.ref_tracks_associated(end))
%             first_track = 1;
%             track_group = 2;
%         elseif (track_flag_list_probe_2(end, 2) < track_stat_ref.ref_tracks_unassociated(end))
%             first_track = find(track_stat_ref.ref_tracks_unassociated == track_flag_list_probe_2(end, 2)) + 1;
%         elseif (track_flag_list_probe_2(end, 2) == track_stat_ref.ref_tracks_unassociated(end))
%             first_track = 1;
%             track_group = 3;
%         elseif (track_flag_list_probe_3(end, 2) < track_stat_cor.cor_tracks_unassociated(end))
%             first_track = find(track_stat_cor.cor_tracks_unassociated == track_flag_list_probe_3(end, 2)) + 1;
%         elseif (track_flag_list_probe_3(end, 2) == track_stat_cor.cor_tracks_unassociated(end))
%             first_track = find(track_stat_cor.cor_tracks_unassociated == track_flag_list_probe_3(end, 2)) + 1;
%             %             display('Manual selection is already done.')
%             %             return
%         end
        
    end
    %fid = fopen(fullfile(track_folder_name, manual_pick_file_name), 'a');
end
if initialize_flag
    %% initialize an output file
    fid = fopen(fullfile(track_folder_name, manual_pick_file_name), 'w');
    fprintf(fid, '%s\n%s\n%s\n%s\n\n', ['ref csv file:', analysis_info.ref_csv_name],...
        ['cor csv file:', analysis_info.cor_csv_name],...
        ['ref image:', analysis_info.ref_image_name],...
        ['cor image:', analysis_info.cor_image_name]);
    fprintf(fid, '%s\t%s\t%s\n', '1 for ref, 2 for cor', 'Track number', '1 for good, 0 for bad');
    fprintf(fid,'%s\n', 'Column 1: Track group. 1 for associated tracks, 2 for unassociated reference tracks, 3 for unassociated corresponding tracks');
    fprintf(fid,'%s\n', 'Column 2: Track number.');
    fprintf(fid,'%s\n', 'Column 3: Decision. 1 for good tracks to be further analyzed. 0 for bad tracks that are discarded.');
    fprintf(fid, '%s\t%s\t%s\n',...
        'Track group', 'Track number', 'Decision');
    fclose(fid)
    %% Track category is 1 for reference tracks whether they have associated corresponding tracks.
    %% Track category is 1 for corresponding tracks that are not associated with reference tracks.
    track_group = 1;
%     first_track = 1; % this may be obsolete now
    first_track_index_1 = 1;
    first_track_index_2 = 1;
    first_track_index_3 = 1;

    ii = 1;
    jj = 1;
    kk = 1;
    
    % group numbers to identify tracks as assoc, cor only, or ref only
    track_flag_list = NaN*ones(asso_track_num + ref_track_num + cor_track_num, 3);
    track_flag_list(1:asso_track_num,1) = 1;
    track_flag_list((asso_track_num + 1):(asso_track_num + ref_track_num),1) = 2;
    track_flag_list((asso_track_num  + ref_track_num +1):(asso_track_num + ref_track_num + cor_track_num),1) = 3;

end

asso_track_num = numel(track_stat_ref.ref_tracks_associated);
ref_track_num = numel(track_stat_ref.ref_tracks_unassociated);
cor_track_num = numel(track_stat_cor.cor_tracks_unassociated);


if track_group == 1
%     for i = first_track:asso_track_num

    

    for i = assocTrackClassificationOrder(first_track_index_1:end) % will need to have this start some way through, e.g. indexed baesd on number of events I've classified

        ref_track = track_stat_ref.ref_tracks_associated(i);
        profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
        
        % profile_length needs to be at least 3 time points.
        if profile_length < 3
            clean_track_flag =  false;
            track_flag_list(i, 2:3) = [ref_track, clean_track_flag];
            %             track_flag_list(i, 4) = 2;
            %             fprintf(fid, '%d\t%d\t%d\n', 1, ref_track, 0);
            dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [1, ref_track, 0], '-append', 'delimiter', '\t');
            continue
        end
        % If part of data is missing for any reason, the track is not
        % considered.
        if any(any(isnan(associated_tracks(ref_track, 1).cor_montage))) ||any(any(isnan(associated_tracks(ref_track, 1).ref_montage)))
            clean_track_flag =  false;
            track_flag_list(i, 2:3) = [ref_track, clean_track_flag];
            dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [1, ref_track, 0], '-append', 'delimiter', '\t');
            continue
        end
      
        
        %% Plot trajectories
        figure('Position', [35,185,1828,695]);
        j = 0;
        for ref_track2 = ROI_track_ref(ref_track, 1).associated_ref_track
            j = j+1;
            tp_sequence = ROI_track_ref(ref_track2, 1).tp_first:ROI_track_ref(ref_track2, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                ref_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)/2];
            else
                ref_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)];
            end
            subplot(6,5,[1,6]);
            plot(xcoord_ref(tp_sequence, ref_track2) + shift_x_nm, ycoord_ref(tp_sequence,ref_track2) + shift_y_nm,...
                'Color', ref_color);hold on
            plot(xcoord_ref(tp_sequence(1), ref_track2) + shift_x_nm, ycoord_ref(tp_sequence(1),ref_track2) + shift_y_nm,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot(xcoord_ref(tp_sequence(end), ref_track2) + shift_x_nm, ycoord_ref(tp_sequence(end),ref_track2) + shift_y_nm,...
                '.','MarkerSize', 10, 'Color', ref_color);
            
            subplot(6,5,[2,7]);
            imshow(fit_contrast(double(ref_stack_3D(:,:,tp_sequence(1))))); hold on;

            plot((xcoord_ref(tp_sequence, ref_track2))/pixel_size,...
                (ycoord_ref(tp_sequence,ref_track2))/pixel_size,...
                'Color', ref_color);hold on
            plot((xcoord_ref(tp_sequence(1), ref_track2))/pixel_size,...
                (ycoord_ref(tp_sequence(1),ref_track2))/pixel_size,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot((xcoord_ref(tp_sequence(end), ref_track2))/pixel_size,...
                (ycoord_ref(tp_sequence(end),ref_track2))/pixel_size,...
                '.','MarkerSize', 10, 'Color', ref_color);
%             axis([0 100 0 100]);
           % show the region near your track
            
           meanX = nanmean(xcoord_ref(tp_sequence, ref_track2)/pixel_size);
           meanY = nanmean(ycoord_ref(tp_sequence,ref_track2)/pixel_size);
           
            axis([meanX-30, meanX+30, meanY-30, meanY+30]);

            
            title(['ref-', analysis_info.ref_channel '-image, tp=', num2str(tp_sequence(1))]);
            
            subplot(6,5,[3:5]);
            plot(tp_sequence, xcoord_ref(tp_sequence, ref_track2) + shift_x_nm,...
                'Color', ref_color);
            hold on
            plot(tp_sequence(1), xcoord_ref(tp_sequence(1), ref_track2) + shift_x_nm,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot(tp_sequence(end), xcoord_ref(tp_sequence(end), ref_track2) + shift_x_nm,...
                '.','MarkerSize', 10, 'Color', ref_color);
            
            subplot(6,5,[8:10]);
            plot(tp_sequence, ycoord_ref(tp_sequence,ref_track2) + shift_y_nm,'Color', ref_color);
            hold on
            plot(tp_sequence(1), ycoord_ref(tp_sequence(1), ref_track2) + shift_y_nm,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot(tp_sequence(end), ycoord_ref(tp_sequence(end), ref_track2) + shift_y_nm,...
                '.','MarkerSize', 10, 'Color', ref_color);
        end
        j = 0;
        for cor_track = ROI_track_ref(ref_track, 1).associated_cor_track
            j = j+1;
            tp_sequence = ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                cor_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)];
            else
                cor_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)/2];
            end
            subplot(6,5,[1,6]);
            plot(xcoord_cor(tp_sequence, cor_track), ycoord_cor(tp_sequence,cor_track), 'Color', cor_color);
            hold on
            plot(xcoord_cor(tp_sequence(1), cor_track), ycoord_cor(tp_sequence(1),cor_track),...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot(xcoord_cor(tp_sequence(end), cor_track), ycoord_cor(tp_sequence(end),cor_track),...
                '.','MarkerSize', 10, 'Color', cor_color);
            subplot(6,5,[2,7]);
            plot((xcoord_cor(tp_sequence, cor_track)-shift_x_nm)/pixel_size,...
                (ycoord_cor(tp_sequence,cor_track)-shift_y_nm)/pixel_size, 'Color', cor_color);
            hold on
            plot((xcoord_cor(tp_sequence(1), cor_track)-shift_x_nm)/pixel_size,...
                (ycoord_cor(tp_sequence(1),cor_track)-shift_y_nm)/pixel_size,...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot((xcoord_cor(tp_sequence(end), cor_track)-shift_x_nm)/pixel_size,...
                (ycoord_cor(tp_sequence(end),cor_track)-shift_y_nm)/pixel_size,...
                '.','MarkerSize', 10, 'Color', cor_color);
            
            
            subplot(6,5,[3:5]);
            plot(tp_sequence, xcoord_cor(tp_sequence, cor_track), 'Color', cor_color);
            hold on
            plot(tp_sequence(1), xcoord_cor(tp_sequence(1), cor_track),...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot(tp_sequence(end), xcoord_cor(tp_sequence(end), cor_track),...
                '.','MarkerSize', 10, 'Color', cor_color);
            subplot(6,5,[8:10]);
            plot(tp_sequence, ycoord_cor(tp_sequence,cor_track),'Color', cor_color);
            hold on
            plot(tp_sequence(1), ycoord_cor(tp_sequence(1), cor_track),...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot(tp_sequence(end), ycoord_cor(tp_sequence(end), cor_track),...
                '.','MarkerSize', 10, 'Color', cor_color);
        end
        subplot(6,5, [1,6]); axis equal; xlabel('x');ylabel('y');
        subplot(6,5, [3:5]); xlabel('time');ylabel('x');
        subplot(6,5, [8:10]); xlabel('time');ylabel('y');
        %%
        subplot(6,5, [11:15]); imshow(fit_contrast(associated_tracks(ref_track, 1).ref_movie));
        hold on
        ROI_x = ROI_track_ref(ref_track,1).ROI_x;
        ROI_y = ROI_track_ref(ref_track,1).ROI_y;
        
        ROI_delx = ROI_x(1,2) - ROI_x(1,1) + 1;
        ROI_dely = ROI_y(1,2) - ROI_y(1,1) + 1;
        j = 0;
        for ref_track2 = ROI_track_ref(ref_track, 1).associated_ref_track
            j = j+1;
            tp_sequence = ROI_track_ref(ref_track2, 1).tp_first:ROI_track_ref(ref_track2, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                ref_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)/2];
            else
                ref_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)];
            end
            
            
            ref_x = xcoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_x(1,1) + 1;
            ref_y = ycoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_y(1,1) + 1;
            del_x_offset = [(ROI_track_ref(ref_track2).tp_first - ROI_track_ref(ref_track).tp_prev_all):(ROI_track_ref(ref_track2).tp_last - ROI_track_ref(ref_track).tp_prev_all)]'*ROI_delx;
            
            plot(ref_x + del_x_offset, ref_y, 'x', 'Color', ref_color)
            
        end
        subplot(6,5, [16:20]); imshow(fit_contrast(associated_tracks(ref_track, 1).cor_movie));
        hold on
        j = 0;
        for cor_track = ROI_track_ref(ref_track, 1).associated_cor_track
            j = j+1;
            tp_sequence = ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                cor_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)];
            else
                cor_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)/2];
            end
            
            cor_x = xcoord_cor(tp_sequence, cor_track)/pixel_size -ROI_track_ref(ref_track, 1).ROI_x(1,1) - shift_x + 1;
            cor_y = ycoord_cor(tp_sequence, cor_track)/pixel_size -ROI_track_ref(ref_track, 1).ROI_y(1,1) - shift_y+ 1;
            del_x_offset = [(ROI_track_cor(cor_track).tp_first - ROI_track_ref(ref_track).tp_prev_all):(ROI_track_cor(cor_track).tp_last - ROI_track_ref(ref_track).tp_prev_all)]'*ROI_delx;
            
            plot(cor_x + del_x_offset, cor_y, 'x','Color', cor_color)
        end
        %% if no corresponding track is associated to the reference track, show "o" at the corresponding
        %% position
        if numel(ROI_track_ref(ref_track, 1).associated_cor_track) == 0
            j = 0;
            for ref_track2 = ROI_track_ref(ref_track, 1).associated_ref_track
                j = j+1;
                tp_sequence = ROI_track_ref(ref_track2, 1).tp_first:ROI_track_ref(ref_track2, 1).tp_last;
                if strcmp(ref_channel, 'RFP')
                    ref_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)/2];
                else
                    ref_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)];
                end
                
                ref_x = xcoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_x(1,1) + 1;
                ref_y = ycoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_y(1,1) + 1;
                del_x_offset = [(ROI_track_ref(ref_track2).tp_first - ROI_track_ref(ref_track).tp_prev_all):(ROI_track_ref(ref_track2).tp_last - ROI_track_ref(ref_track).tp_prev_all)]'*ROI_delx;
                
                %plot(ref_x + del_x_offset + shift_x, ref_y + shift_y, 'o', 'Color', [0 0.1 0.8])
                plot(ref_x + del_x_offset, ref_y, 'o', 'Color', [0 0.1 0.8])
            end
        end
        subplot(6,5, [21:25]);
        imshow(fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage));
        title('ref montage')
        if ~any(any(isnan(associated_tracks(ref_track, 1).cor_montage)))
            subplot(6,5, [26:30]); imshow(fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage));title('cor montage')
        end
        suplabel(['ref track number:', num2str(ref_track)]);
        set(gcf,'color','w');
        
        %%
        clean_track_flag =7;
        while ~any(clean_track_flag == [0,1,2,3,99])
        
        clean_track_flag = input(['ref track ', num2str(i), ', #' num2str(ii) ' of ' num2str(asso_track_num) '.\nType 1 if for a good track and 0 for a bad track.'...
            '\n Type 2 for unusual tracks. \n Type 3 for special cases. \n']);
        
            % prevent no character entry from breaking program 
            if ~isa(clean_track_flag,'numeric') || isempty(clean_track_flag) 
                display('The answer needs to be either 0 1 2 or 3.')
                clean_track_flag = 7; % some number to force you to go again
            end
        
        
            %display('The answer needs to be either 0 or 1.')
            %clean_track_flag = input('Type 1 if for a good track and 0 for a bad track.');
        end
        
        
        % Allow you to exit this stage of the analysis early and move on to the next
        % track group (MA) if you type 99
        
        if clean_track_flag == 99
            
            disp('Moving on to unassociated ref tracks.')
   
            break
            
        end
        
        track_flag_list(i, :) = [1, ref_track, clean_track_flag];
        %         fprintf(fid, '%d\t%d\t%d\n', 1, ref_track, clean_track_flag);
        dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [1, ref_track, clean_track_flag], '-append', 'delimiter', '\t');
        track_flag_list(i, 2:3) = [ref_track, clean_track_flag];
        save(fullfile(track_folder_name, 'manually_cleaned_track_list.mat'), 'track_flag_list');
        close all
        
        ii = ii + 1; % index of tracks (1:end) (independent of track order)
        
    end
    track_group = 2;
%     first_track = 1; % will commenting this make a bug?
end
if track_group == 2
    
%     jj = 1;
%     for i = first_track:ref_track_num
    for i = refOnlyTrackClassificationOrder(first_track_index_2:end)
        ref_track = track_stat_ref.ref_tracks_unassociated(i);
        profile_length = ROI_track_ref(ref_track, 1).tp_post_all - ROI_track_ref(ref_track, 1).tp_prev_all + 1;
        
        % profile_length needs to be at least 3 time points.
        if profile_length < 3
            clean_track_flag =  false;
            track_flag_list(i, :) = [2, ref_track, clean_track_flag];
            %             track_flag_list(i, 4) = 2;
            %             fprintf(fid, '%d\t%d\t%d\n', 2, ref_track, 0);
            dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [2, ref_track, 0], '-append', 'delimiter', '\t');
            continue
        end
        
        % If part of data is missing for any reason, the track is not
        % considered.
        if any(any(isnan(associated_tracks(ref_track, 1).cor_montage))) ||any(any(isnan(associated_tracks(ref_track, 1).ref_montage)))
            clean_track_flag =  false;
            track_flag_list(i, :) = [2,ref_track, clean_track_flag];
            dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [2, ref_track, 0], '-append', 'delimiter', '\t');
            continue
        end

        %% Plot trajectories
        figure('Position', [35,185,1828,695]);
        j = 0;
        for ref_track2 = ROI_track_ref(ref_track, 1).associated_ref_track
            j = j+1;
            tp_sequence = ROI_track_ref(ref_track2, 1).tp_first:ROI_track_ref(ref_track2, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                ref_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)/2];
            else
                ref_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)];
            end
            subplot(6,5,[1,6]);
            plot(xcoord_ref(tp_sequence, ref_track2) + shift_x_nm, ycoord_ref(tp_sequence,ref_track2) + shift_y_nm,...
                'Color', ref_color);hold on
            plot(xcoord_ref(tp_sequence(1), ref_track2) + shift_x_nm, ycoord_ref(tp_sequence(1),ref_track2) + shift_y_nm,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot(xcoord_ref(tp_sequence(end), ref_track2) + shift_x_nm, ycoord_ref(tp_sequence(end),ref_track2) + shift_y_nm,...
                '.','MarkerSize', 10, 'Color', ref_color);
            
            subplot(6,5,[2,7]);
            imshow(fit_contrast(double(ref_stack_3D(:,:,tp_sequence(1))))); hold on;
            plot((xcoord_ref(tp_sequence, ref_track2))/pixel_size,...
                (ycoord_ref(tp_sequence,ref_track2))/pixel_size,...
                'Color', ref_color);hold on
            plot((xcoord_ref(tp_sequence(1), ref_track2))/pixel_size,...
                (ycoord_ref(tp_sequence(1),ref_track2))/pixel_size,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot((xcoord_ref(tp_sequence(end), ref_track2))/pixel_size,...
                (ycoord_ref(tp_sequence(end),ref_track2))/pixel_size,...
                '.','MarkerSize', 10, 'Color', ref_color);
%             axis([0 100 0 100]);
            
            meanX = nanmean(xcoord_ref(tp_sequence, ref_track2)/pixel_size);
            meanY = nanmean(ycoord_ref(tp_sequence,ref_track2)/pixel_size);
           
            axis([meanX-30, meanX+30, meanY-30, meanY+30]);
            
            title(['ref-', analysis_info.ref_channel '-image, tp=', num2str(tp_sequence(1))]);
            
            subplot(6,5,[3:5]);
            plot(tp_sequence, xcoord_ref(tp_sequence, ref_track2) + shift_x_nm,...
                'Color', ref_color);
            hold on
            plot(tp_sequence(1), xcoord_ref(tp_sequence(1), ref_track2) + shift_x_nm,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot(tp_sequence(end), xcoord_ref(tp_sequence(end), ref_track2) + shift_x_nm,...
                '.','MarkerSize', 10, 'Color', ref_color);
            
            subplot(6,5,[8:10]);
            plot(tp_sequence, ycoord_ref(tp_sequence,ref_track2) + shift_y_nm,'Color', ref_color);
            hold on
            plot(tp_sequence(1), ycoord_ref(tp_sequence(1), ref_track2) + shift_y_nm,...
                'x','MarkerSize', 10, 'Color', ref_color);
            plot(tp_sequence(end), ycoord_ref(tp_sequence(end), ref_track2) + shift_y_nm,...
                '.','MarkerSize', 10, 'Color', ref_color);
        end
        j = 0;
        for cor_track = ROI_track_ref(ref_track, 1).associated_cor_track
            j = j+1;
            tp_sequence = ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                cor_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)];
            else
                cor_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)/2];
            end
            subplot(6,5,[1,6]);
            plot(xcoord_cor(tp_sequence, cor_track), ycoord_cor(tp_sequence,cor_track), 'Color', cor_color);
            hold on
            plot(xcoord_cor(tp_sequence(1), cor_track), ycoord_cor(tp_sequence(1),cor_track),...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot(xcoord_cor(tp_sequence(end), cor_track), ycoord_cor(tp_sequence(end),cor_track),...
                '.','MarkerSize', 10, 'Color', cor_color);
            subplot(6,5,[2,7]);
            plot((xcoord_cor(tp_sequence, cor_track)-shift_x_nm)/pixel_size,...
                (ycoord_cor(tp_sequence,cor_track)-shift_y_nm)/pixel_size, 'Color', cor_color);
            hold on
            plot((xcoord_cor(tp_sequence(1), cor_track)-shift_x_nm)/pixel_size,...
                (ycoord_cor(tp_sequence(1),cor_track)-shift_y_nm)/pixel_size,...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot((xcoord_cor(tp_sequence(end), cor_track)-shift_x_nm)/pixel_size,...
                (ycoord_cor(tp_sequence(end),cor_track)-shift_y_nm)/pixel_size,...
                '.','MarkerSize', 10, 'Color', cor_color);
            
            
            subplot(6,5,[3:5]);
            plot(tp_sequence, xcoord_cor(tp_sequence, cor_track), 'Color', cor_color);
            hold on
            plot(tp_sequence(1), xcoord_cor(tp_sequence(1), cor_track),...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot(tp_sequence(end), xcoord_cor(tp_sequence(end), cor_track),...
                '.','MarkerSize', 10, 'Color', cor_color);
            subplot(6,5,[8:10]);
            plot(tp_sequence, ycoord_cor(tp_sequence,cor_track),'Color', cor_color);
            hold on
            plot(tp_sequence(1), ycoord_cor(tp_sequence(1), cor_track),...
                'x','MarkerSize', 10, 'Color', cor_color);
            plot(tp_sequence(end), ycoord_cor(tp_sequence(end), cor_track),...
                '.','MarkerSize', 10, 'Color', cor_color);
        end
        subplot(6,5, [1,6]); axis equal; xlabel('x');ylabel('y');
        subplot(6,5, [3:5]); xlabel('time');ylabel('x');
        subplot(6,5, [8:10]); xlabel('time');ylabel('y');
        %%
        subplot(6,5, [11:15]); imshow(fit_contrast(associated_tracks(ref_track, 1).ref_movie));
        hold on
        ROI_x = ROI_track_ref(ref_track,1).ROI_x;
        ROI_y = ROI_track_ref(ref_track,1).ROI_y;
        
        ROI_delx = ROI_x(1,2) - ROI_x(1,1) + 1;
        ROI_dely = ROI_y(1,2) - ROI_y(1,1) + 1;
        j = 0;
        for ref_track2 = ROI_track_ref(ref_track, 1).associated_ref_track
            j = j+1;
            tp_sequence = ROI_track_ref(ref_track2, 1).tp_first:ROI_track_ref(ref_track2, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                ref_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)/2];
            else
                ref_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)];
            end
            
            
            ref_x = xcoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_x(1,1) + 1;
            ref_y = ycoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_y(1,1) + 1;
            del_x_offset = [(ROI_track_ref(ref_track2).tp_first - ROI_track_ref(ref_track).tp_prev_all):(ROI_track_ref(ref_track2).tp_last - ROI_track_ref(ref_track).tp_prev_all)]'*ROI_delx;
            
            plot(ref_x + del_x_offset, ref_y, 'x', 'Color', ref_color)
            
        end
        subplot(6,5, [16:20]); imshow(fit_contrast(associated_tracks(ref_track, 1).cor_movie));
        hold on
        j = 0;
        for cor_track = ROI_track_ref(ref_track, 1).associated_cor_track
            j = j+1;
            tp_sequence = ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last;
            if strcmp(ref_channel, 'RFP')
                cor_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)];
            else
                cor_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_cor_track)/2];
            end
            
            cor_x = xcoord_cor(tp_sequence, cor_track)/pixel_size -ROI_track_ref(ref_track, 1).ROI_x(1,1) - shift_x + 1;
            cor_y = ycoord_cor(tp_sequence, cor_track)/pixel_size -ROI_track_ref(ref_track, 1).ROI_y(1,1) - shift_y+ 1;
            del_x_offset = [(ROI_track_cor(cor_track).tp_first - ROI_track_ref(ref_track).tp_prev_all):(ROI_track_cor(cor_track).tp_last - ROI_track_ref(ref_track).tp_prev_all)]'*ROI_delx;
            
            plot(cor_x + del_x_offset, cor_y, 'x','Color', cor_color)
        end
        %% if no corresponding track is associated to the reference track, show "o" at the corresponding
        %% position
        if numel(ROI_track_ref(ref_track, 1).associated_cor_track) == 0
            j = 0;
            for ref_track2 = ROI_track_ref(ref_track, 1).associated_ref_track
                j = j+1;
                tp_sequence = ROI_track_ref(ref_track2, 1).tp_first:ROI_track_ref(ref_track2, 1).tp_last;
                if strcmp(ref_channel, 'RFP')
                    ref_color = [1, 0, j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)/2];
                else
                    ref_color = [0, 0.9, 1-j/numel(ROI_track_ref(ref_track, 1).associated_ref_track)];
                end
                
                ref_x = xcoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_x(1,1) + 1;
                ref_y = ycoord_ref(tp_sequence, ref_track2)/pixel_size -ROI_track_ref(ref_track, 1).ROI_y(1,1) + 1;
                del_x_offset = [(ROI_track_ref(ref_track2).tp_first - ROI_track_ref(ref_track).tp_prev_all):(ROI_track_ref(ref_track2).tp_last - ROI_track_ref(ref_track).tp_prev_all)]'*ROI_delx;
                
                %plot(ref_x + del_x_offset + shift_x, ref_y + shift_y, 'o', 'Color', [0 0.1 0.8])
                plot(ref_x + del_x_offset, ref_y, 'o', 'Color', [0 0.1 0.8])
            end
        end
        subplot(6,5, [21:25]);
        imshow(fit_contrast_ignore_zero(associated_tracks(ref_track, 1).ref_montage));
        title('ref montage')
        if ~any(any(isnan(associated_tracks(ref_track, 1).cor_montage)))
            subplot(6,5, [26:30]); imshow(fit_contrast_ignore_zero(associated_tracks(ref_track, 1).cor_montage));title('cor montage')
        end
        suplabel(['ref track number:', num2str(ref_track)]);
        set(gcf,'color','w');
        
        %%
        clean_track_flag = 7;
        while ~any(clean_track_flag == [0,1,2,3,99])
        
        clean_track_flag = input(['ref track ', num2str(i), ', #' num2str(jj) ' of ' num2str(ref_track_num) '.\nType 1 if for a good track and 0 for a bad track.'...
            '\n Type 2 for unusual tracks. \n Type 3 for special cases. \n']);
        
            % prevent no character entry from breaking program 
            if ~isa(clean_track_flag,'numeric') || isempty(clean_track_flag) 
                display('The answer needs to be either 0 1 2 or 3.')
                clean_track_flag = 7; % some number to force you to go again
            end
        
        
            %display('The answer needs to be either 0 or 1.')
            %clean_track_flag = input('Type 1 if for a good track and 0 for a bad track.');
        end
        % commenting out old code
        %{
        clean_track_flag = input(['ref track #', num2str(i), '-from' num2str(ref_track_num) '.\nType 1 if for a good track and 0 for a bad track.'...
            '\n Type 2 for unusual tracks.']);
        while ~any(clean_track_flag == [0,1,2])
            display('The answer needs to be either 0 or 1.')
            clean_track_flag = input('Type 1 if for a good track and 0 for a bad track.');
        end
        %}
        %         fprintf(fid, '%d\t%d\t%d\n', 2, ref_track, clean_track_flag);
        
          % Allow you to exit this stage of the analysis early and move on to the next
        % track group (MA) if you type 99
        
        if clean_track_flag == 99
            disp('Moving on to unassociated cor tracks.')
    
            break
            
        end
                
        track_flag_list(i, :) = [2, ref_track, clean_track_flag];
        dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [2, ref_track, clean_track_flag], '-append', 'delimiter', '\t');
        save(fullfile(track_folder_name, 'manually_cleaned_track_list.mat'), 'track_flag_list');
        close all
        
        jj = jj+1;
    end
    track_group = 3; % move on to next set of categories
%     first_track = 1; % will this override my attempt to start at a later track number for the next category (if i've already classified some)? 
end

if track_group ==3
    
%     kk = 1;
    
    cor_track_num = numel(track_stat_cor.cor_tracks_unassociated);
%     for i = first_track:cor_track_num
     for i = corOnlyTrackClassificationOrder(first_track_index_3:end)
  
        cor_track = track_stat_cor.cor_tracks_unassociated(i);
        
        if any(any(isnan(unassociated_tracks(cor_track, 1).ref_montage))) || any(any(isnan(unassociated_tracks(cor_track, 1).cor_montage)))
            clean_track_flag =  false;
            %             fprintf(fid, '%d\t%d\t%d\n', 3, cor_track, clean_track_flag);
            %             fprintf(fid, '%d\t%d\t%d\n', 3, cor_track, 0);
            dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [3, cor_track, 0], '-append', 'delimiter', '\t');
            continue
        end
        
        
        ROI_x = ROI_track_cor(cor_track,1).ROI_x;
        ROI_y = ROI_track_cor(cor_track,1).ROI_y;
        
        ROI_delx = ROI_x(1,2) - ROI_x(1,1) + 1;
        ROI_dely = ROI_y(1,2) - ROI_y(1,1) + 1;
        
        
        %% Plot trajectories
        figure('Position', [35,185,1828,695]);
        if strcmp(ref_channel, 'RFP')
            cor_color = [0, 0.8, 0.1];
        else
            cor_color = [1, 0, 0.1];
        end
        
        tp_sequence = ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last;
        subplot(6,5,[1,6]);
        plot(xcoord_cor(tp_sequence, cor_track), ycoord_cor(tp_sequence,cor_track), 'Color', cor_color);
        hold on
        plot(xcoord_cor(tp_sequence(1), cor_track), ycoord_cor(tp_sequence(1),cor_track),...
            'x','MarkerSize', 10, 'Color', cor_color);
        plot(xcoord_cor(tp_sequence(end), cor_track), ycoord_cor(tp_sequence(end),cor_track),...
            '.','MarkerSize', 10, 'Color', cor_color);
        subplot(6,5,[2,7]);
        imshow(fit_contrast(double(cor_stack_3D(:,:,tp_sequence(1))))); hold on;
        plot(xcoord_cor(tp_sequence, cor_track)/pixel_size, ycoord_cor(tp_sequence,cor_track)/pixel_size, 'Color', cor_color);
        hold on
        plot(xcoord_cor(tp_sequence(1), cor_track), ycoord_cor(tp_sequence(1),cor_track),...
            'x','MarkerSize', 10, 'Color', cor_color);
        plot(xcoord_cor(tp_sequence(end), cor_track), ycoord_cor(tp_sequence(end),cor_track),...
            '.','MarkerSize', 10, 'Color', cor_color);
        
        meanX = nanmean(xcoord_cor(tp_sequence, cor_track)/pixel_size);
        meanY = nanmean(xcoord_cor(tp_sequence, cor_track)/pixel_size);
        
        axis([meanX-30, meanX+30, meanY-30, meanY+30]);
        
        title(['cor-', analysis_info.cor_channel '-image, tp=', num2str(tp_sequence(1))]);
        
        subplot(6,5,[3:5]);
        plot(tp_sequence, xcoord_cor(tp_sequence, cor_track), 'Color', cor_color);
        hold on
        plot(tp_sequence(1), xcoord_cor(tp_sequence(1), cor_track),...
            'x','MarkerSize', 10, 'Color', cor_color);
        plot(tp_sequence(end), xcoord_cor(tp_sequence(end), cor_track),...
            '.','MarkerSize', 10, 'Color', cor_color);
        subplot(6,5,[8:10]);
        plot(tp_sequence, ycoord_cor(tp_sequence,cor_track),'Color', cor_color);
        hold on
        plot(tp_sequence(1), ycoord_cor(tp_sequence(1), cor_track),...
            'x','MarkerSize', 10, 'Color', cor_color);
        plot(tp_sequence(end), ycoord_cor(tp_sequence(end), cor_track),...
            '.','MarkerSize', 10, 'Color', cor_color);
        subplot(6,5, [1,6]); axis equal; xlabel('x');ylabel('y');
        subplot(6,5, [3:5]); xlabel('time');ylabel('x');
        subplot(6,5, [8:10]); xlabel('time');ylabel('y');
        %%
        subplot(6,5, [11:15]); imshow(fit_contrast(unassociated_tracks(cor_track, 1).ref_movie));
        hold on
        subplot(6,5, [16:20]); imshow(fit_contrast(unassociated_tracks(cor_track, 1).cor_movie));
        hold on
        tp_sequence = ROI_track_cor(cor_track, 1).tp_first:ROI_track_cor(cor_track, 1).tp_last;
        
        cor_x = xcoord_cor(tp_sequence, cor_track)/pixel_size -ROI_track_cor(cor_track, 1).ROI_x(1,1) + 1;
        cor_y = ycoord_cor(tp_sequence, cor_track)/pixel_size -ROI_track_cor(cor_track, 1).ROI_y(1,1) + 1;
        del_x_offset = [(ROI_track_cor(cor_track).tp_first - ROI_track_cor(cor_track).tp_prev):(ROI_track_cor(cor_track).tp_last - ROI_track_cor(cor_track).tp_prev)]'*ROI_delx;
        
        plot(cor_x + del_x_offset, cor_y, 'x','Color', cor_color)
        
        subplot(6,5, [21:25]);
        imshow(fit_contrast_ignore_zero(unassociated_tracks(cor_track, 1).ref_montage));
        title('ref montage')
        subplot(6,5, [26:30]); imshow(fit_contrast_ignore_zero(unassociated_tracks(cor_track, 1).cor_montage));
        title('cor montage')
        suplabel(['cor track number:', num2str(cor_track)]);
        set(gcf,'color','w');
        
        %%
        clean_track_flag = 7;
        while ~any(clean_track_flag == [0,1,2,3,99])
        clean_track_flag = input(['cor track ',  num2str(i), ', #' num2str(kk) ' of ' num2str(cor_track_num) '.\nType 1 if for a good track and 0 for a bad track.'...
            '\n Type 2 for unusual tracks. \n Type 3 for special cases. \n']);
        
            % prevent no character entry from breaking program 
            if ~isa(clean_track_flag,'numeric') || isempty(clean_track_flag) 
                display('The answer needs to be either 0 1 2 or 3.')
                clean_track_flag = 7; % some number to force you to go again
            end
        
        
            %display('The answer needs to be either 0 or 1.')
            %clean_track_flag = input('Type 1 if for a good track and 0 for a bad track.');
        end
        
        % commeting out old code
        %{
        clean_track_flag = input(['cor track #', num2str(i), '-from' num2str(cor_track_num) '.\nType 1 if for a good track and 0 for a bad track.']);
        while clean_track_flag ~= 0 && clean_track_flag ~= 1
            
            
            
            display('The answer needs to be either 0 or 1.')
            clean_track_flag = input('Type 1 if for a good track and 0 for a bad track.');
        end
        %}
        %         fprintf(fid, '%d\t%d\t%d\n', 3, cor_track, clean_track_flag);
        
          % Allow you to exit this stage of the analysis early and move on to the next
        % track group (MA) if you type 99
        
        if clean_track_flag == 99
            disp('Moving on to saving tracks.')
            break
            
        end
        
        
        track_flag_list(i, :) = [3, ref_track, clean_track_flag];
        dlmwrite(fullfile(track_folder_name, manual_pick_file_name), [3, cor_track, clean_track_flag], '-append', 'delimiter', '\t');
        save(fullfile(track_folder_name, 'manually_cleaned_track_list.mat'), 'track_flag_list');
        close all
        
        kk = kk+1;
        
    end
end

%% Convert the text file to a matlab file.
track_flag_list = dlmread(fullfile(track_folder_name, manual_pick_file_name), '\t', 10, 0);
save(fullfile(track_folder_name, 'manually_cleaned_track_list.mat'), 'track_flag_list');
%fclose(fid);

%% ADd variables to the 'associated_tracks_ref_cor.mat' file

save(fullfile(track_folder_name, track_file), 'assocTrackClassificationOrder', 'refOnlyTrackClassificationOrder', 'corOnlyTrackClassificationOrder', 'ii', 'jj', 'kk',  '-append')

%%

disp('Done.')
% [y, Fs] = audioread('notify.wav');
% sound(y, Fs);