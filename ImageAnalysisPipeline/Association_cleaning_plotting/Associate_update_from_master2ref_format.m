% Goes through the given folder. Converts associated_tracks_master_slave.mat to associated_tracks_ref_cor.mat
%%
clear
close all
clc

%% Select folder for all of the data
load('previous_folder.mat')
data_folder = uigetdir(previous_folder, 'Pick folder of all analysis.');
if isdir(data_folder)
    previous_folder = fullfile(data_folder, '');
    save('previous_folder', 'previous_folder')
else
    display('You did not pick a folder.')
    return
end
%% Parameter to use
% align_grid is the size of the grid to be used for alignment. An image
% with a pixel size of align_grid is made.
% number of extra time points before and after the track.
extra_tp = 5;
threshold = 1/350^2; % distance btw two spots are less than 350nm.
montage_size = 3;
peak_region_mask = fspecial('disk', montage_size) > 0.01;
bg_size = 5;
bg_region_mask = (fspecial('disk', bg_size) > 0.001);
peak_index_within_bg = (bg_size-montage_size+1):(bg_size+montage_size+1);
bg_region_mask(peak_index_within_bg, peak_index_within_bg)=bg_region_mask(peak_index_within_bg, peak_index_within_bg).*(1-peak_region_mask);

%% Go through folders to pool data.
folder_ID1 = 0;
list1 = dir(data_folder);
for foldernum1 =3:numel(list1)
    if list1(foldernum1).isdir
        list2 = dir(fullfile(data_folder, list1(foldernum1).name));
        folder_ID1 = folder_ID1 + 1;
        folder_ID2 = 0;
        for foldernum2 = 3:numel(list2)
            if list2(foldernum2).isdir
                
                % load tracking data and manual tag data.
                track_folder_name = fullfile(data_folder, list1(foldernum1).name, list2(foldernum2).name);
                if exist(fullfile(track_folder_name, 'associated_tracks_master_slave.mat'), 'file');
                    load(fullfile(track_folder_name, 'associated_tracks_master_slave.mat'));
                    pixel_size = analysis_info.pixel_size; % in nm
                else
                    display([fullfile(track_folder_name, 'associated_tracks_master_slave.mat'), ' is not found.']);
                    continue
                end
                display(['Creating: ',fullfile(track_folder_name, 'associated_tracks_master_slave.mat')]);
                
                % Retrieve the file names
                
                [ref_file, ref_folder_name]= divide_folder_file_name(analysis_info.master_csv_name);
                [cor_file, cor_folder_name]= divide_folder_file_name(analysis_info.slave_csv_name);
                [ref_img_file, ref_img_folder_name] = divide_folder_file_name(analysis_info.master_image_name);
                [cor_img_file, cor_img_folder_name] = divide_folder_file_name(analysis_info.slave_image_name);
                
                %% if thes file cannot be found in the designated folder, search the file with the same name in the corresponding folder
                if ~exist(fullfile(ref_img_folder_name, ref_img_file), 'file') || ~exist(fullfile(cor_img_folder_name, cor_img_file), 'file')
                    if exist(fullfile(track_folder_name, ref_img_file), 'file') &&  exist(fullfile(track_folder_name, cor_img_file), 'file')
                        ref_img_folder_name = track_folder_name;
                        cor_img_folder_name = track_folder_name;
                    else
                        display([track_folder_name, ' cannot be analyzed due to missing image file(s).']);
                    end
                end
                %% if thes file cannot be found in the designated folder, search the file with the same name in the corresponding folder
                if ~exist(fullfile(ref_folder_name, ref_file), 'file') || ~exist(fullfile(cor_folder_name, cor_file), 'file')
                    if exist(fullfile(track_folder_name, ref_file), 'file') &&  exist(fullfile(track_folder_name, cor_file), 'file')
                        ref_folder_name = track_folder_name;
                        cor_folder_name = track_folder_name;
                    else
                        display([track_folder_name, ' cannot be analyzed due to missing file(s).']);
                    end
                end
                analysis_info = rmfield(analysis_info, 'master_csv_name');
                analysis_info = rmfield(analysis_info, 'slave_csv_name');
                analysis_info = rmfield(analysis_info, 'master_image_name');
                analysis_info = rmfield(analysis_info, 'slave_image_name');
                
                % Update the channel info
                analysis_info.ref_channel = analysis_info.master_channel;
                analysis_info.cor_channel = analysis_info.slave_channel;
                
                analysis_info = rmfield(analysis_info, 'master_channel');
                analysis_info = rmfield(analysis_info, 'slave_channel');
                
                %% read in the csv file for ref.
                % the columns are for Position X, Position Y, Unit, Category, Collection, Time, Parent (Track num), ID
                text_scan_format = '%f%f%s%s%s%f%f%f';
                csv_file_ref = fullfile(ref_folder_name, ref_file);
                fid = fopen(csv_file_ref);
                if strcmp(ref_file(end-12:end), '_Position.csv')
                    ref_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',4);
                    if numel(ref_data{8}) == 0
                        text_scan_format = '%f%f%f%s%s%s%f%f%f';
                        ref_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',4);
                    end
                elseif strcmp(ref_file(end-12:end), '.Position.csv')
                    ref_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',2);
                    if numel(ref_data{8}) == 0
                        text_scan_format = '%f%f%f%s%s%s%f%f%f';
                        ref_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',4);
                    end
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
                if strcmp(cor_file(end-12:end), '_Position.csv')
                    cor_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',4);
                elseif strcmp(cor_file(end-12:end), '.Position.csv')
                    cor_data = textscan(fid,text_scan_format,'delimiter',',','headerlines',2);
                end
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
                %                 h = function_plot_ref_cor_tracks(position_data_ref, position_data_cor);
                %                 figure(h); title('Original track positions.')
                
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
                shift_y_nm = analysis_info.shift_y_nm;
                shift_x_nm = analysis_info.shift_x_nm;
                shift_y = analysis_info.shift_y;
                shift_x = analysis_info.shift_x;
                
                %% select only the complete track
                min_tp = min(position_data_ref(:, 2));
                max_tp = max(position_data_ref(:, 2));
                
                complete_track_flag = ones(numel(track_list_ref),1);
                count = 0;
                for track = track_list_ref'
                    track_id = find(position_data_ref(:,1) == track);
                    count = count + 1;
                    tp_first = min(position_data_ref(track_id, 2));
                    tp_last = max(position_data_ref(track_id, 2));
                    
                    if or(tp_first <= extra_tp, tp_last > (max_tp - extra_tp))
                        complete_track_flag(count) = false;
                        continue
                    end
                end
                track_list_ref = track_list_ref(complete_track_flag==1);
                num_track_ref = max(track_list_ref);
                
                %% Convert position_data matrix to a structure data for ref.
                ROI_track_ref = struct('tp_first', {}, 'tp_last', {}, 'center_coord',{}, 'coord', {}, 'xcoord', {}, 'ycoord', {});
                xcoord_ref = NaN*ones(tp_max, num_track_ref);
                ycoord_ref = NaN*ones(tp_max, num_track_ref);
                
                min_x = min(position_data_ref(:, 3));
                max_x = max(position_data_ref(:, 3));
                min_y = min(position_data_ref(:, 4));
                max_y = max(position_data_ref(:, 4));
                
                
                for track = track_list_ref'
                    track_id = find(position_data_ref(:,1) == track);
                    
                    tp_first = min(position_data_ref(track_id, 2));
                    tp_last = max(position_data_ref(track_id, 2));
                    tp_prev = max(min_tp, tp_first - extra_tp);
                    tp_post = min(max_tp, tp_last + extra_tp);
                    
                    ROI_track_ref(track, 1).tp_first = tp_first;
                    ROI_track_ref(track, 1).tp_last = tp_last;
                    ROI_track_ref(track, 1).tp_prev = tp_prev;
                    ROI_track_ref(track, 1).tp_post = tp_post;
                    
                    %% XY Coordinates of a track
                    track_info_ref = NaN*ones((tp_last - tp_first + 1), 2);
                    track_info_ref(position_data_ref(track_id, 2) - tp_first + 1, :) = position_data_ref(track_id, 3:4);
                    
                    %% Fill in the gaps by averaging the positions before and after the gap
                    track_info_ref = fill_gap(track_info_ref);
                    
                    %% save x and y coordinates of tracks.
                    xcoord_ref(tp_first:tp_last, track) = track_info_ref(:,1);
                    ycoord_ref(tp_first:tp_last, track) = track_info_ref(:,2);
                    
                    ROI_track_ref(track, 1).center_coord = track_info_ref(:,1:2);
                end
                
                track_stat_ref.xcoord = xcoord_ref;
                track_stat_ref.ycoord = ycoord_ref;
                
                clear position_data_ref
                
                %% select only the complete track
                min_tp = min(position_data_cor(:, 2));
                max_tp = max(position_data_cor(:, 2));
                
                complete_track_flag = ones(numel(track_list_cor),1);
                count = 0;
                for track = track_list_cor'
                    track_id = find(position_data_cor(:,1) == track);
                    count = count + 1;
                    tp_first = min(position_data_cor(track_id, 2));
                    tp_last = max(position_data_cor(track_id, 2));
                    
                    if or(tp_first <= extra_tp, tp_last > (max_tp - extra_tp))
                        complete_track_flag(count) = false;
                        continue
                    end
                end
                track_list_cor = track_list_cor(complete_track_flag==1);
                num_track_cor = max(track_list_cor);
                
                %% Convert position_data matrix to a structure data for cor.
                ROI_track_cor= struct('tp_first', {}, 'tp_last', {}, 'xcoord', {}, 'ycoord', {});
                xcoord_cor= NaN*ones(tp_max, num_track_cor);
                ycoord_cor= NaN*ones(tp_max, num_track_cor);
                
                track_stat_cor.extra_tp = extra_tp;
                
                min_x = min(position_data_cor(:, 3));
                max_x = max(position_data_cor(:, 3));
                min_y = min(position_data_cor(:, 4));
                max_y = max(position_data_cor(:, 4));
                
                count = 0;
                for track = track_list_cor'
                    track_id = find(position_data_cor(:,1) == track);
                    count = count + 1;
                    
                    tp_first = min(position_data_cor(track_id, 2));
                    tp_last = max(position_data_cor(track_id, 2));
                    
                    tp_prev = max(min_tp, tp_first - extra_tp);
                    tp_post = min(max_tp, tp_last + extra_tp);
                    
                    ROI_track_cor(track, 1).tp_first = tp_first;
                    ROI_track_cor(track, 1).tp_last = tp_last;
                    ROI_track_cor(track, 1).tp_prev = tp_prev;
                    ROI_track_cor(track, 1).tp_post = tp_post;
                    
                    %% XY Coordinates of a track
                    track_info_cor= NaN*ones((tp_last - tp_first + 1), 2);
                    track_info_cor(position_data_cor(track_id, 2) - tp_first + 1, :) = position_data_cor(track_id, 3:4);
                    
                    %% Fill in the gaps by averaging the positions before and after the gap
                    track_info_cor= fill_gap(track_info_cor);
                    
                    %% save x and y coordinates of tracks.
                    xcoord_cor(tp_first:tp_last, track) = track_info_cor(:,1);
                    ycoord_cor(tp_first:tp_last, track) = track_info_cor(:,2);
                    
                    ROI_track_cor(track, 1).center_coord = track_info_cor(:,1:2);
                end
                track_stat_cor.xcoord = xcoord_cor;
                track_stat_cor.ycoord = ycoord_cor;
                
                clear position_data_cor
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
                
                ref_matrix_size2 = size(xcoord_ref, 2);
                cor_matrix_size2 = size(xcoord_cor, 2);
                
                % ref_track_info = [ref track id, number of reference tracks, number of corresponding trcaks].
                ref_track_info = zeros*ones(numel(track_list_ref), 3);
                count = 0;
                for ref_track = track_list_ref'
                    if sum(ref_track == ref_tracks_taken)
                        continue
                    end
                    ref_primary_track = [ref_primary_track; ref_track];
                    ref_tracks_taken = [ref_tracks_taken; ref_track];
                    
                    dist_score_mat = (xcoord_cor- repmat(xcoord_ref(:, ref_track), [1, cor_matrix_size2]) - shift_x_nm).^2 +...
                        (ycoord_cor- repmat(ycoord_ref(:, ref_track), [1, cor_matrix_size2]) - shift_y_nm).^2;
                    dist_score_mat = 1./dist_score_mat;
                    dist_score = nansum(dist_score_mat,1);
                    overlapping_time = max(10, sum(~isnan(dist_score_mat),1));
                    candidates = find(dist_score > threshold*overlapping_time/2);
                    
                    ROI_x = [round(min(xcoord_ref(:,ref_track))/pixel_size - 5), round(max(xcoord_ref(:,ref_track))/pixel_size + 5)];
                    ROI_y = [round(min(ycoord_ref(:,ref_track))/pixel_size - 5), round(max(ycoord_ref(:,ref_track))/pixel_size + 5)];
                    
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
                        dist_score_mat = (xcoord_ref +shift_y_nm - repmat(xcoord_cor(:, cor_track), [1, ref_matrix_size2])).^2 +...
                            (ycoord_ref+shift_x_nm - repmat(ycoord_cor(:, cor_track), [1, ref_matrix_size2])).^2;
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
                track_stat_cor.cor_tracks_unassociated = setdiff(track_list_cor, cor_tracks_taken);
                
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
                                max(ROI_track_ref(ref_track).ref_coord_extrapol(:,1) + shift_x_nm) <= (max_x - bg_size*pixel_size) &&...
                                min(ROI_track_ref(ref_track).ref_coord_extrapol(:,2) + shift_y_nm) >= (1+bg_size)*pixel_size &&...
                                max(ROI_track_ref(ref_track).ref_coord_extrapol(:,2) + shift_y_nm) <= (max_y - bg_size*pixel_size)
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
                for cor_track = track_stat_cor.cor_tracks_unassociated
                    
                    ROI_x = [round(min(xcoord_cor(:,cor_track))/pixel_size - 5), round(max(xcoord_cor(:,cor_track))/pixel_size + 5)];
                    ROI_y = [round(min(ycoord_cor(:,cor_track))/pixel_size - 5), round(max(ycoord_cor(:,cor_track))/pixel_size + 5)];
                    
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
                    % if the coordinate is not at the boundary, assign the
                    % known coordinate to the unknown coordinate.
                    if min(ROI_track_cor(cor_track).cor_coord_extrapol(:,1) - shift_x_nm) >= (1+bg_size)*pixel_size &&...
                            max(ROI_track_cor(cor_track).cor_coord_extrapol(:,1) - shift_x_nm) <= (max_x - bg_size*pixel_size) &&...
                            min(ROI_track_cor(cor_track).cor_coord_extrapol(:,2) - shift_y_nm) >= (1+bg_size)*pixel_size &&...
                            max(ROI_track_cor(cor_track).cor_coord_extrapol(:,2) - shift_y_nm) <= (max_y - bg_size*pixel_size)
                        ROI_track_cor(cor_track).ref_coord_extrapol(:,1) = ROI_track_cor(cor_track).cor_coord_extrapol(:,1) - shift_x_nm;
                        ROI_track_cor(cor_track).ref_coord_extrapol(:,2) = ROI_track_cor(cor_track).cor_coord_extrapol(:,2) - shift_y_nm;
                    end
                    
                end
                
                
                %% save montages of the corresponding tracks that are not associated with ref.
                unassociated_tracks = struct('ref_movie', {},'cor_movie', {});
                for cor_track = track_stat_cor.cor_tracks_unassociated
                    
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
                analysis_info.extrp_tp = extra_tp;
                
                %%
                save(fullfile(cor_folder_name, 'associated_tracks_ref_cor.mat'),...
                    'position_data_ref_original', 'position_data_cor_original', ...
                    'associated_tracks', 'unassociated_tracks', 'track_stat_ref', 'track_stat_cor', 'ROI_track_ref', 'ROI_track_cor', 'analysis_info');
                
                display(['The association results are saved as -', fullfile(cor_folder_name, 'associated_tracks_ref_cor.mat')]);

            end
        end
    end
end
%%
[y, Fs] = audioread('notify.wav');
sound(y, Fs);