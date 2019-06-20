

function [ROI_track_ref, ROI_track_low, ROI_track_high, ...
    track_stat_ref, track_stat_low, track_stat_high] = ...
    throw_tracks_at_timepoint(track_list_ref, position_data_ref, extra_tp)


%% select only the complete track
min_tp = min(position_data_ref(:, 2));
max_tp = max(position_data_ref(:, 2));

complete_track_flag = zeros(numel(track_list_ref),1);
low_track_flag = zeros(numel(track_list_ref),1);
high_track_flag = zeros(numel(track_list_ref),1);

count = 0;
for track = track_list_ref'
    track_id = find(position_data_ref(:,1) == track);
    count = count + 1;
    tp_first = min(position_data_ref(track_id, 2));
    tp_last = max(position_data_ref(track_id, 2));
    
    complete_track_flag(count) = 1;
    
    
    if (tp_first <= extra_tp)
        low_track_flag(count) = 1;
        
    end
    
    if (tp_last > (max_tp - extra_tp))
        high_track_flag(count) = 1;
        
    end
    
end

high_track_list = track_list_ref(high_track_flag==1);
low_track_list = track_list_ref(low_track_flag==1);
track_list_ref = track_list_ref(complete_track_flag==1);

num_track_ref = max(track_list_ref);
num_track_high = max(high_track_list);
num_track_low = max(low_track_list);

%% Convert position_data matrix to a structure data for ref.
ROI_track_ref = struct('tp_first', {}, 'tp_last', {}, 'center_coord',{}, 'coord', {}, 'xcoord', {}, 'ycoord', {});
ROI_track_high = struct('tp_first', {}, 'tp_last', {}, 'center_coord',{}, 'coord', {}, 'xcoord', {}, 'ycoord', {});
ROI_track_low = struct('tp_first', {}, 'tp_last', {}, 'center_coord',{}, 'coord', {}, 'xcoord', {}, 'ycoord', {});


xcoord_ref = NaN*ones(max_tp, num_track_ref);
xcoord_low = NaN*ones(max_tp, num_track_low);
xcoord_high = NaN*ones(max_tp, num_track_high);

ycoord_ref = NaN*ones(max_tp, num_track_ref);
ycoord_low = NaN*ones(max_tp, num_track_low);
ycoord_high = NaN*ones(max_tp, num_track_high);

min_x = min(position_data_ref(:, 3));
max_x = max(position_data_ref(:, 3));
min_y = min(position_data_ref(:, 4));
max_y = max(position_data_ref(:, 4));
track_stat_ref.minx=min_x;
track_stat_ref.miny=min_y;
track_stat_ref.maxx=max_x;
track_stat_ref.maxy=max_y;


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
track_stat_ref.track_list=track_list_ref;


%% For low tracks

for track = low_track_list'
    track_id = find(position_data_ref(:,1) == track);
    
    tp_first = min(position_data_ref(track_id, 2));
    tp_last = max(position_data_ref(track_id, 2));
    tp_prev = max(min_tp, tp_first - extra_tp);
    tp_post = min(max_tp, tp_last + extra_tp);
    
    ROI_track_low(track, 1).tp_first = tp_first;
    ROI_track_low(track, 1).tp_last = tp_last;
    ROI_track_low(track, 1).tp_prev = tp_prev;
    ROI_track_low(track, 1).tp_post = tp_post;
    
    %% XY Coordinates of a track
    track_info_low = NaN*ones((tp_last - tp_first + 1), 2);
    track_info_low(position_data_ref(track_id, 2) - tp_first + 1, :) = position_data_ref(track_id, 3:4);
    
    %% Fill in the gaps by averaging the positions before and after the gap
    track_info_low = fill_gap(track_info_low);
    
    %% save x and y coordinates of tracks.
    xcoord_low(tp_first:tp_last, track) = track_info_low(:,1);
    ycoord_low(tp_first:tp_last, track) = track_info_low(:,2);
    
    ROI_track_low(track, 1).center_coord = track_info_low(:,1:2);
end

track_stat_low.xcoord = xcoord_low;
track_stat_low.ycoord = ycoord_low;
track_stat_low.track_list=low_track_list;



%% For high tracks

for track = high_track_list'
    track_id = find(position_data_ref(:,1) == track);
    
    tp_first = min(position_data_ref(track_id, 2));
    tp_last = max(position_data_ref(track_id, 2));
    tp_prev = max(min_tp, tp_first - extra_tp);
    tp_post = min(max_tp, tp_last + extra_tp);
    
    ROI_track_high(track, 1).tp_first = tp_first;
    ROI_track_high(track, 1).tp_last = tp_last;
    ROI_track_high(track, 1).tp_prev = tp_prev;
    ROI_track_high(track, 1).tp_post = tp_post;
    
    %% XY Coordinates of a track
    track_info_high = NaN*ones((tp_last - tp_first + 1), 2);
    track_info_high(position_data_ref(track_id, 2) - tp_first + 1, :) = position_data_ref(track_id, 3:4);
    
    %% Fill in the gaps by averaging the positions before and after the gap
    track_info_high = fill_gap(track_info_high);
    
    %% save x and y coordinates of tracks.
    xcoord_high(tp_first:tp_last, track) = track_info_high(:,1);
    ycoord_high(tp_first:tp_last, track) = track_info_high(:,2);
    
    ROI_track_high(track, 1).center_coord = track_info_high(:,1:2);
end

track_stat_high.xcoord = xcoord_high;
track_stat_high.ycoord = ycoord_high;
track_stat_high.track_list=high_track_list;

clear position_data_ref


end