%% This function takes position data of the following format.
%% position_data: [Parent, Time, Position X, Position Y].
%% Tracks that touch the edges or within "boundary_area" distance from the edges are removed.
%% "position_data_cleaned" has only the tracks that are not touching the boundary.


function position_data_cleaned = throw_tracks_at_edges(position_data, boundary_area)

%% Parameters for debugging
show_process = false;
%% Make a list of tracks
position_data = sortrows(position_data, [1,2]);

track_list = sort(unique(position_data(:,1)));

%% Determine the edges. Label tracks at the edge as bad.
min_x1 = min(position_data(:, 3)) + boundary_area;
max_x1 = max(position_data(:, 3)) - boundary_area;
min_x2 = min(position_data(:, 4)) + boundary_area;
max_x2 = max(position_data(:, 4)) - boundary_area;

bad_tracks = (position_data(:, 3) <= min_x1) | (position_data(:, 3) >= max_x1) |...
    (position_data(:, 4) <= min_x2)| (position_data(:, 4) >= max_x2);
bad_track_list = unique(position_data(bad_tracks));
good_track_list = setdiff(track_list, bad_track_list);
if show_process
    %% Plot all of the tracks
    tic
    figure;
    for track_num = track_list'
        track_id = find(position_data(:,1) == track_num);
        plot(position_data(track_id, 3), position_data(track_id, 4), 'b'); hold on
    end
    xlabel('x')
    ylabel('y')
    axis equal
    title('All tracks');
    toc
    %%
    figure;
    tic
    for track_num = good_track_list'
        track_id = find(position_data(:,1) == track_num);
        plot(position_data(track_id, 3), position_data(track_id, 4), 'b'); hold on
    end
    toc
    
    tic
    for track_num = bad_track_list'
        track_id = find(position_data(:,1) == track_num);
        plot(position_data(track_id, 3), position_data(track_id, 4), 'k'); hold on
    end
    toc
    xlabel('x')
    ylabel('y')
    axis equal
    title('Good tracks in blue, tracks at edges in red');
end
%% Make a matrix only of tracks that are not touching the edges.
good_track_mask = false*ones(size(position_data, 1),1);
for track_num = good_track_list'
    track_id = position_data(:,1) == track_num;
    good_track_mask(track_id,1) = true;
end
position_data_cleaned = position_data(find(good_track_mask), :);
% if show_process
%     figure,
%     plot(position_data(:, 3),position_data(:,4), '.r');
%     plot(position_data_cleaned(:, 3),position_data_cleaned(:,4), '.k');
%     title('Red dots are discarded because they belong to tracks at the edge.')
% end
return
