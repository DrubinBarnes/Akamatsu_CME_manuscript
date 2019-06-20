function h = function_plot_ref_cor_tracks(varargin)
% function_plot_master_slave_tracks
% master and slave track data should be given in the following form:
%   position_data: [Parent (track ID), Time, Position X, Position Y].
% Returns figure handler.
%
%
%      Sun Hae Hong
%      David Drubin lab
%      University of  California, Berkeley
%
%      Copyright 2013
%
%%
if nargin == 2
    position_data_master = varargin{1};
    position_data_slave = varargin{2};
    bg_color = 'w';
    show_start = false;
elseif nargin == 3
    position_data_master = varargin{1};
    position_data_slave = varargin{2};
    bg_color = varargin{3};
    show_start = false;
elseif nargin == 4
    position_data_master = varargin{1};
    position_data_slave = varargin{2};
    bg_color = varargin{3};
    show_start = varargin{4}; % Lable the start of the track with crosses
else
    error('Number of input arguments should be between 2 and 4');
    h = -1;
    return
end


%% make a list of slave tracks and master tracks.
position_data_master = sortrows(position_data_master, [1,2]);
track_list_master = sort(unique(position_data_master(:,1)));

position_data_slave = sortrows(position_data_slave, [1,2]);
track_list_slave = sort(unique(position_data_slave(:,1)));

%% Start a figure and plot 

h = figure;
% set(gcf,'color','w');
whitebg(gcf, bg_color')
for track_num = track_list_master'
    track_id = find(position_data_master(:,1) == track_num);
    plot3(position_data_master(track_id, 3), position_data_master(track_id, 4),position_data_master(track_id, 2), 'Color', [0.8 0 0.8]);
    if show_start
        plot3(position_data_master(track_id(1), 3), position_data_master(track_id(1), 4),position_data_master(track_id(1), 2),'x', 'Color', [0.4 0 0.4]);
    end
    hold on;
end

for track_num = track_list_slave'
    track_id = find(position_data_slave(:,1) == track_num);
    plot3(position_data_slave(track_id, 3), position_data_slave(track_id, 4),position_data_slave(track_id, 2), 'Color', [0 0.8 0.2]);
    if show_start
        plot3(position_data_slave(track_id(1), 3), position_data_slave(track_id(1), 4),position_data_slave(track_id(1), 2), 'x','MarkerSize', 3, 'Color', [0 0.4 0.1]);
    end
    hold on;
end
%
title('All tracks')
min_x = min(min(position_data_slave(:, 3)), min(position_data_master(:, 3)));
max_x = max(max(position_data_slave(:, 3)), max(position_data_master(:, 3)));
min_y = min(min(position_data_slave(:, 4)), min(position_data_master(:, 4)));
max_y = max(max(position_data_slave(:, 4)), max(position_data_master(:, 4)));
min_tp = min(min(position_data_slave(:, 2)), min(position_data_master(:, 2)));
max_tp = max(max(position_data_slave(:, 2)), max(position_data_master(:, 2)));

axis equal
axis([min_x, max_x, min_y, max_y, min_tp, max_tp]);
xlabel('x');
ylabel('y');
zlabel('time');

view(0, 90)

return