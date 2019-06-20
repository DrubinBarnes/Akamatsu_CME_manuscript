function h = function_errorbar(varargin)
% This function draws stylish errorbar.
% 
%
%      Sun Hae Hong
%      David Drubin lab
%      University of  California, Berkeley
%
%      Copyright 2014
%
%%
if nargin == 3
    y_data = varargin{1};
    if numel(y_data) ~= max(size(y_data))
        error('The input for function_errorbar needs to be an array.');
        return
    end
    error_low = varargin{2};
    error_high = varargin{3};
    x_data = 1:numel(y_data);
    h = figure;
    plot_color = [0 0 1];
elseif nargin == 4
    x_data = varargin{1};
    y_data = varargin{2};
    if numel(y_data) ~= max(size(y_data))
        error('The input for function_errorbar needs to be an array.');
        return
    end
    error_low = varargin{3};
    error_high = varargin{4};
    h = figure;
    plot_color = [0 0 1];
elseif nargin == 5
    x_data = varargin{1};
    y_data = varargin{2};
    if numel(y_data) ~= max(size(y_data))
        error('The input for function_errorbar needs to be an array.');
        return
    end
    error_low = varargin{3};
    error_high = varargin{4};
    h = varargin{5};
    plot_color = [0 0 1];
elseif nargin == 6
    x_data = varargin{1};
    y_data = varargin{2};
    if numel(y_data) ~= max(size(y_data))
        error('The input for function_errorbar needs to be an array.');
        return
    end
    error_low = varargin{3};
    error_high = varargin{4};
    h = varargin{5};
    plot_color = varargin{6};
end

%%
% Change the shape of the array
if size(x_data,1) > size(x_data,2)
    x_data = x_data';
end
if size(y_data,1) > size(y_data,2)
    y_data = y_data';
end
if size(error_low,1) > size(error_low,2)
    error_low = error_low';
end
if size(error_high,1) > size(error_high,2)
    error_high = error_high';
end
if size(plot_color,1) > size(plot_color,2)
    plot_color = plot_color';
end
figure(h);
line([x_data; x_data], [error_low; error_high], 'Color', min([plot_color+.2; 1,1,1],[], 1));hold on;
plot(x_data, y_data, 'Color', plot_color,'Linewidth', 1); 

return
