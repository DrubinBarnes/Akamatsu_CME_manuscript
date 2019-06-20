     %% Make a GUI to pick which stats to plot so you don't always get a million little pictures

function [] = plotStats_dif_timesteps()
%% different ways that the data can be plotted, change S.ops and corrFuncs to add additional ones
S.ops =  {"Associated corr stats", "Associated ref stats",...
    "Unassociated corr", "Unassociated ref", ...
    "Aligned corr disappear", "Aligned corr start", ...
    "Aligned ref end", "Aligned ref start", ...
    "Alignted ref max", "Aligned corr max",...
    "All profiles", "Lifetimes histogram", ...
    "Associated lifetimes histogram", "Half Dropoff",...
    "Compare Stats"
};

S.funcs = {@plotAssociatedCorr, @plotAssociatedRef,...
    @plotUnassociatedCorr, @plotUnassociatedRef, @plotAlignCorrDisappear, @plotAlignCorrStart, ...
    @plotAlignRefEnd, @plotAlignRefStart, @plotAlignRefMax, @plotAlignCorrMax, @plotProfiles, ...
    @lifetimeHistogram, @AssociatedLifetimeHistogram, @plotHalfDropoff, @compare_stats};
S.plotAll = false;

%positions are determined [left bottom width height]
S.CNT = zeros(1, length(S.ops));
S.fh = figure('units','pixels',...
              'position',[500 200 500 650],...
              'menubar','none',...
              'name','Select Graphs to Plot',...
              'numbertitle','off',...
              'resize','on');

S.done = uicontrol('Style', 'pushbutton', 'String', 'done', ...
    'Position', [60 20 150 20], 'Callback', {@completedSelection, S});

for i = 1:length(S.ops)
    S.bt(i) = uicontrol('Style', 'checkbox', 'String', S.ops{i},...
        'Position', [20 20+(30*i) 300 20]);
end

maxBox = 20 + (30 * numel(S.ops));

% S.ed = uicontrol('style','edit', 'units','pix','position',[60 maxBox + 30 300 20],...
%                  'string',{'Program Name'}, 'horizontalalign','center');
             
uicontrol('style', 'text', 'string', 'What is the time interval between time points?',...
    'position', [60 maxBox + 90 300 20]);

S.cellSelect = uicontrol('style', 'checkbox', 'string', 'Plot all cells?',...
    'position', [250 maxBox + 120 150 20]);

S.saveData = uicontrol('style', 'checkbox', 'string', 'Write to csv?',...
    'position', [375 maxBox + 120 150 20]);
S.save = false;

S.time = uicontrol('Style', 'edit', 'string', '0',...
    'Position', [60 maxBox + 60 150 20]);

uicontrol('style', 'pushbutton', 'string', 'Pick Analysis Folder',...
    'position', [60 maxBox + 120 150 20], 'callback', {@pickFolder, S});

% set all the callbacks once S finished
assign(S);
end

function [] = getSaveData(~, ~, S)
S.save = ~S.save;
assign(S);
end

function [] = setTimePoint(~, ~, S)
S.freq = get(S.time, 'string');
assign(S);
end

function [] = plotAllCells(~, ~, S)
S.plotAll = ~S.plotAll;
assign(S);
end

function [] = pickFolder(~, ~, S)
folder = uigetdir('Pick folder of all analysis.');
S.folder = folder;
assign(S);
end

function [] = addGraph(varargin)
% Callback for popupmenu.
index = varargin{4};
S = varargin{3};
S.CNT(index) = ~S.CNT(index);
% update the S in all of the callbacks - is there a better way to do this
assign(S);
end

function completedSelection(varargin)
S = varargin{3};
delete(S.fh);
disp(S);
plotFuncs_dif_timesteps(S);
end

function assign(S)
for j = 1:length(S.bt)
    set(S.bt(j), 'callback', {@addGraph, S, j});
end
set(S.done, 'callback', {@completedSelection, S});
set(S.time,'callback',{@setTimePoint, S});
set(S.cellSelect, 'callback', {@plotAllCells, S});
set(S.saveData, 'callback', {@getSaveData, S});
end


