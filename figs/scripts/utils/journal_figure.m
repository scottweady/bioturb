
function fig = journal_figure(sz, varargin)
% JOURNAL_FIGURE Create a figure formatted for journal publication

% Get scale factor
try
    sc = varargin{1};
catch
    sc = 1;
end

% Get interpreter
try 
    tex = varargin{2};
catch
    tex = 'latex';
end

% Create figure
fig = figure;

% Set defaults
set(fig,'defaultlinelinewidth', 1.25 * sc)
set(fig,'defaultAxesLineWidth', 1 * sc)
set(fig,'defaultAxesFontSize', 9 * sc)
set(fig,'defaulttextinterpreter',tex); 
set(fig,'defaultAxesTickLabelInterpreter',tex); 
set(fig,'defaultLegendInterpreter',tex);
set(fig,'defaultLegendFontSize', 9 * sc);

% Set figure size
fig.Units = 'inches';
fig.PaperSize = sz * sc;
fig.PaperPosition = [0 0 sz] * sc;
fig.Position = [2 2 sz * sc];
fig.PaperPositionMode = 'auto';
