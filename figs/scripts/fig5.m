
close all
addpath('utils')
fig = journal_figure([6.75 3.325], 2);

%% Configure

% Volume fraction
nu = '1.0';

% Preallocate for storing l2 norm
l2norm = [];

% Span of box sizes
Lspan = 25 : 25 : 150;

% Number of box size simulations
NL = length(Lspan);

% Colormap
colors = cmocean('phase', 5);
colors = colors(2 : end - 1, :);

%% Plot correlation functions
for i = 2 : 4

  col = colors(i - 1, :);
  cmap = 1 - linearrgbmap(1 - col, NL + 1);
  cmap = cmap(2 : end, :);
  
  % Subplot 1, 2 : evolution of l2 norm and correlation function
  for nL = 1 : NL
  
    % Get box size
    L = Lspan(nL);
    
    % Plot correlation functions
    plotCorrelationFunction(nu, L, cmap(nL, :), i, 50, 1);
  
  end

end

%% Format
subplot_list = cell(3, 3);

for i = 2 : 4

  subplot_list{i-1} = subplot(3, 4, 7 + i);
  
  plot([0 20], [0 0], 'k--', 'HandleVisibility', 'off')
  xlabel('$r / \ell$')

  switch i

    case 2
        subplot_list{i - 1, 3} = ylabel('Corr[$c''$]');
        ylim([-0.005 0.02])
        yticks(0 : 0.01 : 0.02)
        yticklabels({'0', '0.01', '0.02'})
        
    case 3
        subplot_list{i - 1, 3} = ylabel('Corr[$\mathbf{n}$]');
        ylim([-0.01 0.04])
        yticks(0 : 0.02 : 0.04)
        yticklabels({'0', '0.02', '0.04'})
        
    case 4
        subplot_list{i - 1, 3} = ylabel('Corr[$\mathbf{Q}$]');
        ylim([-0.01 0.09])
        yticks(0 : 0.04 : 0.08)
        yticklabels({'0','0.04','0.08'})

  end

  xlim([0 15])
  xticks(0 : 5 : 15)
  
  col = colors(i - 1, :);
  
  cmap = 1 - linearrgbmap(1 - col, NL + 1);
  cmap = cmap(2 : end, :);

  colormap(subplot_list{i - 1}, cmap)
  subplot_list{i - 1, 2} = colorbar;
  clim([min(Lspan) max(Lspan)])
  
  subplot_list{i - 1, 2}.TickLabelInterpreter = 'latex';
  subplot_list{i - 1, 2}.FontSize = 18;
  subplot_list{i - 1, 2}.Ticks = [min(Lspan) max(Lspan)];

end

formatCorrelationFunctions(fig, subplot_list)
plotSnapshots(nu)