
close all
addpath('utils')
fig = journal_figure([6.75 3.25], 2);

%% Configure

% Volume fraction
nu = '0.125';

% Preallocate for storing l2 norm
correlationLength = [];

% Span of box sizes
Lspan = 25 : 25 : 175;

% Initialize L2 norm of velocity field storage
L2norm = [];

% Number of box size simulations
NL = length(Lspan);

% Make colormap
cmap = 1 - linearrgbmap(1 - [60 0 143] / 256, NL + 1);
cmap = cmap(2 : end, :);

%% Plot correlation functions and l2 norms
for nL = 1 : NL

  % Get box size
  L = Lspan(nL);
  
  % Compute l2 norm and plot correlation functions
  l2norm = plotVelocityCorrelation(nu, L, cmap(nL, :));

  % Store correlation length and l2 norm
  L2norm = [L2norm; L, l2norm];

end

%% Load in images and plot
try
  vel50 = resizeVelocity(imread('../snapshots/velocity/nu0.125_L50_1.5.png'));
  vel100 = resizeVelocity(imread('../snapshots/velocity/nu0.125_L100_2.5.png'));
  vel150 = resizeVelocity(imread('../snapshots/velocity/nu0.125_L150_3.75.png'));
catch
  vel50 = [];
  vel100 = [];
  vel150 = [];
end

sp50 = subplot(3, 4, 1);

  imshow(vel50)

sp100 = subplot(3, 4, 5);

  imshow(vel100)

sp150 = subplot(3, 4, [2 6]);

  imshow(vel150)
  velocityColorbarHandle = colorbar;
  colormap(magma)
  velocityColorbarHandle.Units = 'inches';
  velocityColorbarHandle.Ticks = [0 1];
  velocityColorbarHandle.TickLabels = {'$0$', '$|\mathbf{u}|_{\rm max}$'};
  velocityColorbarHandle.TickLabelInterpreter = 'latex';
  velocityColorbarLabel = ylabel(velocityColorbarHandle, '$|\mathbf{u}|$', 'interpreter', 'latex', 'rotation', 0, 'FontSize', 18, 'units', 'inches');

%% Format
sp4 = subplot(3, 4, [3 4 7 8]);
  
  % Plot C = 0 for reference
  plot([0 1], [0 0], 'k--', 'HandleVisibility', 'off')
  r = linspace(0, 0.5, 128);

  xlabel('$r / L$')
  ylabel('$\overline{{\rm Corr}[\mathbf{u}]}$')
  
  xlim([0 0.5])
  ylim([-4e-5 4e-4])
  xticks(0 : 0.1 : 0.5)
  yticks((-1 : 4) * 1e-4)
  colormap(sp4, cmap)
  lengthColorbarHandle = colorbar;
  clim([25 175])
  
  lengthColorbarHandle.TickLabelInterpreter = 'latex';
  lengthColorbarHandle.Ticks = [25 175];
  lengthColorbarHandle.FontSize = 18;
  lengthColorbarLabelHandle = ylabel(lengthColorbarHandle, '$L / \ell$', 'interpreter', 'latex', 'rotation', 0);

sp5 = subplot(3, 4, [9 10]);

  legend('location', 'northwest')
  xlabel('$t$')
  ylabel('$||\mathbf{u}||_2$')
  
  ax = gca;
  xlim([0 300])
  ylim([0 0.025])

sp6 = subplot(3, 4, [11 12]);
  
  plot(L2norm(:, 1), L2norm(:, 2) ./ L2norm(:, 1), 'o-', 'Color',  [0.1375 0.0685 0.3191])
  ylim([0 0.025])
  
  hold on
  
  yspan = ylim;
  fill([0 75 75 0 0], [yspan(1) yspan(1) yspan(2) yspan(2) yspan(1)], 'k:', 'FaceAlpha', 0.05, 'LineWidth', 2)
  
  xlabel('$L / \ell$')
  ylabel('$\overline{||\mathbf{u}||_2}$')

%% Position subplots
sp50.Units = 'inches';
sp100.Units = 'inches';
sp150.Units = 'inches';

sp50.Position([3 4]) = 1.1;

sp100.Position(1) = 0.7;
sp100.Position([3 4]) = 2 * sp50.Position([3 4]);
sp100.Position(2) = fig.PaperSize(2) - sp100.Position(4) - sp50.Position(4) - 0.2;

sp150.Position(1) = sum(sp100.Position([1 3]))-0.2;
sp150.Position(2) = sp100.Position(2);
sp150.Position([3 4]) = 3 * sp50.Position([3 4]);

sp50.Position(2) = sum(sp100.Position([2 4]));
sp50.Position(1) = sum(sp100.Position([1 3])) - sp50.Position(3) - 0.25 * sp100.Position(3);

sp4.Units = 'inches';
sp5.Units = 'inches';
sp6.Units = 'inches';

sp4.FontSize = 18;
sp5.FontSize = 18;
sp6.FontSize = 18;

sp4.Position(1) = sum(sp150.Position([1 3])) + 1.75;
sp4.Position(2) = sp150.Position(2) + 0.4;
sp4.Position(4) = 0.8 * sp150.Position(4);
sp4.Position(3) = 0.35 * fig.PaperSize(1);

sp5.Position(1) = sp50.Position(1) - 0.3;
sp5.Position(2) = 0.7;
sp5.Position(3) = 0.6 * (sum(sp50.Position([1 3])) + sum(sp150.Position([1 3])));
sp5.Position(4) = 0.6 * (fig.PaperSize(2) - sp150.Position(4));

sp6.Position([1 3]) = sp4.Position([1 3]);
sp6.Position([2 4]) = sp5.Position([2 4]);

lengthColorbarHandle.Units = 'inches';
lengthColorbarHandle.Position(1) = sum(sp4.Position([1 3])) + 0.1;

lengthColorbarLabelHandle.Units = 'inches';
lengthColorbarLabelHandle.Position(1) = 0.6;
lengthColorbarLabelHandle.Position(2) = 0.5 * lengthColorbarHandle.Position(4) + 0.15;

scaleBar = annotation('rectangle', 'units', 'inches', 'FaceColor', 0.2 * [1 1 1], 'EdgeColor', 'none');
scaleBar.Position(3) = 0.1 * sp50.Position(3);
scaleBar.Position(4) = 0.5 * (88.5564 / 80) * sp50.Position(4);
scaleBar.Position(1) = sp100.Position(1) - 2 * scaleBar.Position(3);
scaleBar.Position(2) = sp100.Position(2) + 0.2;
scaleBarLabel = annotation('textbox', 'units', 'inches', 'fontsize', ax.FontSize, 'string', '$50\ell$', 'interpreter', 'latex', 'edgecolor', 'none');
scaleBarLabel.Position(1) = scaleBar.Position(1) - 0.5;
scaleBarLabel.Position(2) = scaleBar.Position(2) - 0.2;

velocityColorbarHandle.Position(4) = 0.6 * sp150.Position(4);
velocityColorbarHandle.Position(2) = sp150.Position(2) + 0.2 * sp150.Position(4);
velocityColorbarLabel.Position(1) = 0.5;

A = subplotLabel('(a)', gca);
B = subplotLabel('(b)', sp4, 'northwest', [-1 -0.4]);
C = subplotLabel('(c)', sp5, 'northwest', [-1 -0.4]);
D = subplotLabel('(d)', sp6, 'northwest', [-1 -0.4]);

A.Position(1) = 0;
A.Position(2) = B.Position(2);

%% Plotting function
function l2norm = plotVelocityCorrelation(nu, L, col)

  % Get file name
  source = sprintf('../../data/processed/correlation_functions/nu%s_L%d', nu, L);

  r = load(strcat(source, '/r.dat'));
  r = r(1, :);

  C = load(strcat(source, '/u.dat'));
  Nt = size(C, 1);
    
  % Preallocate for storing l2 norm
  l2Evolution = sqrt(C(:, 1));

  % Store time
  t = 0 : (Nt - 1);
  
  % Averaging interval
  Nspan = min(200, Nt) : Nt;

  Cbar = mean(C(Nspan, :), 1);

  % Subplot 1 : Evolution of the L2 norm
  subplot(3, 4, [9 10])
  plot(t, l2Evolution / L, '-', 'Color', col, 'DisplayName', sprintf('$L = %d$', L), 'HandleVisibility', 'off')
  hold on

  % Subplot 2 : Correlation functions
  subplot(3, 4, [3 4 7 8])
  plot(r(2 : end) / L, Cbar(2 : end) / L^2, '-', 'Color', col, 'DisplayName', 'simulation')
  hold on

  l2norm = mean(l2Evolution(Nspan));

end
