
close all

% Initialize figure
fig = journal_figure([6.75 3.25], 2);

% Volume fraction
nu = '1.0';

% Field to plot
i = 5; %2 - concentration, 3 - polarity, 4 - nematic, 5 - velocity

% Preallocate for storing l2 norm
correlationLength = [];

% Span of box sizes
Lspan = 25 : 25 : 150;

% Initialize L2 norm of velocity field storage
L2norm = [];

% Number of box size simulations
NL = length(Lspan);

% Make colormap
cmap = 1 - linearrgbmap(1 - [60 0 143] / 256, NL + 1);
cmap = cmap(2 : end, :);

% Subplot 1, 2 : evolution of l2 norm and correlation function
for nL = 1 : NL

  % Get box size
  L = Lspan(nL);
  
  % Compute l2 norm and plot correlation functions
  [xi, l2norm] = render(nu, L, cmap(nL, :));

  % Store correlation length and l2 norm
  correlationLength = [correlationLength; L, xi];
  L2norm = [L2norm; L, l2norm];

end

% Load in images
vel50 = recolorImage(imread('../snapshots/velocity/nu1.0_L50_6.png'));
vel100 = recolorImage(imread('../snapshots/velocity/nu1.0_L100_7.5.png'));
vel150 = recolorImage(imread('../snapshots/velocity/nu1.0_L150_8.png'));

% Plot particle simulations
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

% Subplot 2 : Correlation function; create labels
sp4 = subplot(3, 4, [3 4 7 8]);
  
  % Plot C = 0 for reference
  plot([0 1], [0 0], 'k--', 'HandleVisibility', 'off')
  r = linspace(0, 0.5, 128);
  color = 'k';

  xlabel('$r / L$')
  ylabel('$\overline{{\rm Corr}[\mathbf{u}]}$')
  
  
  xlim([0 0.5])
  ylim([-1e-3 2e-2])
  xticks(0 : 0.1 : 0.5)
  colormap(sp4, cmap)
  lengthColorbarHandle = colorbar;
  clim([25 150])
  
  lengthColorbarHandle.TickLabelInterpreter = 'latex';
  lengthColorbarHandle.Ticks = [25 150];
  lengthColorbarHandle.FontSize = 18;
  lengthColorbarLabelHandle = ylabel(lengthColorbarHandle, '$L / \ell$', 'interpreter', 'latex', 'rotation', 0);

sp5 = subplot(3, 4, [9 10]);

  legend('location', 'northwest')
  xlabel('$t$')
  ylabel('$||\mathbf{u}||_2$')
  
  ax = gca;
  xlim([0 200])
  ylim([0 0.2])

sp6 = subplot(3, 4, [11 12]);
  
  plot(L2norm(:, 1), L2norm(:, 2) ./ L2norm(:, 1), 'o-', 'Color',  [0.1375 0.0685 0.3191])
  ylim([0 0.2])
  
  hold on
  
  Lgrid = linspace(0, 150, 128);
  plot(Lgrid, L2norm(end, 2) / L2norm(end, 1) * (Lgrid(end)./Lgrid), '--', 'Color', [0.1375 0.0685 0.3191])
  
  xlabel('$L / \ell$')
  ylabel('$\overline{||\mathbf{u}||_2}$')

sp50.Units = 'inches';
sp100.Units = 'inches';
sp150.Units = 'inches';

sp50.Position([3 4]) = 1.1;

sp100.Position(1) = 0.7;
sp100.Position([3 4]) = 2 * sp50.Position([3 4]);
sp100.Position(2) = fig.PaperSize(2) - sp100.Position(4) - sp50.Position(4) - 0.2;

sp150.Position(1) = sum(sp100.Position([1 3])) - 0.2;
sp150.Position(2) = sp100.Position(2);
sp150.Position([3 4]) = 3 * sp50.Position([3 4]);

sp50.Position(2) = sum(sp100.Position([2 4]));
sp50.Position(1) = sum(sp100.Position([1 3])) - sp50.Position(3) - 0.25 * sp100.Position(3);

sp4.Units = 'inches';
sp5.Units = 'inches';
sp6.Units = 'inches';

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

sp4.FontSize = 18;
sp5.FontSize = 18;
sp6.FontSize = 18;

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
scaleBarLabel = annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', '$50\ell$', 'interpreter', 'latex', 'edgecolor', 'none');
scaleBarLabel.Position(1) = scaleBar.Position(1) - 0.5;
scaleBarLabel.Position(2) = scaleBar.Position(2) - 0.2;

velocityColorbarHandle.Position(4) = 0.6 * sp150.Position(4);
velocityColorbarHandle.Position(2) = sp150.Position(2) + 0.2 * sp150.Position(4);
velocityColorbarLabel.Position(1) = 0.5;

label = @(s) annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', s, 'interpreter', 'latex', 'edgecolor', 'none');
A = label('(a)');
B = label('(b)');
C = label('(c)');
D = label('(d)');

A.Position(1) = 0;

B.Position(1) = sp4.Position(1) - 1;
B.Position(2) = sum(sp4.Position([2 4])) - 0.4;

A.Position(2) = B.Position(2);

C.Position(1) = A.Position(1);
C.Position(2) = sum(sp5.Position([2 4])) - 0.4;

D.Position(1) = B.Position(1);
D.Position(2) = C.Position(2);

asymptote = annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', '$\sim \ell / L$', 'interpreter', 'latex', 'edgecolor', 'none');
asymptote.Position(1) = sp6.Position(1) + 1.25;
asymptote.Position(2) = sp6.Position(2) + sp6.Position(4) * 0.5;

% Process function
function [correlationLength, l2norm] = render(nu, L, col)

  % Get file name
  source = sprintf('../data/processed/correlation/nu%s_L%d', nu, L);
  filenames = sortfiles(source, 'corr', []);
  Nfiles = length(filenames);

  % Preallocate for computing average correlation function
  Cbar = 0;

  % Preallocate for storing l2 norm
  l2Evolution = [];

  % Store time
  t = 1 : Nfiles;

  % Loop over files
  for nf = 1 : Nfiles
    C = load(strcat(source, '/', filenames{nf}));
    l2Evolution = [l2Evolution; sqrt(abs(C(5, 1)))];
  end

  % Averaging interval
  Nspan = min(50, Nfiles) : Nfiles;

  % Compute average
  for nf = Nspan
    C = load(strcat(source, '/', filenames{nf}));
    Cbar = Cbar + C(5, :);
  end

  % Normalize
  Cbar = Cbar / length(Nspan);

  % Get radial coordinate
  r = C(1, :);

  % Subplot 1 : Evoluation of the L2 norm
  subplot(3, 4, [9 10])
  plot(t, l2Evolution / L, '-', 'Color', col, 'DisplayName', sprintf('$L = %d$', L), 'HandleVisibility', 'off')
  hold on

  % Subplot 2 : Correlation functions
  subplot(3, 4, [3 4 7 8])

  if L == 150
    plot(r(2 : end) / L, Cbar(2 : end) / L^2, '-', 'Color', col, 'DisplayName', 'simulation')
  else
    plot(r(2 : end) / L, Cbar(2 : end) / L^2, '-', 'Color', col, 'HandleVisibility', 'off')
  end
  hold on

  [~,id] = find(Cbar < 0, 1);
  correlationLength = r(id);
  l2norm = mean(l2Evolution(Nspan));

end

function im = recolorImage(im)

  N1 = size(im, 1);
  N2 = size(im, 2);
  
  sp = 2;

  cx = 100;
  cy = 120;

  im = im((cx + 1) : sp : N1, (cy + 1) : sp : (N2 - cy), :);
  N1 = size(im, 1);
  N2 = size(im, 2);

  im = reshape(im, [N1 * N2 3]);
  im = reshape(im, [N1 N2 3]);
  
end