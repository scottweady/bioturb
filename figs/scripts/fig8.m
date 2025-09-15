
close all
addpath('utils')
fig = journal_figure([3.375 2], 2);

%% Compute velocity norm over parameter space and plot as a contour
nuspan = {'0.0625', '0.125', '0.25', '0.5', '0.75', '1.0'};
Lspan = 25 : 25 : 150;

Nnu = length(nuspan);
NL = length(Lspan);

U = zeros(Nnu + 1, NL);

for nnu = 1 : Nnu
  for nL = 1 : NL
    
    nu = str2double(nuspan{nnu});
    U(nnu + 1, nL) = (nu / 0.2) * velocityNorm(nuspan{nnu}, Lspan(nL));

  end
end

% Reform nu as double and create grid
nuspan = 5 * [0 0.0625 0.125 0.25 0.5 0.75 1.0];
[L, nu] = meshgrid(Lspan, nuspan);
contourf(nu, L, U, 16, 'EdgeColor', 'none'), hold on
nu = nu(2 : end, :); L = L(2 : end, :); U = U(2 : end, :);
scatter(nu(:), L(:), 500, U(:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2), hold on
plot([1 1], [25 150], 'w:')

%% Format
box on
axis tight

% Axis labels
xlabelHandle = xlabel('$\nu$');
ylabel('$L / \ell$')

xlim(5 * [0 1 + 1e-3])
ylim([25-1e-3 150+1e-3])

xticks(0 : 1 : 5)
yticks(25 : 25 : 150)

% Colormap
colormap(flipud(cmocean('matter')))
colorbarHandle = colorbar;
colorbarHandle.TickLabelInterpreter = 'latex';
colorbarHandle.Ticks = [0 0.7];
colorbarLabel = ylabel(colorbarHandle, '$\overline{||\mathbf{u}||_2}$', 'interpreter', 'latex', 'rotation', 0, 'FontSize', 18);

clim([0 0.7])

% Position
ax = gca;
ax.Units = 'inches';

ax.Position(1) = 0.8;
ax.Position(2) = 0.5;
ax.Position(3) = 0.7 * fig.PaperSize(1);
ax.Position(4) = 0.8 * fig.PaperSize(2);

xlabelHandle.Units = 'inches';
xlabelHandle.Position(2) = -0.25;

colorbarHandle.Units = 'inches';
colorbarHandle.Position(1) = sum(ax.Position([1 3])) + 0.1;
colorbarLabel.Units = 'inches';
colorbarLabel.Position(1) = colorbarHandle.Position(3) + 0.5;

nu_star = annotation('textbox', 'FontSize', 18, 'interpreter', 'latex', 'edgecolor', 'none', 'String', '$\nu^*$', 'units', 'inches', 'Color', 'w');
nu_star.Position(1) = 1.8;
nu_star.Position(2) = 3.15;

ax.ClippingStyle = 'rectangle';

%% Function to compute velocity norm
function l2norm_bar = velocityNorm(nu, L)

  % Get file name
  source = sprintf('../../data/processed/correlation_functions/nu%s_L%d', nu, L);
  C = load(strcat(source, '/u.dat'));
  l2norm = sqrt(C(:, 1));
  Nt = length(l2norm);

  % Averaging interval
  Nspan = min(Nt - 50, 150) : Nt;
  l2norm_bar = mean(l2norm(Nspan)) / L;


end