
close all

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

% Reform as double
nuspan = 5 * [0 0.0625 0.125 0.25 0.5 0.75 1.0];
[L, nu] = meshgrid(Lspan, nuspan);

fig = journal_figure([6.75/2 2], 2);

el = 1; %particle length
b = 0.2; %particle diameter
mu = 1; %viscosity
Us = 1; %swimming speed
eta = 4 * pi * mu / log(2 * el / b); %slender body drag coefficient
sigma = eta * Us * el^2 / 8; %dipole strength

contourf(nu, L, U, 16, 'EdgeColor', 'none')
hold on
nu = nu(2 : end, :);
L = L(2 : end, :);
U = U(2 : end, :);
scatter(nu(:), L(:), 500, U(:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2)

hold on
plot([1 1], [25 150], 'w:')
box on
axis tight

xlabelHandle = xlabel('$\nu$');
ylabel('$L / \ell$')

colormap(flipud(cmocean('matter')))
colorbarHandle = colorbar;
colorbarHandle.TickLabelInterpreter = 'latex';
colorbarHandle.Ticks = [0 0.7];
colorbarLabel = ylabel(colorbarHandle, '$\overline{||\mathbf{u}||_2}$', 'interpreter', 'latex', 'rotation', 0, 'FontSize', 18);

clim([0 0.7])
xlim(5 * [0 1 + 1e-3])
ylim([25-1e-3 150+1e-3])

xticks(0 : 1 : 5)
yticks(25 : 25 : 150)

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

set(gca,'ClippingStyle','rectangle');

% Process function
function l2norm = velocityNorm(nu, L)

  % Get file name
  source = sprintf('../data/processed/correlation/nu%s_L%d', nu, L);
  filenames = sortfiles(source, 'corr', 'dat');
  Nfiles = length(filenames);

  if Nfiles == 0
    l2norm = nan;
    return
  end

  % Preallocate for storing l2 norm
  l2Evolution = zeros(Nfiles, 1);

  % Loop over files
  for nf = 1 : Nfiles
    C = load(strcat(source, '/', filenames{nf}));
    l2Evolution(nf) = sqrt(abs(C(5, 1)));
  end

  % Averaging interval
  Nspan = min(175, Nfiles) : Nfiles;

  if isempty(Nspan)
    l2norm = nan;
  else
    l2norm = mean(l2Evolution(Nspan)) / L;
  end


end