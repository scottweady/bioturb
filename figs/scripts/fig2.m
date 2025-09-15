
close all

% Initialize figure
fig = journal_figure([6.75 3.325], 2);

% Volume fraction
nu = '0.125';

% Preallocate for storing l2 norm
l2norm = [];

% Span of box sizes
Lspan = 25 : 25 : 175;

% Number of box size simulations
NL = length(Lspan);

% Color map
colors = cmocean('phase', 5);
colors = colors(2 : end - 1, :);

for i = 2 : 4

  col = colors(i - 1, :);
  cmap = 1 - linearrgbmap(1 - col, NL + 1);
  cmap = cmap(2 : end, :);
  
  % Subplot 1, 2 : evolution of l2 norm and correlation function
  for nL = 1 : NL
  
    % Get box size
    L = Lspan(nL);
    
    % Compute l2 norm and plot correlation functions
    render(nu, L, cmap(nL, :), i);
  
  end

end

sp = cell(3, 3);

for i = 2 : 4

  sp{i - 1} = subplot(3, 4, 7 + i);
  
  plot([0 1], [0 0], 'k--', 'HandleVisibility', 'off')
  xlabel('$r / L$')
  
  switch i

    case 2
        sp{i - 1, 3} = ylabel('Corr[$c''$]');
        ylim([-0.005 0.02])
        yticks(0 : 0.01 : 0.02)
        yticklabels({'0', '0.01', '0.02'})
        
    case 3
        sp{i - 1, 3} = ylabel('Corr[$\mathbf{n}$]');
        ylim([-0.01 0.04])
        yticks(0 : 0.02 : 0.04)
        yticklabels({'0', '0.02', '0.04'})
        
    case 4
        sp{i - 1, 3} = ylabel('Corr[$\mathbf{Q}$]');
        ylim([-0.01 0.09])
        yticks(0 : 0.04 : 0.08)
        yticklabels({'0','0.04','0.08'})
    
  end

  
  xlim([0 0.5])
  xticks(0 : 0.1 : 0.5)
  
  col = colors(i - 1, :);
  
  cmap = 1 - linearrgbmap(1 - col, NL + 1);
  cmap = cmap(2 : end, :);

  colormap(sp{i - 1}, cmap)
  sp{i - 1, 2} = colorbar;
  clim([25 175])
  
  sp{i - 1, 2}.TickLabelInterpreter = 'latex';
  sp{i - 1, 2}.FontSize = 18;
  sp{i - 1, 2}.Ticks = [25 175];

end

for i = 1 : 3
  for j = 1 : 2
    sp{i, j}.Units = 'inches';
  end
end

sp{1, 1}.Position(3) = 0.21 * fig.PaperSize(1);
sp{1, 1}.Position(1) = fig.PaperSize(1) - sp{1, 1}.Position(3) - 1;

sp{2, 1}.Position([1 3]) = sp{1, 1}.Position([1 3]);
sp{3, 1}.Position([1 3]) = sp{1, 1}.Position([1 3]);

axes(sp{1, 1})
xlabel([])
xticklabels([])

axes(sp{2, 1})
xlabel([])
xticklabels([])

sp{1, 2}.Title.String = '$L/\ell$';
sp{1, 2}.Title.Interpreter = 'latex';

sp{3, 1}.Position(2) = 0.65;
sp{3, 1}.Position(4) = 0.26 * fig.PaperSize(2);

sp{2, 1}.Position([1 3 4]) = sp{3, 1}.Position([1 3 4]);
sp{1, 1}.Position([1 3 4]) = sp{3, 1}.Position([1 3 4]);

sp{2, 1}.Position(2) = sum(sp{3, 1}.Position([2 4])) + 0.225;
sp{1, 1}.Position(2) = sum(sp{2, 1}.Position([2 4])) + 0.225;

sp{1, 2}.Position(1) = sum(sp{1, 1}.Position([1 3])) + 0.1;
sp{1, 2}.Position([2 4]) = sp{1, 1}.Position([2 4]);
sp{2, 2}.Position([2 4]) = sp{2, 1}.Position([2 4]);
sp{3, 2}.Position([2 4]) = sp{3, 1}.Position([2 4]);

sp{2, 2}.Position(1) = sp{1, 2}.Position(1);
sp{3, 2}.Position(1) = sp{1, 2}.Position(1);

sp{1, 1}.FontSize = 18;
sp{2, 1}.FontSize = 18;
sp{3, 1}.FontSize = 18;

sp{1, 3}.Position(1) = sp{3, 3}.Position(1);
sp{2, 3}.Position(1) = sp{3, 3}.Position(1);

label = @(s) annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', s, 'interpreter', 'latex', 'edgecolor', 'none');

C = label('(c)');
D = label('(d)');
E = label('(e)');

C.Position(1) = sum(sp{1, 1}.Position([1 3])) - 0.45;
C.Position(2) = sum(sp{1, 1}.Position([2 4])) - 0.7;

D.Position(1) = sum(sp{2, 1}.Position([1 3])) - 0.45;
D.Position(2) = sum(sp{2, 1}.Position([2 4])) - 0.7;

E.Position(1) = sum(sp{3, 1}.Position([1 3])) - 0.45;
E.Position(2) = sum(sp{3, 1}.Position([2 4])) - 0.7;

% Load images
polar = @(L) recolorImage(imread(sprintf('../snapshots/polar/nu0.125_L%d.png', L)));
nematic = @(L) recolorImage(imread(sprintf('../snapshots/nematic/nu0.125_L%d.png', L)));

sp_config = cell(2, 3);

sp_config{1, 1} = subplot(3, 4, 1);
imshow(polar(50))
sp_config{1, 2} = subplot(3, 4, 2);
imshow(polar(100))
sp_config{1, 3} = subplot(3, 4, 3);
imshow(polar(150))

sp_config{2, 1} = subplot(3, 4, 5);
imshow(nematic(50))
sp_config{2, 2} = subplot(3, 4, 6);
imshow(nematic(100))
sp_config{2, 3} = subplot(3, 4, 7);
imshow(nematic(150))

for i = 1 : 2
  for j = 1 : 3
    sp_config{i, j}.Units = 'inches';
    sp_config{i, j}.Position([3 4]) = 2.95 * [2096 2160] / 2160;
  end
end


for i = 1 : 2
  sp_config{i, 1}.Position(1) = 0.125;
  sp_config{i, 2}.Position(1) = sum(sp_config{i, 1}.Position([1 3]));
  sp_config{i, 3}.Position(1) = sum(sp_config{i, 2}.Position([1 3]));
end

for j = 1 : 3
  sp_config{2, j}.Position(2) = 0.5;
  sp_config{1, j}.Position(2) = sum(sp_config{2, j}.Position([2 4])) + 0.125;
end

A = label('(a)');
A.Position(1) = sp_config{1, 1}.Position(1);
A.Position(2) = sum(sp_config{1, 1}.Position([2 4])) - 0.5;

B = label('(b)');
B.Position(1) = sp_config{2, 1}.Position(1);
B.Position(2) = sum(sp_config{2, 1}.Position([2 4])) - 0.5;

boxsize = @(L) annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', strcat('$L/\ell = ', num2str(L), '$'), 'interpreter', 'latex', 'edgecolor', 'none');

L50 = boxsize(50);
L50.Position(1) = sp_config{2, 1}.Position(1) + 0.3 * sp_config{2, 1}.Position(3);
L50.Position(2) = -0.1;

L100 = boxsize(100);
L100.Position(1) = sp_config{2, 2}.Position(1) + 0.3 * sp_config{2, 2}.Position(3);
L100.Position(2) = -0.1;

L150 = boxsize(150);
L150.Position(1) = sp_config{2, 3}.Position(1) + 0.3 * sp_config{2, 3}.Position(3);
L150.Position(2) = -0.1;

% Process function
function render(nu, L, col, i)

  % Get file name
  source = sprintf('../data/processed/correlation/nu%s_L%d', nu, L);
  filenames = sortfiles(source, 'corr', []);
  Nfiles = length(filenames);

  % Initialize for computing average correlation function
  Cbar = 0;

  % Averaging interval
  Nspan = min(150, Nfiles) : Nfiles;

  % Compute average
  for nf = Nspan
    C = load(strcat(source, '/', filenames{nf}));
    Cbar = Cbar + C(i, :);
  end

  % Normalize
  Cbar = Cbar / length(Nspan);

  % Get radial coordinate
  r = C(1, :);

  % If concentration, remove mean
  Cbar = Cbar - (i == 2);

  % Subplot 2 : Correlation functions
  subplot(3, 4, 7 + i)
  plot(r(2 : end) / L, Cbar(2 : end), 'o-', 'Color', col, 'DisplayName', sprintf('$L = %d$', L))
  hold on

end

function im = recolorImage(im)

  N1 = size(im, 1);
  N2 = size(im, 2);
  
  sp = 1;

  cx = 100;
  cy = 120;

  im = im((cx + 1) : sp : N1, (cy + 1) : sp : (N2 - cy), :);
  N1 = size(im, 1);
  N2 = size(im, 2);

  im = reshape(im, [N1 * N2 3]);
  im(sum(im, 2) == 0, :) = 255;
  im = reshape(im, [N1 N2 3]);
  im = uint8(round(255 * (double(im) / 180).^2));

end
