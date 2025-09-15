
close all
fig = journal_figure([6.75 2], 2);

% Load in continuum data
load('../data/fig4/unorm.mat');
load('../data/fig4/correlation.mat', 'r', 'C');

tc = 2.3453; % characteristic time scale
umag = load('../data/fig4/umag.dat');
umag = reshape(umag, [256 256 256]) / tc;

% Load in discrete data
load('../data/fig4/L175.mat');

% Color scheme
c_continuum = [0.1375 0.0685 0.3191];
c_discrete = [0.71 0.21 0.48];

sp1 = subplot(1, 3, 1);

  N = size(umag, 1);
  
  sliceHandle = slice(umag, 1 : 4 : N, 1 : 4 : N, 1 : 4 : N);
  
  for n = 1 : length(sliceHandle)
  
    sliceHandle(n).AlphaData = (0.5 + (sliceHandle(n).CData / max(umag(:))).^1.5) / 1.5;
    sliceHandle(n).FaceAlpha = 'interp';
  
  end
  
  shading interp
  axis equal tight
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  colorbarHandle = colorbar;
  colormap(magma)
  view(60, 15)

  clim([0 0.03])
  
  xticks([]), yticks([]), zticks([])

sp2 = subplot(1, 3, 2);

  id = 1 : 10 : length(unorm);
  plot(unorm(id, 1), unorm(id, 2), 'Color', c_continuum, 'DisplayName', 'continuum')
  hold on
  plot(t_175, l2_175, 'Color', c_discrete, 'DisplayName', 'discrete')
  hold on

  xlim([0 300])
  ylim([0 0.025])
  
  xlabel('$t$')
  ylabel('$||\mathbf{u}||_2$')

  legend('FontSize', 18, 'location', 'southeast')

sp3 = subplot(1, 3, 3);

  plot(r, C / C(1), 'Color', c_continuum, 'DisplayName', 'continuum'), hold on
  plot(r_175, C_175 / C_175(1), 'Color', c_discrete, 'DisplayName', 'discrete')
  plot([0 0.5], [0 0], 'k--', 'HandleVisibility', 'off')
  
  xlim([0 0.5])
  ylim([-0.1 1])
  
  xlabel('$r / L$')
  ylabel('$\overline{{\rm Corr}[\mathbf{u}]}$')

  legend('FontSize', 18, 'location', 'southwest')

sp1.Units = 'inches';
sp2.Units = 'inches';
sp3.Units = 'inches';

sp1.Position(1) = 0.25;
sp2.Position(1) = sum(sp1.Position([1 3])) + 2.25;
sp3.Position(1) = sum(sp2.Position([1 3])) + 1.5;

sp2.Position(2) = 0.75;
sp2.Position(3) = 0.26 * fig.PaperSize(1);
sp2.Position(4) = 0.7 * fig.PaperSize(2);

sp3.Position(2 : 4) = sp2.Position(2 : 4);

colorbarHandle.Units = 'inches';
colorbarHandle.TickLabelInterpreter = 'latex';
colorbarHandle.Limits = [0 0.03];
colorbarHandle.Ticks = 0 : 0.01 : 0.03;
colorbarHandle.FontSize = 18;
colorbarHandle.Title.String = '$|\mathbf{u}|$';
colorbarHandle.Title.Interpreter = 'latex';

colorbarHandle.Position(1) = sum(sp1.Position([1 3])) + 0.25;
colorbarHandle.Position(4) = (2 / 3) * sp1.Position(4);
colorbarHandle.Position(2) = sp1.Position(2) + (1 / 6) * sp1.Position(4);

ax = gca;

label = @(s) annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', s, 'interpreter', 'latex', 'edgecolor', 'none');
a = label('(a)');
b = label('(b)');
c = label('(c)');

a.Position(1) = 0;
b.Position(1) = sp2.Position(1) - 1;

c.Position(1) = sp3.Position(1) - 0.75;
c.Position(2) = sum(sp3.Position([2 4])) - 0.25;

a.Position(2) = c.Position(2);
b.Position(2) = c.Position(2);
