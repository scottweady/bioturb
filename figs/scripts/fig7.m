
close all
addpath('utils')
fig = journal_figure([3.375 2.6], 2);

%% Plot mean square displacement and orientation correlation
Lspan = [25 50 75 100 125 150];
NL = length(Lspan);

% Color map
cmap = 1 - linearrgbmap(1 - [60 0 143] / 256, NL + 1);
cmap = cmap(2 : end, :);

% Loop over box sizes
for nL = 1 : NL

    L = Lspan(nL);

    % Dilute
    nu = '0.125';
    msd = load(sprintf('../../data/processed/msd/nu%s_L%d/meanSquareDisplacement.dat', nu, L));
    t = 0 : length(msd) - 1;

    sp1 = subplot(2, 2, 1);

      plot(t, msd(1, :) / L^2, '-', 'Color', cmap(nL, :), 'HandleVisibility','off'), hold on
      xlim([0 100])
      ylim([0 0.8])
      yticks(0 : 0.2 : 0.8)
      xticklabels([])
      ylabel1 = ylabel('$\langle |\mathbf{x}_n(t) - \mathbf{x}_n(t_0)|^2\rangle$');

    sp3 = subplot(2, 2, 3);

      semilogy(t, msd(2, :) + msd(6, :) + msd(10, :), '-', 'Color', cmap(nL, :)), hold on
      xlim(sp1.XLim)
      ylim([1e-2 1])

      xlabel('$t - t_0$', 'FontSize', 18)
      ylabel3 = ylabel('$\langle \mathbf{p}_n(t) \cdot \mathbf{p}_n(t_0)\rangle$');

    % Concentrated
    nu = '1.0';
    msd = load(sprintf('../../data/processed/msd/nu%s_L%d/meanSquareDisplacement.dat', nu, L));
    t = 0 : length(msd) - 1;

    sp2 = subplot(2, 2, 2);

      plot(t, msd(1, :) / L^2, '-', 'Color', cmap(nL, :), 'HandleVisibility','off'), hold on
      xlim([0 20])
      xticklabels([])
      ylim(sp1.YLim)
      yticklabels([])

    sp4 = subplot(2, 2, 4);
      semilogy(t, msd(2, :) + msd(6, :) + msd(10, :), '-', 'Color', cmap(nL, :)), hold on
      xlabel('$t - t_0$', 'FontSize', 18)
      xlim(sp2.XLim);
      ylim(sp3.YLim);
      yticklabels([])

end

%% Format

sp1.Units = 'inches';
sp2.Units = 'inches';
sp3.Units = 'inches';
sp4.Units = 'inches';

sp1.FontSize = 18;
sp2.FontSize = 18;
sp3.FontSize = 18;
sp4.FontSize = 18;

colorbarHandle = colorbar;
colormap(cmap)
clim([min(Lspan) max(Lspan)])

sp1.Position(1) = 1;
sp3.Position(1) = sp1.Position(1);
sp3.Position(2) = 0.75;
sp2.Position(1) = sum(sp1.Position([1 3])) + 0.35;
sp4.Position(1) = sp2.Position(1);

sp4.Position(2 : 4) = sp3.Position(2 : 4);

sp1.Position(2) = sum(sp3.Position([2 4])) + 0.4;
sp2.Position(2 : 4) = sp1.Position(2 : 4);

colorbarHandle.Ticks = [25 150];

colorbarLabelHandle = ylabel(colorbarHandle, '$L / \ell$', 'interpreter', 'latex', 'rotation', 0, 'units', 'inches', 'FontSize', 18);
colorbarHandle.TickLabelInterpreter = 'latex';
colorbarHandle.Units = 'inches';
colorbarHandle.FontSize = 18;

colorbarHandle.Position(1) = sum(sp4.Position([1 3])) + 0.1;
colorbarHandle.Position(2) = sp4.Position(2);
colorbarHandle.Position(4) = sum(sp2.Position([2 4])) - sp4.Position(2); 

colorbarLabelHandle.Position(1) = 0.5;
colorbarLabelHandle.Position(2) = 2.1;

ylabel1.Position(1) = ylabel3.Position(1);

subplotLabel('(a)', sp1, 'northwest', [0 -0.55]);
subplotLabel('(b)', sp2, 'northwest', [0 -0.55]);
subplotLabel('(c)', sp3, 'southwest', [0 -0.15]);
subplotLabel('(d)', sp4, 'southwest', [0 -0.15]);

title(sp1, '$\nu = 0.625$', 'FontSize', 18)
title(sp2, '$\nu = 5$', 'FontSize', 18)
