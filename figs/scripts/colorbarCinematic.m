
close all

journal_figure([2 2], 2);

axis off
velocityColorbarHandle = colorbar('eastoutside');
colormap(magma)
velocityColorbarHandle.Units = 'inches';
velocityColorbarHandle.Ticks = [0 1];
% velocityColorbarHandle.TickLabels = {'$0$', '$|\mathbf{u}|_{\rm max}$'};
velocityColorbarHandle.TickLabels = {'$0$', '$0.025$'};
% velocityColorbarHandle.TickLabels = {'$0$', '$0.05$'};
velocityColorbarHandle.TickLabelInterpreter = 'latex';
velocityColorbarLabel = ylabel(velocityColorbarHandle, '$|\mathbf{u}|$', 'interpreter', 'latex', 'rotation', 0, 'FontSize', 18, 'units', 'inches');
velocityColorbarHandle.YAxisLocation = 'right';

velocityColorbarHandle.Position(1) = 2;
velocityColorbarHandle.Position(3) = 0.25;

velocityColorbarLabel.Position(1) = 0.5;