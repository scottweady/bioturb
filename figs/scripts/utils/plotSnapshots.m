
function plotSnapshots(nu)

% Load images
% polar = @(L) resizeOrderParameter(imread(sprintf('../snapshots/polar/nu%s_L%d.png', nu, L)));
% nematic = @(L) resizeOrderParameter(imread(sprintf('../snapshots/nematic/nu%s_L%d.png', nu, L)));

subplot_list = cell(2, 3);

subplot_list{1, 1} = subplot(3, 4, 1);
imshow(polar(nu, 50))
subplot_list{1, 2} = subplot(3, 4, 2);
imshow(polar(nu, 100))
subplot_list{1, 3} = subplot(3, 4, 3);
imshow(polar(nu, 150))

subplot_list{2, 1} = subplot(3, 4, 5);
imshow(nematic(nu, 50))
subplot_list{2, 2} = subplot(3, 4, 6);
imshow(nematic(nu, 100))
subplot_list{2, 3} = subplot(3, 4, 7);
imshow(nematic(nu, 150))

for i = 1 : 2
  for j = 1 : 3
    subplot_list{i, j}.Units = 'inches';
    subplot_list{i, j}.Position([3 4]) = 2.95 * [2096 2160] / 2160;
  end
end


for i = 1 : 2
  subplot_list{i, 1}.Position(1) = 0.125;
  subplot_list{i, 2}.Position(1) = sum(subplot_list{i, 1}.Position([1 3]));
  subplot_list{i, 3}.Position(1) = sum(subplot_list{i, 2}.Position([1 3]));
end

for j = 1 : 3
  subplot_list{2, j}.Position(2) = 0.5;
  subplot_list{1, j}.Position(2) = sum(subplot_list{2, j}.Position([2 4])) + 0.125;
end

subplotLabel('(a)', subplot_list{1, 1}, 'northwest', [0 -0.55]);
subplotLabel('(b)', subplot_list{2, 1}, 'northwest', [0 -0.55]);

boxsize = @(L) annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', strcat('$L/\ell = ', num2str(L), '$'), 'interpreter', 'latex', 'edgecolor', 'none');

L50 = boxsize(50);
L50.Position(1) = subplot_list{2, 1}.Position(1) + 0.3 * subplot_list{2, 1}.Position(3);
L50.Position(2) = -0.1;

L100 = boxsize(100);
L100.Position(1) = subplot_list{2, 2}.Position(1) + 0.3 * subplot_list{2, 2}.Position(3);
L100.Position(2) = -0.1;

L150 = boxsize(150);
L150.Position(1) = subplot_list{2, 3}.Position(1) + 0.3 * subplot_list{2, 3}.Position(3);
L150.Position(2) = -0.1;

end

function im = polar(nu, L)
  try
    im = resizeOrderParameter(imread(sprintf('../snapshots/polar/nu%s_L%d.png', nu, L)));
  catch
    im = [];
  end
end

function im = nematic(nu, L)
  try
    im = resizeOrderParameter(imread(sprintf('../snapshots/nematic/nu%s_L%d.png', nu, L)));
  catch
    im = [];
  end
end