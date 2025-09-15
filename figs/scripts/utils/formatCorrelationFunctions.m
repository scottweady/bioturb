function formatCorrelationFunctions(fig, subplot_list)

for i = 1 : 3
  for j = 1 : 2
    subplot_list{i, j}.Units = 'inches';
    subplot_list{i, j}.FontSize = 18;
  end
end


subplot_list{1, 1}.Position(3) = 0.21 * fig.PaperSize(1);
subplot_list{1, 1}.Position(1) = fig.PaperSize(1) - subplot_list{1, 1}.Position(3) - 1;

subplot_list{2, 1}.Position([1 3]) = subplot_list{1, 1}.Position([1 3]);
subplot_list{3, 1}.Position([1 3]) = subplot_list{1, 1}.Position([1 3]);

axes(subplot_list{1, 1})
xlabel([])
xticklabels([])

axes(subplot_list{2, 1})
xlabel([])
xticklabels([])

subplot_list{1, 2}.Title.String = '$L/\ell$';
subplot_list{1, 2}.Title.Interpreter = 'latex';

subplot_list{3, 1}.Position(2) = 0.65;
subplot_list{3, 1}.Position(4) = 0.26 * fig.PaperSize(2);

subplot_list{2, 1}.Position([1 3 4]) = subplot_list{3, 1}.Position([1 3 4]);
subplot_list{1, 1}.Position([1 3 4]) = subplot_list{3, 1}.Position([1 3 4]);

subplot_list{2, 1}.Position(2) = sum(subplot_list{3, 1}.Position([2 4])) + 0.225;
subplot_list{1, 1}.Position(2) = sum(subplot_list{2, 1}.Position([2 4])) + 0.225;

subplot_list{1, 2}.Position(1) = sum(subplot_list{1, 1}.Position([1 3])) + 0.1;

for i = 1 : 3
  subplot_list{i, 2}.Position([2 4]) = subplot_list{i, 1}.Position([2 4]);
end

subplot_list{2, 2}.Position(1) = subplot_list{1, 2}.Position(1);
subplot_list{3, 2}.Position(1) = subplot_list{1, 2}.Position(1);

subplot_list{1, 3}.Position(1) = subplot_list{3, 3}.Position(1);
subplot_list{2, 3}.Position(1) = subplot_list{3, 3}.Position(1);

subplotLabel('(c)', subplot_list{1, 1}, 'northeast', [-0.45 -0.7]);
subplotLabel('(d)', subplot_list{2, 1}, 'northeast', [-0.45 -0.7]);
subplotLabel('(e)', subplot_list{3, 1}, 'northeast', [-0.45 -0.7]);