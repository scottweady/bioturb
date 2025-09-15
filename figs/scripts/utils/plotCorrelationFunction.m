% Process function
function plotCorrelationFunction(nu, L, col, i, N0, R)

  % Get file name
  source = sprintf('../../data/processed/correlation_functions/nu%s_L%d', nu, L);

  r = load(strcat(source, '/r.dat'));
  r = r(1, :);

  switch i
    case 2
      C = load(strcat(source, '/c.dat'));
    case 3
      C = load(strcat(source, '/n.dat'));
    case 4
      C = load(strcat(source, '/Q.dat'));
  end

  Nt = size(C, 1);

  % Averaging interval
  Nspan = min(N0, Nt) : Nt;
  Cbar = sum(C(Nspan, :), 1) / length(Nspan) - (i == 2);

  % Subplot 2 : Correlation functions
  subplot(3, 4, 7 + i)
  plot(r(2 : end) / R, Cbar(2 : end), 'o-', 'Color', col, 'DisplayName', sprintf('$L = %d$', L))
  hold on

end