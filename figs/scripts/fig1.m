
% 3D flow around a slender force-free swimmer 
% See fig. 1 in manuscript

%% Setup

close all
fig = journal_figure([3.375 3.375], 2);

N = 16; %number of grid points in each direction
R = 1.5; %domain radius

xspan = 2 * R * linspace(-1, 1, N);
[xx, yy, zz] = meshgrid(xspan, xspan, xspan);

% Parameters
L = 1; % distance between point forces
a = 0.2; % radius of sphere
U = zeros(N, N, N, 3); % velocity field
p = [1; 0; 0]; % direction of swimming

%% Streamlines

% Quadrature points and weights for slender body integration
s = [-1 -0.5 0.5 1];
w = [1 8 8 1] / 9;
M = length(s);

% Compute velocity field
for i = 1 : N
  for j = 1 : N
    for k = 1 : N
      
      x = [xx(i, j, k) yy(i, j, k) zz(i, j, k)];
      M_ij = 0;

      for m = 1 : M
        sm = s(m);
        M_ij = M_ij + (sign(sm) * w(m)) * rpy(x, sm * p', a);
      end

      U(i, j, k, :) = reshape(M_ij * p, [1 1 3]);

    end
  end
end

hold on

% Streamline seeds
sx = R * linspace(-1, 1, 4);
sy = R * linspace(-1, 1, 4);
sz = 1.* R * [-1 1];

% Compute streamlines
[sx, sy, sz] = meshgrid(sx, sy, sz);
xyz = stream3(xx, yy, zz, U(:, :, :, 1), U(:, :, :, 2), U(:, :, :, 3), sx, sy, sz);
X = xyz{5};

% Plot streamlines as tubes
for i = 1 : length(xyz)
  X = xyz{i};
  if length(X) > 3 
    tubeplot(X', 0.05, 32);
  end
end

%% Swimmer body

s = linspace(0, pi); %centerline
i = ones(1, 16); %points around cross-section

% Swimmer surface points
x = [-L * i, L * i, L + a * cos(s - pi/2), L * i, -L * i, -L + a * cos(s + pi/2), -L * i]; 
z = [-a * i, -a * i, a * sin(s - pi/2), a * i, a * i, a * sin(s + pi/2), -a * i];
theta = linspace(0, 2 * pi, 128); %angular points

% Create surface points
x1 = z' * cos(theta);
x2 = z' * sin(theta);
x3 = x' * ones(size(theta));

% Plot swimmer surface
sh = surf(x3, x1, x2);
sh.EdgeColor = 'none';

%% Surface velocity

m = 3; %half number of points on the swimmer
i = ones(1, m); 
x = linspace(0, L, 2 * m);
z = [-a * i -a * i];

% Angular points
theta = linspace(0, 2 * pi, 17); 
theta = theta(1 : end - 1);
c = 1.01; %scale factor to move arrows off surface

% Surface points
x1 = c * z' * cos(theta);
x2 = c * z' * sin(theta);
x3 = c * x' * ones(size(theta));
x1 = x1(:); x2 = x2(:); x3 = x3(:);

% Surface velocity
us = 0 * x1(:); vs = 0 * x1(:); ws = 0 * x1(:) + 1;

% Plot
h = quiver3(x3, x2, x1, ws, vs, us, 1, 'Color', 0.95 * [1 1 1], 'LineWidth', 1, 'MaxHeadSize', 0.1);
shading interp

% Rotate arrows to point correctly
T = hgtransform('Parent', gca); %transform to rotate 90 degrees on x Axis
Rot = makehgtform('xrotate', pi/2); %set the transformation Matrix
set(T, 'Matrix', Rot); %set the parent of quiver to the new transform
set(h, 'Parent', T); %arrows should point correctly now

%% Swimmer arrow

% Render swimming direction arrow
arrow3d([-1.2 -1.4], [0 0], [0 0], 0.7, 0.025, 0.05, 0.2 * [1 1 1]);

%% Format

axis equal off
axis(0.97 * R * [-1 1 -1 1 -1 1])
box off
view(-20, 10)

% Colormap
load('fig1_cmap.mat', 'cmap')
colormap(cmap)

%% Functions
function M = rpy(x,y,a)
% Rotne-Prager-Yamakawa tensor
% x, y are 1x3 vectors
% a is the sphere radius

  rvec = x - y;
  rabs = sqrt(sum(rvec.^2,2));
  rhat = rvec ./ rabs;

  M = zeros(3);

  for i = 1 : 3
    for j = 1 : 3
      if rabs > 2 * a
        M(i,j) = (1./rabs + (2*a^2/3)./rabs.^3).*(i==j) + ...
          (1./rabs - 2*a^2./rabs.^3).*rhat(i).*rhat(j);
      else
        M(i,j) = (4/(3*a) - 3*rabs/(8*a^2)).*(i==j) + ...
          (rabs./(8*a^2)).*rhat(i).*rhat(j);
      end
    end
  end
  
end