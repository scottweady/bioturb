
function [r, C] = correlationFunction(U, L)

  N = size(U, 1);
  U_h = fft3(U);

  if ndims(U) == 3
    Pu_h = U_h .* conj(U_h);
  else
    Pu_h = sum(U_h .* conj(U_h), 4 : ndims(U));
  end
  
  nshift = N;

  Pu = ifft3(Pu_h);
  Pu = pad(Pu);
  Pu = ifftshift(Pu);

  count = zeros(3 * N, 1);
  r = zeros(3 * N, 1);
  C = zeros(3 * N, 1);

  nmax = 0;

  for nx = 1 : 2 * N
    for ny = 1 : 2 * N
      for nz = 1 : 2 * N

        nr = 1 + round(sqrt( (nx - nshift)^2 + (ny - nshift)^2 + (nz - nshift)^2 ));

        r(nr) = (nr - 1) * L / N;
        C(nr) = C(nr) + Pu(nx, ny, nz) / N^3;
        count(nr) = count(nr) + 1;
        nmax = max(nmax, nr);

      end
    end
  end

  r = r(1 : nmax);
  C = (C(1 : nmax) ./ count(1 : nmax));

end


function P_pad = pad(P)

  shp = size(P);
  Nx = shp(1);
  Ny = shp(2);
  Nz = shp(3);

  P_pad = zeros(2 * shp);
  
  i = 1 : Nx;
  j = 1 : Ny;
  k = 1 : Nz;
  
  P_pad(i, j, k) = P;
  P_pad(Nx + i, j, k) = P;
  P_pad(i, Ny + j, k) = P;
  P_pad(i, j, Nz + k) = P;

  P_pad(Nx + i, Ny + j, k) = P;
  P_pad(Nx + i, j, Nz + k) = P;
  P_pad(i, Ny + j, Nz + k) = P;
  P_pad(Nx + i, Ny + j, Nz + k) = P;
  
end

function U_h = fft3(U)
    U_h = fft(fft(fft(U, [], 1), [], 2), [], 3);
end

function U = ifft3(U_h)
    U = real(ifft(ifft(ifft(U_h, [], 1), [], 2), [], 3));
end

