function im = resizeOrderParameter(im)
% Resize, scale, and recolor order parameter image

  N1 = size(im, 1);
  N2 = size(im, 2);
  
  sp = 2;

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