function g = test_fibered_elasticity_square45_g_nmnn (x, y, ind)

  g = zeros ([2, size(x)]);
  switch (ind)
    case 2
      g(1,:,:) = 0.*x;
      g(2,:,:) = 0.*y;
    case 4
      g(1,:,:) = 0.*x;
      g(2,:,:) = 0.*y;
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end