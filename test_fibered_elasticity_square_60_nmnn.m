% TEST_FIBERED_ELASTICITY_SQUARE_60_NMNN: data function for Neumann boundary condition.

function g = test_fibered_elasticity_square_60_nmnn (x, y, ind)

  g = zeros ([2, size(x)]);
  switch (ind)
    case 2
      g(1,:,:) = cos(x) .* cos(y);
      g(2,:,:) = cos(x-y);
    case 4
      g(1,:,:) = -sin(x) .* sin(y);
      g(2,:,:) = -cos(x-y);
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end