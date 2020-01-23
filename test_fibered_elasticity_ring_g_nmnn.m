% TEST_FIBERED_ELASTICITY_RING_G_NMNN: data function for Neumann boundary condition.
function g = test_fibered_elasticity_ring_g_nmnn (x, y, ind, graduex, a, lambda_lame, mu_lame, Ef)

  
  gradu = reshape(graduex(x, y), [2, 2, numel(x)]);
  
  strain = 0.5 * (gradu + permute(gradu, [2, 1, 3]));
  
  lambda = reshape(lambda_lame(x, y), [1, 1, numel(x)]);
  mu = reshape(repmat(reshape(mu_lame(x, y), [1, numel(x)]), [4, 1]), [2, 2, numel(x)]);
  
  lambda_div = (gradu(1, 1, :) + gradu(2, 2, :)) .* lambda;
  
  
  stress = 2.0 * mu .* strain;
  stress(1, 1, :) = stress(1, 1, :) + lambda_div;
  stress(2, 2, :) = stress(2, 2, :) + lambda_div;
  
  A = repmat(a * a', [1, 1, numel(x)]);
  J4 = sum(A .* strain, [1, 2]);
  
  stress(1, 1, :) = stress(1, 1, :) + Ef * J4 .* A(1, 1, :); 
  stress(1, 2, :) = stress(1, 2, :) + Ef * J4 .* A(1, 2, :); 
  stress(2, 1, :) = stress(2, 1, :) + Ef * J4 .* A(2, 1, :); 
  stress(2, 2, :) = stress(2, 2, :) + Ef * J4 .* A(2, 2, :); 
  
  n = zeros([2, size(x)]);
  
  switch (ind)
    case 3
      n(2, :) = -1;
    case 4
      n(1, :) = -1;
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  
  g = zeros([2, size(x)]);
  
  g(1, :) = sum(reshape(stress(1, :, :), [2, numel(x)]) .* n, 1); 
  g(2, :) = sum(reshape(stress(2, :, :), [2, numel(x)]) .* n, 1); 


end