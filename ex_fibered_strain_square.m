% EX_PLANE_STRAIN_SQUARE: solve the strain problem on a square.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrb4surf([0 0], [2*pi 0], [0 2*pi], [2*pi 2*pi]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Physical terms of fibered material
problem_data.Ef = 1e5;
problem_data.a = [1/2; sqrt(3)/2];

% Source and boundary terms
dmdx = @(x,y) -(1+sqrt(3))*sin(x).*cos(y) + 3*sin(x-y);
dmdy = @(x,y) -(1+sqrt(3))*cos(x).*sin(y) - 3*sin(x-y);
fx = @(x, y) -(-(problem_data.lambda_lame(x,y)+problem_data.mu_lame(x,y)).*cos(x).*sin(y)-...
    2*problem_data.mu_lame(x,y).*sin(x).*cos(y)+...
    problem_data.Ef*(dmdx(x,y)+sqrt(3)*dmdy(x,y))/16);
fy = @(x, y) -(-(problem_data.lambda_lame(x,y)+problem_data.mu_lame(x,y)).*sin(x).*cos(y)-...
    2*problem_data.mu_lame(x,y).*sin(x-y)+...
    problem_data.Ef*(sqrt(3)*dmdx(x,y)+3*dmdy(x,y))/16);
problem_data.f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
hx = @(x, y, ind) sin(2*pi)*cos(y)*(ind==2)+sin(x)*(ind==3)+sin(x)*cos(2*pi)*(ind==4);
hy = @(x, y, ind) sin(-y)*(ind==1)+sin(2*pi-y)*(ind==2)+sin(x)*(ind==3)+sin(x-2*pi)*(ind==4);
problem_data.h       = @(x, y, ind) cat(1, ...
                reshape (hx (x,y,ind), [1, size(x)]), ...
                reshape (hy (x,y,ind), [1, size(x)]));
% Exact solution (optional)
uxex = @(x,y) sin(x).*cos(y);
uyex = @(x,y) sin(x-y);
problem_data.uex = @(x, y) cat(1, ...
                reshape (uxex (x,y), [1, size(x)]), ...
                reshape (uyex (x,y), [1, size(x)]));
% Gradient of the exact solution (optional)
graduex11 = @(x,y) cos(x).*cos(y);
graduex21 = @(x,y) -sin(x).*sin(y);
graduex12 = @(x,y) cos(x - y);
graduex22 = @(x,y) -cos(x - y);
problem_data.graduex = @(x, y) cat(1, ...
                reshape (graduex11 (x,y), [1, size(x)]), ...
                reshape (graduex12 (x,y), [1, size(x)]),...
                reshape (graduex21 (x,y), [1, size(x)]), ...
                reshape (graduex22 (x,y), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [2 2].^3;     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_elasticity_fiber_material (problem_data, method_data);

% 4) POST-PROCESSING. 
% 4.1) Export to Paraview
output_file = 'plane_strain_square_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
    problem_data.lambda_lame, problem_data.mu_lame)

% 4.2) Plot in Matlab. Comparison with the exact solution.
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure
subplot (1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution'), axis equal tight
subplot (1,2,2)
eu2 = problem_data.uex (X, Y);
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
title ('Exact solution'), axis equal tight

figure
subplot (1,2,1)
surf (X,squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution')
subplot (1,2,2)
eu2 = problem_data.uex (X, Y);
surf(X, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
title ('Exact solution')

error_l2 = sp_l2_error (space, msh, u, problem_data.uex)
error_h1 = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)
%!demo
%! ex_plane_strain_square

%!test
%! problem_data.geo_name = nrb4surf([0 0], [1 0], [0 1], [1 1]);
%! problem_data.nmnn_sides   = [];
%! problem_data.press_sides  = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.symm_sides   = [];
%! E  =  1; nu = .3; 
%! problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
%! problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));
%! fx = @(x, y) -(-(problem_data.mu_lame(x, y)*3 + problem_data.lambda_lame(x, y)).*sin(2*pi*x).*sin(2*pi*y) + ...
%!      (problem_data.mu_lame(x, y) + problem_data.lambda_lame(x, y)).*cos(2*pi*x).*cos(2*pi*y))*(2*pi)^2;
%! fy = fx;
%! problem_data.f = @(x, y) cat(1, ...
%!                 reshape (fx (x,y), [1, size(x)]), ...
%!                 reshape (fy (x,y), [1, size(x)]));
%! problem_data.h       = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
%! uxex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
%! uyex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
%! problem_data.uex = @(x, y) cat(1, ...
%!                 reshape (uxex (x,y), [1, size(x)]), ...
%!                 reshape (uyex (x,y), [1, size(x)]));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [9 9];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_linear_elasticity (problem_data, method_data);
%! error_l2 = sp_l2_error (space, msh, u, problem_data.uex);
%! assert (msh.nel, 81)
%! assert (space.ndof, 288)
%! assert (error_l2, 2.60376176743492e-04, 1e-17)