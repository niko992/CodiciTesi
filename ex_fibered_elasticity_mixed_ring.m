% EX_FIBERED_RING: solve the Poisson problem in one quarter of a ring, discretized with B-splines (non-isoparametric approach).

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_ring.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [3 4];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [1 2];
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Physical terms of fibered material
problem_data.Ef = 1e10;
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
hx = @(x, y, ind) sin(x).*cos(y).*(ind==1)+sin(x).*cos(y).*(ind==2);
hy = @(x, y, ind) sin(x-y)*(ind==1)+sin(x-y)*(ind==2);
problem_data.h       = @(x, y, ind) cat(1, ...
                reshape (hx (x,y,ind), [1, size(x)]), ...
                reshape (hy (x,y,ind), [1, size(x)]));
problem_data.g       = @(x, y, ind)...
    test_fibered_elasticity_ring_g_nmnn (x, y, ind);
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
method_data.degree     = [1 1];       % Degree of the splines
method_data.regularity = [0 0];       % Regularity of the splines
method_data.nsub       = [2 2];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_fibered_elasticity_mixed1 (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Ring_BSP_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

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

% Display errors of the computed solution in the L2 and H1 norm
[error_h1, error_l2] = ...
           sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)