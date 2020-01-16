% EX_PLANE_STRAIN_SQUARE: solve the plane-strain problem on a square
% for anisotropic material with fibers at 45 degree.
clear
clc
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrb4surf([0 0], [1 0], [0 1], [1 1]);

% Type of boundary conditions
problem_data.nmnn_sides   = [2 4];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
problem_data.symm_sides   = [1 3];

% Physical parameters
E  =  1; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Physical terms of fibered material
problem_data.Ef = 1e5;
problem_data.a = [sqrt(2)/2; sqrt(2)/2];

% Source and boundary terms
fx = @(x, y) -1+0.*x;
fy = @(x, y) 1+0.*y;
problem_data.f       = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
hx = @(x, y, ind) 0.*x*(ind==1);
hy = @(x, y, ind) 0.*y*(ind==3);
problem_data.h       = @(x, y, ind) cat(1, ...
                reshape (hx (x,y,ind), [1, size(x)]), ...
                reshape (hy (x,y,ind), [1, size(x)]));
problem_data.g = @(x, y, ind) test_fibered_elasticity_square45_g_nmnn(x, y, ind);
%%
% Reference solution (optional)
uxex = @(x,y) 0.*x;
uyex = @(x,y) 0.*y;
uex = @(x, y) cat(1, ...
                reshape (uxex (x,y), [1, size(x)]), ...
                reshape (uyex (x,y), [1, size(x)]));
graduex11 = @(x,y) 0.*x;
graduex21 = @(x,y) 0.*x;
graduex12 = @(x,y) 0.*x;
graduex22 = @(x,y) 0.*x;
graduex = @(x, y) cat(1, ...
                reshape (graduex11 (x,y), [1, size(x)]), ...
                reshape (graduex12 (x,y), [1, size(x)]),...
                reshape (graduex21 (x,y), [1, size(x)]), ...
                reshape (graduex22 (x,y), [1, size(x)]));
%%
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [1 1];     % Degree of the bsplines
method_data.regularity = [0 0];     % Regularity of the splines
method_data.nsub       = [32 32];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

method_data2.degree     = [3 3];     % Degree of the bsplines
method_data2.regularity = [0 0];     % Regularity of the splines
method_data2.nsub       = [128 128];     % Number of subdivisions
method_data2.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry_ref, msh_ref, space_ref, u_ref] = solve_elasticity_fiber_material (problem_data, method_data2);
[geometry, msh, space, u] = solve_fibered_elasticity_mixed1 (problem_data, method_data);
[geometry2, msh2, space2, u2] = solve_fibered_elasticity_mixed2 (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Ring_BSP_Deg3_Reg2_Sub9';

%% Plot
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

[eu2, F2] = sp_eval (u2, space2, geometry2, vtk_pts);
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[eu_ref, F_ref] = sp_eval (u_ref, space_ref, geometry_ref, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
[X2, Y2]  = deal (squeeze(F2(1,:,:)), squeeze(F2(2,:,:)));
[X_ref,Y_ref]  = deal (squeeze(F_ref(1,:,:)), squeeze(F_ref(2,:,:)));

figure
subplot (1,3,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution with locking phenomena'), axis equal tight
subplot (1,3,2)
quiver (X2, Y2, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
title ('Numerical solution without locking phenomena'), axis equal tight
subplot (1,3,3)
quiver (X_ref, Y_ref, squeeze(eu_ref(1,:,:)), squeeze(eu_ref(2,:,:)))
title ('Reference solution'), axis equal tight

figure
subplot (1,3,1)
surf (X,squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution with locking phenomena')
subplot (1,3,2)
surf (X2,squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
title ('Numerical solution with locking phenomena')
subplot (1,3,3)
surf (X_ref,squeeze(eu_ref(1,:,:)), squeeze(eu_ref(2,:,:)))
title ('Exact solution')

figure
surf (X,squeeze(eu(1,:,:)-eu_ref(1,:,:)), squeeze(eu(2,:,:)-eu_ref(2,:,:)))