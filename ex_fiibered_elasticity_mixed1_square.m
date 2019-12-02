% EX_FIBERED_ELASTICITY_MIXED1_SQUARE: solve the fibered elasticity problem with the mixed form1.
clear;
clc;
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrb4surf([0 0], [1 0], [0 1], [1 1]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
E  =  1; nu = 0.3; 
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
hx = @(x, y, ind) sin(1)*cos(y)*(ind==2)+sin(x)*(ind==3)+sin(x)*cos(1)*(ind==4);
hy = @(x, y, ind) sin(-y)*(ind==1)+sin(1-y)*(ind==2)+sin(x)*(ind==3)+sin(x-1)*(ind==4);
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
graduex12 = @(x,y) -sin(x).*sin(y);
graduex21 = @(x,y) cos(x - y);
graduex22 = @(x,y) -cos(x - y);
problem_data.graduex = @(x, y) cat(1, ...
                reshape (graduex11 (x,y), [1, size(x)]), ...
                reshape (graduex21 (x,y), [1, size(x)]),...
                reshape (graduex12 (x,y), [1, size(x)]), ...
                reshape (graduex22 (x,y), [1, size(x)]));
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [1 1]; % Degree of the bsplines
method_data.regularity = [0 0]; % Regularity of the splines
method_data.nsub       = [16 16]; % Number of subdivisions
method_data.nquad      = [8 8]; % Points for the Gaussian quadrature rule


% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = ...
                    solve_fibered_elasticity_mixed1 (problem_data, method_data);
                
vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
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

figure
surf (X,squeeze(abs(eu(1,:,:)-eu2(1,:,:))), squeeze(abs(eu(2,:,:)-eu2(2,:,:))))

error_l2 = sp_l2_error (space, msh, u, problem_data.uex)