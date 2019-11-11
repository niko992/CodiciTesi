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
method_data.degree     = [1 1]; % Degree of the bsplines
method_data.regularity = [0 0]; % Regularity of the splines
method_data.nsub       = [2 2]; % Number of subdivisions
method_data.nquad      = [4 4]; % Points for the Gaussian quadrature rule


% 3) CALL TO THE SOLVER
[geometry, msh, space, sp_mul, u] = ...
                    solve_fibered_elasticity_mixed1 (problem_data, method_data);