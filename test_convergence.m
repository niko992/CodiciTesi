% EX_PLANE_STRAIN_SQUARE: solve the plane-strain problem on a square.
% clear
% clc
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
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [0 0];     % Regularity of the splines
method_data.nsub       = [128 128];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry_ref, msh_ref, space_ref, u_ref] = solve_fibered_elasticity_mixed1 (problem_data, method_data);
[euref, Xref, Yref] = test_refinement(u_ref, space_ref, geometry_ref);
figure(1)
surf (Xref,squeeze(euref(1,:,:)), squeeze(euref(2,:,:)))
%%
subdivision = [1,2,3,4,5];

for i = 1:length(subdivision)
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [1 1];     % Degree of the bsplines
method_data.regularity = [0 0];     % Regularity of the splines
method_data.nsub       = [2 2].^subdivision(i);     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_elasticity_fiber_material (problem_data, method_data);

[eu, X, Y] = test_refinement(u, space, geometry);
figure(i+1)
surf (X,squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))

% error_l2(i) = sp_l2_error(space_ref,msh_ref,u_ref-u_fine,problem_data.u_ex);
% error(i) = sum((eu-eu_ref).^2, [1 2 3])/sum(eu_ref.^2,[1,2,3]);
end

%% 
% figure
% semilogy(subdivision,error)
% legend('error in norm l2')
% title(['Convergence method P1 '])
% 
% fprintf('Convergence rate of the error in the l2 norm: %f\n', (log(error(end)) - log(error(1))) / (log(2^(-subdivision(end))) - log(2^(-subdivision(1)))));
