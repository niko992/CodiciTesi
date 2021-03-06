% SOLVE_FIBERED_ELASTICITY_MIXED1: Solve the fibered elasticity problem with a mixed formulation, and a B-spline discretization.
%
% Example to solve the problem
%
%    \int \psi(u,e) + \int \lambda:(\epsilon(u)-e) = f  in Omega = F((0,1)^n)
%                   sigma_an(u) x n = g                       on Gamma_N
%                                 u = h                       on Gamma_D
%
% with the variational mixed formulation and dividing \psi in deviatoric,
% anisotropic and volumetric parts
%
%    \int d_u(\psi_dev(u)) v + \int lambda \epsilon(v) = f   \forall v \in H^1(\Omega),
%    \int d_e(\psi_aniso(e)) v - \int \lambda g = 0          \forall g \in (L^2(\Omega)^(2x2)),  
%    \int \eta (\epsilon(u) - e) = 0                         \forall \eta \in (L^2(\Omega)^(2x2)).
%
% USAGE:
%
%  [geometry, msh, space, sp_mul, eigv, eigf] = 
%                  solve_maxwell_eig_mixed1 (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - lambda_lame:  first Lame' parameter
%    - mu_lame:      second Lame' parameter
%    - f:            source term
%    - h:            function for Dirichlet boundary condition
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - Ef:           inextensibility coefficient of the fibers
%    - a:            direction of the fibers
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space_u:  space object that defines the discrete functions (see sp_vector)
%  u:        displacement computed with the mixed formulation


function [geometry, msh, space_u, u] = ...
              solve_fibered_elasticity_mixed1 (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% load geometry
geometry = geo_load (geo_name);
% degelev  = max (degree - (geometry.nurbs.order-1), 0);
% nurbs    = nrbdegelev (geometry.nurbs, degelev);
[~, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
%

% Compute the mesh structure using the finest mesh
rule       = msh_gauss_nodes (nquad);
[qn, qw]   = msh_set_quad_nodes (zeta, rule);
msh        = msh_cartesian (zeta, qn, qw, geometry);

% Compute the space structures
[space_u, space_lam] = space_fibered_elasticity (geometry.nurbs.knots, nsub, degree, regularity, msh);
ndof_lam = nsub(1)*nsub(2)*3;
ndof_u = space_u.ndof;

% Assemble the matrices
A = op_su_ev_tp(space_u, space_u, msh, lambda_lame, mu_lame);
B = op_lam_ev_tp (space_lam, space_u, msh, ndof_lam);
C = op_lam_eta_tp (space_lam, space_lam, msh, ndof_lam);
D = op_j4lam_j4eta_tp (space_lam, space_lam, msh, Ef, a, ndof_lam);
F = op_f_v_tp (space_u, msh, f);

nu = size(A, 1);
nl = size(B, 2);

K = [A, B, sparse(nu, nl);
     B', sparse(nl, nl), -C;
    sparse(nl, nu), -C', D];

n = size(K, 1);

% Apply symmetry conditions
u = zeros (space_u.ndof, 1);
symm_dofs = [];
for iside = symm_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  for idim = 1:msh.rdim
    normal_comp(idim,:) = reshape (msh_side.normal(idim,:,:), 1, msh_side.nqn*msh_side.nel);
  end

  parallel_to_axes = false;
  for ind = 1:msh.rdim
    ind2 = setdiff (1:msh.rdim, ind);
    if (all (all (abs (normal_comp(ind2,:)) < 1e-10)))
      symm_dofs = union (symm_dofs, space_u.boundary(iside).comp_dofs{ind});
      parallel_to_axes = true;
      break
    end
  end
  if (~parallel_to_axes)
    error ('solve_linear_elasticity: We have only implemented the symmetry condition for boundaries parallel to the axes')
  end

end

% Apply Dirichlet boundary conditions
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space_u, msh, h, drchlt_sides);
u = zeros(n, 1);
u(drchlt_dofs) = u_drchlt;

rhs = -K * u;

% Apply Neumann boundary conditions
for iside = nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
  gside = @(varargin) g(varargin{:},iside);
  dofs = space_u.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) +...
      op_f_v_tp (space_u.boundary(iside), msh.boundary(iside), gside);
end

rhs(1:nu) = rhs(1:nu) + F; 

int_dofs = setdiff (1:n, [drchlt_dofs;symm_dofs]);

u(int_dofs) = K(int_dofs, int_dofs) \ rhs(int_dofs);

u = u(1:nu);

end
