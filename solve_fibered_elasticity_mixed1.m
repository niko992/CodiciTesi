%   TO BE MODIFIED


% SOLVE_FIBERED_ELASTICITY_MIXED1: Solve the fibered elasticity problem with a mixed formulation, and a B-spline discretization.
%
% Example to solve the problem
%
%    \int \psi(u,e) + \int \lambda:(\epsilon(u)-e) = f  in Omega = F((0,1)^n)
%                   gradu x n = g                       on Gamma_N
%                       u x n = h                       on Gamma_D
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
%    - Ef:
%    - a:
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
%  space:    space object that defines the discrete functions (see sp_vector)
%  sp_mul:   space object for the multiplier (see sp_scalar)
%  eigv:     the computed eigenvalues
%  eigf:     degrees of freedom of the associated eigenfunctions
%
% See also EX_MAXWELL_EIG_MIXED1_SQUARE for an example
%
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
% nurbs    = nrbdegelev (geometry.nurbs, degelev_u);
[~, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
% 
% nurbs    = nrbkntins (nurbs, nknots);
% geometry = geo_load (nurbs);

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
  rhs(dofs) = rhs(dofs) + op_f_v_tp (space_u.boundary(iside), msh.boundary(iside), gside);
end

rhs(1:nu) = rhs(1:nu) + F; 



% rhs_dir  = -A(int_dofs, drchlt_dofs)*u_drchlt;
% 
% rhs = [F(int_dofs) + rhs_dir;...
%     zeros(ndof_lam,1);...
%     zeros(ndof_lam,1)];
% K = [A(int_dofs, int_dofs),B(int_dofs,:),zeros(nintdofs, ndof_lam);...
%     B(int_dofs,:)',zeros(ndof_lam, ndof_lam),-C;...
%     zeros(ndof_lam, nintdofs), -C', D];
% 
% K = [A;  B; zeros(size(A,1), ndof_lam);
%      B';


int_dofs = setdiff (1:n, drchlt_dofs);

u(int_dofs) = K(int_dofs, int_dofs) \ rhs(int_dofs);

u = u(1:nu);

end
