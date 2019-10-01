% OP_J4U_J4V_TP: assemble the matrix A = [a(i,j)], a(i,j) = (J4 (u_j), J4 (v_i)), exploiting the tensor product structure.
%
%   mat = op_j4u_j4v_tp (space_u, space_v, msh, a, Ef);
%   [rows, cols, values] = op_j4u_j4v_tp (space_u, space_v, msh, a, Ef);
%
% INPUT:
%    
%   space_u:     object representing the space of trial functions (see sp_vector)
%   space_v:     object representing the space of test functions (see sp_vector)
%   msh:         object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   a:           constant vector of length = ndir representing the
%                direction of the fiber in the material
%   Ef:          constant representing strength of the fiber
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
function varargout = op_j4u_j4v_tp(space_u, space_v, msh, a, Ef)
  A = spalloc (space_v.ndof, space_u.ndof, 5*space_u.ndof);
  a_kron_a = kron(a',a);
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space_u, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);
    sp2_col = sp_evaluate_col (space_v, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    
    A = A + Ef*op_j4u_j4v(sp1_col, sp2_col, msh_col, a_kron_a);
  end

  if (nargout == 1)
    varargout{1} = A;
    
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end
end