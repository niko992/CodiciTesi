% OP_SU_EV_TP: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)), exploiting the tensor product structure.
%
%   mat = op_eu_ev_tp (spu, spv, msh, lambda, mu);
%   [rows, cols, values] = op_eu_ev_tp (spu, spv, msh, lambda, mu);
%
% INPUT:
%    
%   spu:     object representing the space of trial functions (see sp_vector)
%   spv:     object representing the space of test functions (see sp_vector)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   lambda, mu: function handles to compute the Lame' coefficients
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 


function varargout = op_eu_ev_tp (space1, space2, msh, mu)

  A = spalloc (space2.ndof, space1.ndof, 5*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', true, ...
                               'gradient', true, 'divergence', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A = A + op_eu_ev (sp1_col, sp2_col, msh_col, mu (x{:}));
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
