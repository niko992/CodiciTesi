% OP_P_DIVV_TP: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)), exploiting the tensor product structure.
%
%   mat = op_su_ev_tp (spu, spv, msh, lambda, mu);
%   [rows, cols, values] = op_su_ev_tp (spu, spv, msh, lambda, mu);
%
% INPUT:
%    
%   spp:     object representing the space of trial functions (see sp_vector)
%   spv:     object representing the space of test functions (see sp_vector)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2011 Rafael Vazquez
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

function varargout = op_p_divv_tp (space1, space2, msh)

  A = spalloc (space2.ndof, space1.ndof, 5*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', true, ...
                               'gradient', false);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A = A + op_p_divv (sp1_col, sp2_col, msh_col);
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
