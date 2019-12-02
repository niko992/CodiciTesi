function varargout = op_j4lam_j4eta_tp(space_lam, space_eta, msh, Ef, a, ndof_lam)
  A = spalloc (ndof_lam, ndof_lam, 5*ndof_lam);
  a_kron_a = kron(a',a);
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space_lam, msh_col, 'value', true);
    sp2_col = sp_evaluate_col (space_eta, msh_col, 'value', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    lam_connectivity = create_connectivity_matrix(length(msh_col.elem_list), iel);
    A = A + Ef*op_j4lam_j4eta(sp1_col, sp2_col, msh_col, a_kron_a, lam_connectivity, ndof_lam);
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