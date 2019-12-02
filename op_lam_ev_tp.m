function varargout = op_lam_ev_tp (space1, space2, msh, ndof_lam)
A = spalloc (space2.ndof, msh.nel_dir(1)*msh.nel_dir(2)*3, 5*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', true, ...
                               'gradient', false);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', true, ...
                               'gradient', true, 'divergence', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    lam_connectivity = create_connectivity_matrix(length(msh_col.elem_list), iel);
    A = A + op_lam_ev (sp1_col, sp2_col, msh_col, lam_connectivity, ndof_lam);
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