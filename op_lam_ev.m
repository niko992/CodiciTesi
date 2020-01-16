function varargout = op_lam_ev (sp_lam, sp_u, msh)

  grad_u = reshape (sp_u.shape_function_gradients, sp_u.ncomp, [], msh.nqn, sp_u.nsh_max, msh.nel);
  lam = sp_lam.space_functions;
  ndir = size (grad_u, 2);
  rows = zeros (msh.nel * sp_lam.dim * sp_u.nsh_max, 1);
  cols = zeros (msh.nel * sp_lam.dim * sp_u.nsh_max, 1);
  values = zeros (msh.nel * sp_lam.dim * sp_u.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      lam_iel = reshape (lam(:,:,:,1:sp_lam.dim,iel), sp_u.ncomp, ndir, msh.nqn, sp_lam.dim);
      lam_iel = reshape (lam_iel, [sp_u.ncomp*ndir, msh.nqn, 1, sp_lam.dim]);

      
      grad_u_iel = reshape (grad_u(:,:,:,1:sp_u.nsh(iel),iel), sp_u.ncomp, ndir, msh.nqn, sp_u.nsh(iel));
      eps_u_iel = (grad_u_iel + permute (grad_u_iel, [2 1 3 4]))/2;
      eps_u_iel = reshape (eps_u_iel, [sp_u.ncomp*ndir, msh.nqn, sp_u.nsh(iel), 1]);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

      jacdet_u = bsxfun (@times, jacdet_iel, eps_u_iel);
      aux_val1 = sum(bsxfun(@times,jacdet_u,lam_iel),1);
      values(ncounter+(1:(sp_u.nsh(iel)*sp_lam.dim))) = reshape(sum(aux_val1,2),sp_u.nsh(iel),sp_lam.dim);

      [rows_loc, cols_loc] = ndgrid (sp_u.connectivity(:,iel), sp_lam.connectivity(:,iel));
      rows(ncounter+(1:(sp_lam.dim*sp_u.nsh(iel)))) = rows_loc;
      cols(ncounter+(1:(sp_lam.dim*sp_u.nsh(iel)))) = cols_loc;
      ncounter = ncounter + sp_lam.dim*sp_u.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_p_ev: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), sp_u.ndof, sp_lam.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_p_ev: wrong number of output arguments')
  end

end