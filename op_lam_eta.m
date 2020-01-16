function varargout = op_lam_eta (sp_lam, sp_eta, msh)

  lam = sp_lam.space_functions;
  eta = sp_eta.space_functions;
  ndir = msh.rdim;
  rows = zeros (msh.nel * sp_lam.dim * sp_eta.dim, 1);
  cols = zeros (msh.nel * sp_lam.dim * sp_eta.dim, 1);
  values = zeros (msh.nel * sp_lam.dim * sp_eta.dim, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))

      lam_iel = reshape (lam(:,:,:,1:sp_lam.dim,iel), msh.ndim, ndir, msh.nqn, sp_lam.dim);
      eta_iel = reshape (eta(:,:,:,1:sp_eta.dim,iel), msh.ndim, ndir, msh.nqn, sp_eta.dim);
      eta_iel = reshape (lam_iel, [msh.ndim*ndir, msh.nqn, 1, sp_eta.dim]);
      lam_iel = reshape (lam_iel, [msh.ndim*ndir, msh.nqn, sp_lam.dim, 1]);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

      jacdet_lam = bsxfun (@times, jacdet_iel, lam_iel);
      aux_val1 = sum(bsxfun(@times, jacdet_lam, eta_iel),1);

      values(ncounter+(1:(sp_lam.dim*sp_eta.dim))) = reshape(sum(aux_val1,2),sp_lam.dim,sp_eta.dim);

      [rows_loc, cols_loc] = ndgrid (sp_lam.connectivity(:,iel), sp_eta.connectivity(:,iel));
      rows(ncounter+(1:(sp_lam.dim*sp_eta.dim))) = rows_loc;
      cols(ncounter+(1:(sp_lam.dim*sp_eta.dim))) = cols_loc;
      ncounter = ncounter + sp_lam.dim*sp_eta.dim;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_lam_eta: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), sp_lam.ndof, sp_eta.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_lam_eta: wrong number of output arguments')
  end
end
