function varargout = op_j4lam_j4eta(sp_lam, sp_eta, msh, a_kron_a, lam_connectivity, ndof_lam)
  base = sp_lam.shape_functions;
  add = zeros(size(base,1),1,msh.nel);
  lam11 = cat(1,base,add,add);
  lam12 = cat(1,add,base,add);
  lam22 = cat(1,add,add,base);
  lam = cat(2,lam11,lam12,lam12,lam22);
  lam = permute(lam,[2,1,3]);
  lam_dim = msh.ndim*(msh.ndim+1)/2;
  lam = reshape(lam,msh.rdim,msh.ndim,msh.nqn,lam_dim,msh.nel);
  ndir = msh.rdim;

  rows = zeros (msh.nel * lam_dim * lam_dim, 1);
  cols = zeros (msh.nel * lam_dim * lam_dim, 1);
  values = zeros (msh.nel * lam_dim * lam_dim, 1);
  
  jacdet_weights = msh.jacdet .* msh.quad_weights;

  a_kron_a = reshape(a_kron_a, [ndir * ndir,1]);
      
  a_kron_a_rep = reshape(repmat(a_kron_a,[msh.nqn,msh.nel]),[ndir*ndir,msh.nqn,msh.nel]);

%   jacdet_weights_rep = permute(repmat(jacdet_weights,1,1,ndir*ndir),[3,1,2]);
%   a_kron_a_weights = bsxfun(@times,a_kron_a_rep,jacdet_weights_rep);
    
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))

      lam_iel = reshape (lam(:,:,:,1:lam_dim,iel), msh.ndim, ndir, msh.nqn, lam_dim);
      lam_iel = reshape (lam_iel, [msh.ndim*ndir, msh.nqn, lam_dim, 1]);
      eta_iel = lam_iel;
      
%       a_kron_a_weights_iel = a_kron_a_weights(:,:,iel);
      a_kron_a_iel = a_kron_a_rep(:,:,iel);
      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
      a_kron_a_weights_iel = bsxfun (@times, jacdet_iel, a_kron_a_iel);

%       a_kron_a_weight_epsv = sum(bsxfun (@times, epsv_iel,a_kron_a_weights_iel),1);
      a_kron_a_weights_e = sum(bsxfun (@times, lam_iel, a_kron_a_weights_iel),1);
      a_kron_a_eta = sum(bsxfun(@times, eta_iel, a_kron_a_iel));
      a_kron_a_eta = reshape(a_kron_a_eta,size(a_kron_a_eta,1),size(a_kron_a_eta,2),1,size(a_kron_a_eta,3));
      
      aux_val = sum(bsxfun(@times,a_kron_a_weights_e, a_kron_a_eta),1);
      values(ncounter+(1:(lam_dim*lam_dim))) = reshape (sum(aux_val,2), lam_dim, lam_dim);

      [rows_loc, cols_loc] = ndgrid (lam_connectivity(:,iel), lam_connectivity(:,iel));
      rows(ncounter+(1:(lam_dim*lam_dim))) = rows_loc;
      cols(ncounter+(1:(lam_dim*lam_dim))) = cols_loc;
      ncounter = ncounter + lam_dim*lam_dim;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_lam_ev: singular map in element number %d', iel)
    end
  end
  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), ndof_lam, ndof_lam);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_su_ev: wrong number of output arguments')
  end
end