function varargout = op_lam_eta (sp_lam, sp_eta, msh, lam_connectivity)

  base = sp_lam.shape_functions;
  add = zeros(size(base,1),1,msh.nel);
  lam11 = cat(1,base,add,add);
  lam12 = cat(1,add,base,add);
  lam22 = cat(1,add,add,base);
  lam = cat(2,lam11,lam12,lam12,lam22);
  lam = permute(lam,[2,1,3]);
  lam = reshape(lam,msh.rdim,msh.ndim,msh.nqn,msh.ndim*(msh.ndim+1)/2,msh.nel);
  ndir = msh.rdim;
% sp_lam.nsh = 3???
  rows = zeros (msh.nel * 3 * 3, 1);
  cols = zeros (msh.nel * 3 * 3, 1);
  values = zeros (msh.nel * 3 * 3, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
%       p_iel = reshape (p(:,:,:,1:spp.nsh(iel),iel), spp.ncomp, ndir, msh.nqn, spp.nsh(iel));
%       epsu_iel = repmat (epsu_iel, [1,1,sp_u.nsh(iel),1]);
      lam_iel = reshape (lam(:,:,:,1:3,iel), msh.ndim, ndir, msh.nqn, 3);
      lam_iel = reshape (lam_iel, [msh.ndim*ndir, msh.nqn, 3, 1]);

      
%       grad_u_iel = reshape (grad_u(:,:,:,1:sp_u.nsh(iel),iel), sp_u.ncomp, ndir, msh.nqn, sp_u.nsh(iel));
%       eps_u_iel = (grad_u_iel + permute (grad_u_iel, [2 1 3 4]))/2;
%       eps_u_iel = reshape (eps_u_iel, [sp_u.ncomp*ndir, msh.nqn, sp_u.nsh(iel), 1]);
%       eps_u_iel = repmat (eps_u_iel, [1,1,1,spu.nsh(iel)]);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

      jacdet_lam = bsxfun (@times, jacdet_iel, lam_iel);
      aux_val1 = prod_space(jacdet_lam,lam_iel);

      values(ncounter+(1:(3*3))) = reshape(aux_val1,3,3);

      [rows_loc, cols_loc] = ndgrid (lam_connectivity(:,iel), lam_connectivity(:,iel));
      test = 10;
      rows(ncounter+(1:(3*3))) = rows_loc;
      cols(ncounter+(1:(3*3))) = cols_loc;
      ncounter = ncounter + 3*3;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_p_ev: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), 12, 12);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_p_ev: wrong number of output arguments')
  end
end
