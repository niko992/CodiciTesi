function varargout = op_j4u_j4v(spu, spv, msh, a_kron_a)
%   gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
%   gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
  gradu = spu.shape_function_gradients;
  gradv = spv.shape_function_gradients;
  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
%   jacdet_weights_mu = jacdet_weights .* mu;
%   jacdet_weights_lambda = jacdet_weights .* lambda;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      gradu_iel = reshape (gradu(:,:,:,1:spu.nsh(iel),iel), spu.ncomp, ndir, msh.nqn, spu.nsh(iel));
      epsu_iel = (gradu_iel + permute (gradu_iel, [2 1 3 4]))/2;
      epsu_iel = reshape (epsu_iel, [spu.ncomp*ndir, msh.nqn, 1, spu.nsh(iel)]);
%       epsu_iel = repmat (epsu_iel, [1,1,spv.nsh(iel),1]);
      
      gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp, ndir, msh.nqn, spv.nsh(iel));
      epsv_iel = (gradv_iel + permute (gradv_iel, [2 1 3 4]))/2;
      epsv_iel = reshape (epsv_iel, [spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1]);
%       epsv_iel = repmat (epsv_iel, [1,1,1,spu.nsh(iel)]);

      a_kron_a = reshape(a_kron_a, [1,ndir * ndir])';
      vec = ones(1,msh.nqn);
      
      a_kron_a_iel = kron(vec,a_kron_a);
      
      
      a_kron_a_epsu = bsxfun (@times, a_kron_a_iel, epsu_iel);
%       a_kron_a_epsu = reshape (a_kron_a_epsu, [spu.ncomp*ndir, msh.nqn, 1, spu.nsh(iel)]);

      a_kron_a_epsv = bsxfun (@times, a_kron_a_iel, epsv_iel);
%       a_kron_a_epsv = reshape (a_kron_a_epsv, [spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1]);
      
      aux_val = sum (bsxfun (@times, a_kron_a_epsu, a_kron_a_epsv), 1);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (aux_val, 2), spv.nsh(iel), spu.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc;
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
  end
  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_su_ev: wrong number of output arguments')
  end
end