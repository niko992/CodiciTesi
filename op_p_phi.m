function varargout = op_p_phi (spp, spphi, msh)
  p = reshape(spp.shape_functions, spp.ncomp, msh.nqn, msh.nel);
  phi = reshape(spphi.shape_functions, spphi.ncomp, msh.nqn, msh.nel);
  ndir = msh.rdim;

  rows = zeros (msh.nel * spphi.nsh_max * spp.nsh_max, 1);
  cols = zeros (msh.nel * spphi.nsh_max * spp.nsh_max, 1);
  values = zeros (msh.nel * spphi.nsh_max * spp.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
%       gradu_iel = reshape (gradu(:,:,:,1:spu.nsh(iel),iel), spu.ncomp, ndir, msh.nqn, spu.nsh(iel));
%       epsu_iel = (gradu_iel + permute (gradu_iel, [2 1 3 4]))/2;
%       epsu_iel = reshape (epsu_iel, [spu.ncomp*ndir, msh.nqn, 1, spu.nsh(iel)]);
%       epsu_iel = repmat (epsu_iel, [1,1,spv.nsh(iel),1]);
%       shpq_iel = reshape (spq.shape_functions(:, 1:spq.nsh(iel), iel), msh.nqn, spq.nsh(iel), 1);

      p_iel = reshape(p(:,:, iel), spp.ncomp, msh.nqn, 1);
      phi_iel = reshape(phi(:,:,iel), spphi.ncomp, msh.nqn, 1);
      
%       gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp, ndir, msh.nqn, spv.nsh(iel));
%       epsv_iel = (gradv_iel + permute (gradv_iel, [2 1 3 4]))/2;
%       epsv_iel = reshape (epsv_iel, [spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1]);
% %       epsv_iel = repmat (epsv_iel, [1,1,1,spu.nsh(iel)]);
% 
%       divu_iel = reshape (spu.shape_function_divs(:,1:spu.nsh(iel),iel), [msh.nqn, 1, spu.nsh(iel)]);
%       divu_iel = repmat (divu_iel, [1, spv.nsh(iel), 1]);
%       divv_iel = reshape (spv.shape_function_divs(:,1:spv.nsh(iel),iel), [msh.nqn, spv.nsh(iel), 1]);
%       divv_iel = repmat (divv_iel, [1, 1, spu.nsh(iel)]);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
%       jacdet_lambda_iel = reshape (jacdet_weights_lambda(:,iel), [msh.nqn,1,1]);

%       aux_val1 = 2 * sum (bsxfun (@times, jacdet_epsu, epsv_iel), 1);
      aux_val1 = bsxfun (@times, jacdet_iel, p_iel .* phi_iel);
      values(ncounter+(1:spp.nsh(iel)*spphi.nsh(iel))) = reshape (sum (aux_val1, 2), spp.nsh(iel), spphi.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spphi.connectivity(:,iel), spp.connectivity(:,iel));
      rows(ncounter+(1:spp.nsh(iel)*spphi.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spp.nsh(iel)*spphi.nsh(iel))) = cols_loc;
      ncounter = ncounter + spp.nsh(iel)*spphi.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), msh.nel*msh.nel, msh.nel*msh.nel);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_eu_ev: wrong number of output arguments')
  end
end