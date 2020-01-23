function varargout = op_p_divv (spp, spv, msh)
%gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  p = reshape(spp.shape_functions, spp.ncomp, msh.nqn, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
  ndir = msh.rdim

  rows = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spp.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
%       gradu_iel = reshape (gradu(:,:,:,1:spu.nsh(iel),iel), spu.ncomp, ndir, msh.nqn, spu.nsh(iel));
%       epsu_iel = (gradu_iel + permute (gradu_iel, [2 1 3 4]))/2;
%       epsu_iel = reshape (epsu_iel, [spu.ncomp*ndir, msh.nqn, 1, spu.nsh(iel)]);
%       epsu_iel = repmat (epsu_iel, [1,1,spv.nsh(iel),1]);

      p_iel = reshape(p(:,:,iel), spp.ncomp, msh.nqn, spp.nsh(iel));
      
      gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp, ndir, msh.nqn, spv.nsh(iel));
      epsv_iel = (gradv_iel + permute (gradv_iel, [2 1 3 4]))/2;
      epsv_iel = reshape (epsv_iel, [spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1]);
%       epsv_iel = repmat (epsv_iel, [1,1,1,spu.nsh(iel)]);

%       divu_iel = reshape (spu.shape_function_divs(:,1:spu.nsh(iel),iel), [msh.nqn, 1, spu.nsh(iel)]);
%       divu_iel = repmat (divu_iel, [1, spv.nsh(iel), 1]);
      divv_iel = reshape (spv.shape_function_divs(:,1:spv.nsh(iel),iel), [msh.nqn, spv.nsh(iel), 1]);
      divv_iel = repmat (divv_iel, [1, 1, spv.nsh(iel)]);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
%       jacdet_lambda_iel = reshape (jacdet_weights(:,iel), [msh.nqn,1,1]);

%       aux_val1 = 2 * sum (bsxfun (@times, jacdet_epsu, epsv_iel), 1);
      aux_val1 = bsxfun (@times, jacdet_iel, p_iel .* divv_iel);
      values(ncounter+(1:spp.nsh(iel)*spv.nsh(iel))) = reshape (sum (aux_val1, 2), spv.nsh(iel), spp.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spp.connectivity(:,iel));
      rows(ncounter+(1:spp.nsh(iel)*spv.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spp.nsh(iel)*spv.nsh(iel))) = cols_loc;
      ncounter = ncounter + spp.nsh(iel)*spv.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spp.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_eu_ev: wrong number of output arguments')
  end

end

%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% 
% function mat = op_su_ev (spu, spv, msh, lambda, mu)
%   
%   mat = spalloc (spv.ndof, spu.ndof, 1);
%   
%   gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
%   gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
% 
%   ndir = size (gradu, 2);
% 
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%       mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
%       for idof = 1:spv.nsh(iel)
%         ishg  = gradv(:,:,:,idof,iel);
%         ishgt = permute (ishg, [2, 1, 3]);
%         ieps  = reshape(ishg + ishgt, spv.ncomp * ndir, [])/2;
%         idiv  = spv.shape_function_divs(:, idof, iel);
%         for jdof = 1:spu.nsh(iel) 
%           jshg  = gradu(:,:,:,jdof,iel);
%           jshgt = permute (jshg, [2, 1, 3]);
%           jeps  = reshape(jshg + jshgt, spu.ncomp * ndir, [])/2;
%           jdiv  = spu.shape_function_divs(:, jdof, iel);
%  % The cycle on the quadrature points is vectorized         
%           mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%               sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
%                    (2 * sum (ieps .* jeps, 1).' .* mu(:,iel)  + ...
%                     (idiv .* jdiv) .* lambda(:,iel)));
%         end
%       end
%       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
%         mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
%     end
%   end
% 
% end
