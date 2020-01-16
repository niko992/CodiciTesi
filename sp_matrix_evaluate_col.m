function sp = sp_matrix_evaluate_col(space, msh_col, msh)

    if space.degree == [0,0]
        sp_lam = sp_evaluate_col (space, msh_col, 'value', true);
        base = sp_lam.shape_functions;
        add = zeros(size(base,1),1,msh_col.nel);
        lam11 = cat(1,base,add,add);
        lam12 = cat(1,add,base,add);
        lam22 = cat(1,add,add,base);
        lam = cat(2,lam11,lam12,lam12,lam22);
        lam = permute(lam,[2,1,3]);
        lam_dim = msh_col.ndim*(msh_col.ndim+1)/2;
        lam = reshape(lam,msh_col.rdim,msh_col.ndim,msh_col.nqn,lam_dim,msh_col.nel);
        sp.space_functions = lam;
        sp.dim = lam_dim;
        sp.connectivity = connectivity_matrix(msh_col, msh);
        sp.ndof = msh.nel_dir(1)*msh.nel_dir(2)*sp.dim;
    end
end