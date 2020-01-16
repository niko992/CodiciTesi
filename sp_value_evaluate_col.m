function sp = sp_value_evaluate_col(space, msh_col, msh)

    if space.degree == [0,0]
        sp_val = sp_evaluate_col (space, msh_col, 'value', true);
        base = sp_val.shape_functions;
        sp.space_functions = base;
        sp.dim = 1;
        sp.connectivity = connectivity_matrix(msh_col, msh);
        sp.ndof = msh.nel_dir(1)*msh.nel_dir(2)*sp.dim;
    end
end