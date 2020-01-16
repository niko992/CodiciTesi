function errl2 = sp_l2_error_refined (space, space_ref, msh_ref, u, u_ref)

space = space.constructor(msh_ref);

errl2 = 0;

for iel = 1:msh_ref.nel_dir(1)
    msh_ref_col = msh_evaluate_col (msh_ref, iel);
    sp_ref_col = sp_evaluate_col (space_ref, msh_ref_col, 'value', true, 'gradient', false);
    sp_col = sp_evaluate_col (space, msh_ref_col, 'value', true, 'gradient', false);

    valu = sp_eval_msh (u, sp_col, msh_ref_col);
    valu = reshape (valu, sp_col.ncomp, msh_ref_col.nqn, msh_ref_col.nel);

    valu_ref = sp_eval_msh (u_ref, sp_ref_col, msh_ref_col);
    valu_ref = reshape (valu_ref, sp_col.ncomp, msh_ref_col.nqn, msh_ref_col.nel);

    w_ref = msh_ref_col.quad_weights .* msh_ref_col.jacdet;

    errl2 = errl2 + sum(sum (reshape (sum ((valu - valu_ref).^2, 1), [msh_ref_col.nqn, msh_ref_col.nel]) .* w_ref));
end

errl2 = sqrt (errl2);

end