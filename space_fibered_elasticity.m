function [sp_u, sp_lam] = space_fibered_elasticity (knots, nsub_u, degree_u, regularity_u, msh)
% knots_lam = kntrefine (knots, nsub_u-1, [1, 1], regularity_u);
% sp_lam = sp_bspline (knots_lam,[1, 1],msh);
sp_lam = sp_bspline(knots,[0, 0],msh);
knots_u = kntrefine (knots,nsub_u-1,degree_u, regularity_u);
scalar_space = sp_bspline (knots_u, degree_u, msh);
    for idim = 1:msh.ndim
      scalar_spaces{idim} = scalar_space;
    end
sp_u = sp_vector (scalar_spaces, msh);
clear scalar_spaces scalar_space
end