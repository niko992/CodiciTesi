function u_fine = project_into_finer_space(c_space, f_space, c_msh,f_msh,c_geometry,f_geometry, u_coarse)
msh=msh_precompute(c_msh);
msh_ref = msh_precompute(f_msh);
space=sp_precompute(c_space,msh);
space_ref = sp_precompute(f_space, msh_ref);
px = cell2mat(f_geometry.nurbs.knots(1));
px = px(2:end-1);
py = cell2mat(f_geometry.nurbs.knots(2));
py = py(2:end-1);
u_reshaped = reshape(u_coarse, space.ncomp, space.ndof_dir(1,1), space.ndof_dir(2,2));
geom = nrbmak(u_reshaped,c_geometry.nurbs.knots);
geom_f = nrbkntins(geom, {setdiff(f_geometry.nurbs.knots{1},c_geometry.nurbs.knots{1}), setdiff(f_geometry.nurbs.knots{2},c_geometry.nurbs.knots{2})});
u_fine = geom_f.coefs(1:2,:,:);
u_fine = reshape(u_fine,geom_f.number(1)*geom_f.number(2)*2,1);
end