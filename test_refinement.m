function [eu,X,Y] = test_refinement(u, space, geometry)
    vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
    [eu, F] = sp_eval (u, space, geometry, vtk_pts);
    [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
end