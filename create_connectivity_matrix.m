function M = create_connectivity_matrix(nsub, iel)
    for i = 1:nsub
        M(:,i)= (iel + (i-1)*nsub) * ones (3,1) + [0:(nsub*nsub):(2*nsub*nsub)]';
    end
end