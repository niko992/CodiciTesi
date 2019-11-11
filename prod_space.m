function M = prod_space(sp1,sp2)
size1 = size(sp1);
size2 = size(sp2);
for i = 1 : size1(end)
    for j = 1 : size2(end)
        M(i,j) = sum(sum (sp1(:,:,i) .* sp2(:,:,j),1),2);
    end
end

end