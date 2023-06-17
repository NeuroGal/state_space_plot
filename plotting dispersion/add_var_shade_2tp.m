function add_var_shade_2tp(pt1, pt2, r, col)
if numel(r) == 1; r = repmat(r,2,1);end
vec = pt2-pt1;
if length(pt1) == 3
    theta = (0:0.1:2*pi)';L = length(theta); 
    gen_circ = @(center, r, theta, dir1, dir2) center + r.*cos(theta).*dir1 + r.*sin(theta).*dir2;
    ort_mat = calc_ort_basis(vec)';
    ort_vec1 = ort_mat(2,:);ort_vec2 = ort_mat(3,:);
    pt1_var = gen_circ(pt1, r(1), theta, ort_vec1, ort_vec2); 
    pt2_var = gen_circ(pt2, r(2), theta, ort_vec1, ort_vec2); 
    f = (1:L)'; F = [f, circshift(f,L-1), circshift(f+L,L-1),f+L];
else
    ort_mat = calc_ort_basis(vec)'; ort_vec =  ort_mat(2,:);
    gen_var_line = @(center, r, dir) center + r.*[dir; -dir];
    pt1_var = gen_var_line(pt1, r(1), ort_vec);
    pt2_var = gen_var_line(pt2, r(2), ort_vec);
    F = [1 2 4 3];
end
patch('Faces',F,'Vertices',[pt1_var;pt2_var],'FaceAlpha',0.4,'EdgeColor' ,'none','FaceColor',col);
end