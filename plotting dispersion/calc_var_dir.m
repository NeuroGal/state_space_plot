function std_vec = calc_var_dir(pts_standard_basis, vec)
% pts_standard_basis should be n_pts x 2\3
% calculate the std of the points in the direction(s) orthogonal to vec 
ort_mat = calc_ort_basis(vec);
new_basis_rep = ort_mat \ pts_standard_basis';
std_vec = std(sum(new_basis_rep(2:end,:),1));
end