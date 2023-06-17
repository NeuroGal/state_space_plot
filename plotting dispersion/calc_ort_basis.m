function ort_mat = calc_ort_basis(vec)
% vec is the first in the basis, calculates the full orthonormal matrix (as
% columns)
normalize_vec = @(vec) vec/norm(vec);

if length(vec) == 3
    ort_vec1 = [0 vec(3) -vec(2)];ort_vec1 = normalize_vec(ort_vec1);
    ort_vec2 = cross(vec,ort_vec1);ort_vec2 = normalize_vec(ort_vec2);
    ort_mat = [vec; ort_vec1; ort_vec2]';
else
    ort_vec = [vec(2) -vec(1)];ort_vec = normalize_vec(ort_vec);
    ort_mat = [vec; ort_vec]';
end
end