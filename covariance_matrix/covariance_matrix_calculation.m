
% Calculate the area covariance matrix
function C = covariance_matrix_calculation(Q, P, d, region_size)

P_repmat = repmat(P, [1,d]);
P_rep_perm = permute(P_repmat, [2,1]);
P_mult = P_repmat.*P_rep_perm;

C = (1/(region_size-1))*(double(Q) - (1/(region_size))* double(P_mult));

end