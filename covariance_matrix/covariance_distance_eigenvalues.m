
% Distance metric - generalized eigenvalues method

function d = covariance_distance_eigenvalues(C1,C2)

min = 1e-6; 

% calculate the generalized eigenvalues
[~,D] = eig(C1 + diag(min*ones(9,1)),C2 + diag(min*ones(9,1)));

if(chol(A)~= 0)
    disp(' Matrix is not positive definite symetric');
    return; 
end

d = sum(log(diag(abs(D))).^2);
end