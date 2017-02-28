
% Distance metric - log-Euclidian method
function d = covariance_distance_log_euclidian(C1,C2)
    % perform SVD
    [U1,D1,V1] = svd(C1);
    [U2,D2,V2] = svd(C2);
    
    %take log of D matrix
    D1_log = log(diag(D1));
    D2_log = log(diag(D2));


    D1_log_matrix = diag(D1_log);
    D2_log_matrix = diag(D2_log);
    
    C1_log = U1*D1_log_matrix*V1';
    C2_log = U2*D2_log_matrix*V2';
    
    % calculate distacne
    diff = abs(C1_log -C2_log).^2;
    d = sum(sum(diff));
    
end