
% Implementation of the image integrals to compute the covaraince matrix
% for different region of the image. 
function covariance = covariance_region(F)

[d,n,m] = size(F);

% P is a W × H × d tensor of the integral images
P = int64(zeros(d,n+1,m+1));
for i=1:d
    P(i,:,:) = integralImage(squeeze(int64(F(i,:,:))));
end
P = P(:,2:end,2:end);

% Q is a W × H × d × d tensor of the second order integral images
F2 = int64(zeros(d,d,n,m));
for i=1:d
    for j=1:d
        F2(i,j,:,:) = F(i,:,:).*F(j,:,:);
    end
end
Q = int64(zeros(d,d,n+1,m+1));
for i=1:d
    for j=1:d
        Q(i,j,:,:) = integralImage(squeeze(F2(i,j,:,:)));
    end
end
Q = Q(:,:,2:end,2:end);

% calculate the covariance matrix for the entire image
Q1 = Q(:,:,n,m);
P1 = P(:,n,m);
P1_repmat = repmat(P1, [1,d]);
P1_rep_perm = permute(P1_repmat, [2,1]);
P1_mult = P1_repmat.*P1_rep_perm;
C1 = (1/(size(F,2)*size(F,3)-1))*(double(Q1) - (1/(size(F,2)*size(F,3)))* double(P1_mult));

% calculate the covariance matrix for the first m/4 x m/4 region 
Q1 = Q(:,:,m/4,n/4);
P1 = P(:,n/4,m/4);
P1_repmat = repmat(P1, [1,d]);
P1_rep_perm = permute(P1_repmat, [2,1]);
P1_mult = P1_repmat.*P1_rep_perm;
region_size = (m/4)*(n/4);
C1 = 1/(region_size-1)*(double(Q1) - (1/(region_size)* double(P1_mult)));
covariance(1) = {C1};


k = 0;
% calculate different regions covariance matrix
for i = m/2:m/4:m
    for j = n/2:n/4:n
        k = k + 1;
        Q1 = Q(:,:,(i-m/4),(j-n/4));
        P1 = P(:,(i-m/4),(j-n/4));
        region_size1 = (i-m/4)*(j-n/4);
        if (i-m/4)== 1 || (j-n/4) ==1
            C1 = zeros(9,9);
            Q1 = zeros(9,9);
            P1 = zeros(9,1);
        else
            C1 = covariance_matrix_calculation(Q1,P1,d,region_size1);
            covariance(k) = {C1};
            k = k+1;
        end
        
        Q2 = Q(:,:,i,(j-n/4));
        P2 = P(:,i,(j-n/4) );
        region_size2 = i * (j-n/4);
        if i == 1 || (j-n/4+1) == 1
            C2 = zeros(9,9);
            Q2 = zeros(9,9);
            P2 = zeros(9,1);
        else
            C2 = covariance_matrix_calculation(Q2,P2,d,region_size2);
            covariance(k) = {C2};
            k = k+1;          
        end
        
        Q3 = Q(:,:,(i-m/4),j);
        P3 = P(:,(i-m/4),j);
        region_size3 = (i-m/4) * j;
        if (i-m/4) == 1 || j ==1
            C3 = zeros(9,9);
            Q3 = zeros(9,9);
            P3 = zeros(9,1);
        else
            C3 = covariance_matrix_calculation(Q3,P3,d,region_size3);
            covariance(k) = {C3};
            k = k+1;
        end
        
        Q4 = Q(:,:,i,j);
        P4 = P(:,i,j);
        region_size4 = i * j;
        C4 = covariance_matrix_calculation(Q4,P4,d,region_size4);
        covariance(k) = {C4};
        k = k+1;
        
        Q5 = Q4 + int64(Q1) - int64(Q2) - int64(Q3);
        P5 = P4 + int64(P1) - int64(P2) - int64(P3);
        region_size5 = (i * j) - ((i-m/4+1)*(j-n/4+1)) + 1;
        C5 = covariance_matrix_calculation(Q5,P5,d,region_size5);
        if sum(sum(C4)) ~= sum(sum(C5))
            covariance(k) = {C5};
        end
    end
end
k = k+1;
% % calculate image left half covariance
% Q2 = Q(:,:,m,n/2);  
% P2 = P(:,m,n/2);
% P2_repmat = repmat(P2, [1,d]);
% P2_rep_perm = permute(P2_repmat, [2,1]);
% P2_mult = P2_repmat.*P2_rep_perm;
% C2 = (1/((m/2*n)-1))*(double(Q2) - (1/(m/2*n))* double(P2_mult));
% covariance(k) = {C2};
% k = k+1;

% % calculate image upper half covariance matrix
% Q3 = Q(:,:,m/2, n);  
% P3 = P(:,m/2, n);
% P3_repmat = repmat(P3, [1,d]);
% P3_rep_perm = permute(P3_repmat, [2,1]);
% P3_mult = P3_repmat.*P3_rep_perm;
% C3 = (1/((m/2*n)-1))*(double(Q3) - (1/(m/2*n))* double(P3_mult));
% covariance(k) = {C3};
% k = k+1;

% calculate image right half covariance

Q4 = Q(:,:, m,n) - Q2;
P4 = P(:,m,n) - P2;
P4_repmat = repmat(P4, [1,d]);
P4_rep_perm = permute(P4_repmat, [2,1]);
P4_mult = P4_repmat.*P4_rep_perm;
C4 = (1/((m/2*n)-1))*(double(Q4) - (1/(m/2*n))* double(P4_mult));
covariance(k) = {C4};
k = k+1;

% calculate image bottom half covariance matrix
Q5 = Q(:,:, m, n) - Q3;
P5 = P(:, m, n) - P3;
P5_repmat = repmat(P5, [1,d]);
P5_rep_perm = permute(P5_repmat, [2,1]);
P5_mult = P5_repmat.*P5_rep_perm;
C5 = (1/((m/2*n)-1))*(double(Q5) - (1/(m/2*n))* double(P5_mult));
covariance(k) = {C5};
k = k+1;

end