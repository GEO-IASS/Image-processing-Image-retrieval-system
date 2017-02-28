
% Covariance Matrix method for image base retrival system using
% log-Euclidian distace
% Extracts a 9 dimensional feature vector
% Performs image integrals for covariance matrix calculations
% Calculates distance between images using generalized eigenvalues
% Calculates confusion matrices.

% Silvia G. Ionescu
% 12-15-2016

clear; close all; clc;

% import the picture folders or folder 
label = {'baseball', 'beach', 'bridge', 'cars', 'documents', 'family', 'food', 'landscape', 'office', 'pets', 'wedding'};

for i = 1:length(label)
    disp(i)
    cd(['/ad/eng/users/s/i/sionescu/Documents/MATLAB/cov_test_305/', char(label(i))])
    imagefiles = dir('*.JPEG');
    imagefiles = imagefiles(arrayfun(@(x) ~strcmp(x.name(1),'.'),imagefiles));
    for j = 1: length(imagefiles)
        disp(j)
        tic
        image_name = imagefiles(j).name;
        I = imread(image_name);
        
        if size(I,3) == 3
            
            I = imresize(I, [500 500]);
            unorganized_images(j) = {I};
            % calculates the feature vectors for each image
            F = feature_vector(I); 
            % calculates the covariance matrices for each image
            cov = covariance_region(F);
            cov_unorganized(j,i) = {cov};
            folder_label(j,i) = i;
            
         toc
        end
    end
end

cov_matrix = cov_unorganized(:);
image_label = folder_label(:);

empty = find(image_label == 0); 
diff = setdiff(1:length(cov_matrix),empty);
cov_matrix = cov_matrix(diff);
image_label = image_label(diff); 

% calculates the distance between each test image and the rest of the
% database images
for i = 1:length(cov_matrix)  
    disp(i)
    test_cov = cov_matrix{i};  
    for j = 1:length(cov_matrix)
        folder_sample = cov_matrix{j};
        for k = 1:49
            distance(k,j) = covariance_distance_log_euclidian(test_cov{1,k}, folder_sample{1,k});              
        end
    end
    dis = sum(distance,1);
    max_distance = max(distance,[],1); 
    final_distance = dis-max_distance; 
    [sort_dis,order] = sort(final_distance);
    predicted_label(:,i) = image_label(order(2:6));
end

% calculate confusion matices for the 5 closest images
confusion_matrix_d1 = confusionmat(image_label,predicted_label(1,:));
confusion_matrix_d2 = confusionmat(image_label,predicted_label(2,:));
confusion_matrix_d3 = confusionmat(image_label,predicted_label(3,:));
confusion_matrix_d4 = confusionmat(image_label,predicted_label(4,:));
confusion_matrix_d5 = confusionmat(image_label,predicted_label(5,:));


for i=1:11
    ccr1(i)=(confusion_matrix_d1(i,i))./(sum(confusion_matrix_d1(i,:)));
end