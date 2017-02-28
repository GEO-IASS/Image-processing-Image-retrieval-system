
% Calculate the feature vector for (x,y) position, R,G,B and first and 
% second order derivatives

function F = feature_vector(I)
% I must be a threedimensional color image
I_gray = rgb2gray(I);
I = double(I);
I_gray = double(I_gray);

[n,m,~] = size(I);

F = zeros(9,n,m);
x = 1:size(I,1);
y = 1:size(I,2);
[X,Y] = meshgrid(x,y);

F(1,:,:) = X';
F(2,:,:) = Y';
F(3:5,:,:) = permute(I,[3 1 2]);



% [FX,FY] = gradient(I_gray);

h_x = [-1 0 1];
h_y = h_x';

FX = imfilter(I_gray,h_x, 'symmetric');
FY = imfilter(I_gray,h_y, 'symmetric');



% hh_x = 1/sqrt(6)*[-1 0 1; -1 0 1; -1 0 1];
% hh_y = hh_x';
% FX = imfilter(I_gray,hh_x, 'symmetric');
% FY = imfilter(I_gray,hh_y, 'symmetric');



F(6,:,:) = FX;
F(7,:,:) = FY;

% h_xx = [-1 2 -1];
% h_yy = h_xx';
% 
% FXX = imfilter(I_gray,h_xx, 'symmetric');
% FYY = imfilter(I_gray,h_yy, 'symmetric');
% figure(2)
% subplot(1,2,1)
% imshow(FXX);
% 
% subplot(1,2,2)
% imshow(FYY);
 
[FXX, FYY] = gradient(gradient(I_gray));

F(8,:,:) = FXX;
F(9,:,:) = FYY;

F = int64(F);
end