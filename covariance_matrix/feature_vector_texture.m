
% Calculate the feature vector for (x,y) position, first and 
% second order derivatives
 
function F = feature_vector_texture(I)
% I must be a threedimensional color image
I_gray = rgb2gray(I);
I = double(I);
I_gray = double(I_gray);

[n,m,~] = size(I);

F = zeros(9,n,m);
x = 1:size(I,1);
y = 1:size(I,1);
[X,Y] = meshgrid(x,y);

F(1,:,:) = X;
F(2,:,:) = Y;
F(3,:,:) = I_gray;

[FX,FY] = gradient(I_gray);
F(4,:,:) = FX;
F(5,:,:) = FY;

[FXX, FYY] = gradient(gradient(I_gray));
F(6,:,:) = FXX;
F(7,:,:) = FYY;

F = int64(F);
end