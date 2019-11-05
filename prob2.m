%
%   ECEn 671 - Fall 2010
%   Professor Neal K. Bangerter
%   Homework #3
%
%   Problem #2
%

clear;
clc;
close all;

%   Load in the matrices needed for this problem
load('prob2.mat');

for i=1:20
    a(:,i) = reshape(A(:,:,i),65536,1);
end
x_comp = reshape(x,1,65536);

%compute cross correlation vector
p = zeros(20,1);
for i = 1:20
    p(i)=sum(x_comp*a(:,i));
    for j = 1:20
        R(j,i) = sum(a(:,j).*a(:,i));
    end
end

%find coefficients
c = inv(R)*p;
x_hat = c'*a';
e = x_comp-x_hat;

% for i=1:20
%     A(:,:,i) = reshape(a(i,:,:),256,256,1);
% end
% x = reshape(x_comp, 256,256,1);

%   There are two matrices, one containing 20
%   basis images, and one containing the image
%   x that we are trying to decode.

%   Display the 20 basis images in matrix A
% figure;
% for kk = 1:20
%     subplot(4,5,kk);
%     imshow(A(:,:,kk),[]);
% end

%   And display the image x
% figure;
% imshow(x,[]);
e = e/max(e);
error = reshape(e, 256,256);
figure;
imshow(error,[]);

%   The image you are trying to decode is hidden
%   in x.  In order to decode it, you need to
%   decompose x = x_hat + e, where x_hat is an
%   orthogonal projection of x onto the space
%   spanned by the 20 images in A, and e is the
%   error vector between x and its projection.
%   The hidden image will be the error vector e.

%   Some useful hints:

%   (1) To deal with images as vectors, you should
%       "flatten" each image into a column vector that
%       is 256*256 = 65,536 elements long.
%       The "reshape" command in Matlab will come
%       in handy here.
%
%   (2) Once you are done with your vector math, you
%       will probably want to reshape your result
%       back to a 256 x 256 image so you can display it.
%
%   (3) Make sure your images are either scaled from
%       0 to 1 if you are calling "imshow(image)", or
%       that you use the form "imshow(image, [])"
%       where the "[]" tells imshow to scale black
%       and white values between the minimum and
%       maximum values in the image.
