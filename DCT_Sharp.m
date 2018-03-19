function [FusedDCTSharp,FusedDCTSharp_CV]= DCT_Sharp(im1, im2)

% C)Mostafa Amin-Naji, Babol Noshirvani University of Technology,
% My Official Website: www.Amin-Naji.com
% My Email: Mostafa.Amin.Naji@Gmail.com

% PLEASE CITE THE BELOW PAPERS IF YOU USE THIS CODE

% M. A. Naji and A. Aghagolzadeh, “A new multi-focus image fusion technique
% based on variance in DCT domain,” in 2015 2nd International Conference on
% Knowledge-Based Engineering and Innovation (KBEI), 2015, pp. 478-484. 
% https://doi.org/10.1109/KBEI.2015.7436092 

% Inputs:
%       im1	:	First source image 
%       im2	:	Second source image
%               
% Outputs:
%       FusedDCTSharp	:	Fused image as the result of "DCT+Sharp" method
%       FusedDCTSharp_CV	:	Fused image as the result of "DCT+Sharp+CV" method
% 
% 
% Sample use:
% im1 = imread('pepsi1.tif');
% im2 = imread('pepsi2.tif');
% [FusedDCTSharp,FusedDCTSharp_CV]= DCT_Sharp(im1, im2);



if nargin ~= 2	% Check the correct number of arguments
    error('There should be two input images!')
end

if size(im1,3) == 3     % Check if the images are grayscale
    im1 = rgb2gray(im1);
end
if size(im2,3) == 3
    im2 = rgb2gray(im2);
end

if size(im1) ~= size(im2)	% Check if the input images are of the same size
    error('Size of the source images must be the same!')
end

% Obtain sharpened images by passing the images through a 
% high-pass filter (Unsharp Filter) 

h = fspecial('unsharp');                                                                                                                                                                                                                                       
im1_sharp=imfilter(im1,h,'symmetric');
im2_sharp=imfilter(im2,h,'symmetric');

% Get input image size
[m,n] = size(im1);
FusedDCTSharp = zeros(m,n);
FusedDCTSharp_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	

% Level shifting
im1 = double(im1)-128;
im2 = double(im2)-128;
im1_sharp = double(im1_sharp)-128;
im2_sharp = double(im2_sharp)-128;

% Divide source images into 8*8 blocks and perform the fusion process
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
       im1s = im1(8*i-7:8*i,8*j-7:8*j);
        im2s = im2(8*i-7:8*i,8*j-7:8*j);
        im1Sub = im1_sharp(8*i-7:8*i,8*j-7:8*j);
        im2Sub = im2_sharp(8*i-7:8*i,8*j-7:8*j);
        % Compute the 2-D DCT of 8*8 blocks
        im1sdct = dct2(im1s);
        im2sdct = dct2(im2s);
        im1SubDct = dct2(im1Sub);
        im2SubDct = dct2(im2Sub);
        % Calculate normalized transform coefficients
        im1Norm = im1SubDct ./ 8;
        im2Norm = im2SubDct ./ 8;
        
        % Mean value of 8*8 block of images (Measure for surrounding lumminance)
        im1Mean = im1Norm(1,1);
        im2Mean = im2Norm(1,1);
        
        % Variance of 8*8 block of images
        im1Var = sum(sum(im1Norm.^2)) - im1Mean.^2;
        im2Var = sum(sum(im2Norm.^2)) - im2Mean.^2;
        

          t=60;
        if im1Var > (im2Var+t)
            dctSub = im1sdct;
            Map(i,j) =-1;	% Consistency verification 
       
        end
if            im1Var < (im2Var-t)
            dctSub = im2sdct;
            Map(i,j) = +1;    % Consistency verification 
end
% 
        if im1Var < (im2Var+t)&& im1Var > (im2Var-t)
            Map(i,j)=0;

% dctVarSub=((im1Var)/(im1Var+im2Var)).*im1sdct+((im2Var)/(im1Var+im2Var)).*im2sdct;
            dctSub = (im1sdct+im2sdct)./2;

        end
        

        % Compute the 2-D inverse DCT of 8*8 blocks and construct fused image
        FusedDCTSharp(8*i-7:8*i,8*j-7:8*j) = idct2(dctSub);	% DCT+Variance method
        
    end
end

% Concistency verification using a Majority Filter
fi = ones(3)/9;	% Filter kernel 7*7

cvMapFiltered = imfilter(Map, fi,'symmetric');	% Filtered index map

for i = 1:m/8
    for j = 1:n/8
        % DCT+Variance+CV method
       

       
     FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) =(((1-cvMapFiltered(i,j))/2)*im1(8*i-7:8*i,8*j-7:8*j))+(((1+cvMapFiltered(i,j))/2)*im2(8*i-7:8*i,8*j-7:8*j));


   
    end
end

for i = 1:m/8
    for j = 1:n/8
        % DCT+Variance+CV method
        if cvMapFiltered(i,j) < -0.06
            FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) = im1(8*i-7:8*i,8*j-7:8*j);

        end
        if cvMapFiltered(i,j) > 0.06
            FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) = im2(8*i-7:8*i,8*j-7:8*j);
       
        end
        
    end
end

% inverse shift 
im1 = uint8(double(im1)+128);
im2 = uint8(double(im2)+128);
FusedDCTSharp = uint8(double(FusedDCTSharp)+128);
FusedDCTSharp_CV = uint8(double(FusedDCTSharp_CV)+128);
figure, imshow(im1), title('source image 1');
figure, imshow(im2), title('source image 2');
figure, imshow(FusedDCTSharp), title('dct+variance');
figure, imshow(FusedDCTSharp_CV), title('dct+variance+cv');
