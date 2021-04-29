function [bkg,SMLM_img,I,SM_GT,angle_GT] = generate_estimation_image

addpath(genpath('C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\PSF-optimization'));
addpath(genpath('C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\data of PSF-optimization'));
%%
addpath(genpath('/Users/tingtingwu/Documents/OneDrive - Washington University in St. Louis/github/PSF-optimization'));
addpath(genpath('/Users/tingtingwu/Documents/OneDrive - Washington University in St. Louis/github/data of PSF-optimization'));
%% Imaging system parameters 
%addpath(genpath('/Users/tingtingwu/Documents/OneDrive - Washington University in St. Louis/github/PSF-optimization'));
%addpath(genpath('C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\PSF-optimization'));

signal=1000;
background=2;

Microscopy=struct();
Microscopy.wavelength = 610e-9; 
Microscopy.n1 = 1.518;
Microscopy.n2 = 1.334;
Microscopy.nh = 1.518;
Microscopy.NA = 1.4;
Microscopy.pixelSizeUpsampling = 1;
Microscopy.upsampling = 1;
Microscopy.pix_size=6500/Microscopy.pixelSizeUpsampling;
Microscopy.bfp_radius = 80*Microscopy.upsampling;
Microscopy.Magnitude = 111.1111;
dx = (Microscopy.pix_size*10^-9)/Microscopy.Magnitude;
Microscopy.sampling_size = round(Microscopy.wavelength*Microscopy.bfp_radius/dx/Microscopy.NA)*Microscopy.pixelSizeUpsampling;
if rem(Microscopy.sampling_size,2)==0
    Microscopy.sampling_size=Microscopy.sampling_size+1;
end
Microscopy.image_size =41*Microscopy.pixelSizeUpsampling; %the final image is a 75 by 75 image, odd number
Microscopy.z = 0; %(unit:m) focal length
Microscopy.z2 = 0; %(unit:m) SM's z position
Microscopy.zh = 0;
Microscopy.xy_ind = 0;

imgSz = Microscopy.image_size;
% import opt PSF    
pmask = NaN;
%coordinate are same as the microscop5y coordinate
Microscopy.mask='pixOL_v12.bmp'; 
%Microscopy.mask='opt_v11.bmp'; 
%Microscopy.mask='clear plane.bmp'; 
%Microscopy.mask='DH_biasAdjusted_stan_rescaled-resized_PupilRadius_80_Rot_90_xyShift_0_0.bmp'; 
%Microscopy.mask='trispot-3pi.bmp'; 
%Microscopy.mask='VIPR_perfect_pixOL.bmp'; 

%%
%SM1
Microscopy.rot=0; Microscopy.z=-450*10^-9; zf_init=0; Microscopy.z2 = 0*10^-9;
[basis_matrix_opt1,mask_opt,BFP] = basis_matrix_pixel_based_v3_in(Microscopy,pmask);
loc = [0,0];
imgSz = Microscopy.image_size;
theta =0; phi = 45; gamma=1; 
[muxx,muyy,muzz,muxy,muxz,muyz] = Quickly_rotating_matrix_angleD_gamma_to_M(theta,phi,gamma);
I1 = basis_matrix_opt1*[muxx,muyy,muzz,muxy,muxz,muyz].';
I1 = signal*I1;
I1 = reshape(I1,imgSz,imgSz*2);
I1x = I1(:,1:imgSz); I1y = fliplr(I1(:,imgSz+1:end));
I1x = imtranslate(I1x,[1,1]); I1y = imtranslate(I1y,[-1,-1]);
I1 = imtranslate([I1x,I1y],loc);
SM1_GT = [loc*58.5,Microscopy.z2*10^9,signal,[muxx,muyy,muzz,muxy,muxz,muyz]];
angle_GT1 = [theta,phi,gamma];

basis_matrix_opt2 = basis_matrix_opt1;
% SM2
Microscopy.rot=0; Microscopy.z=-450*10^-9; zf_init=0; Microscopy.z2 = 350*10^-9;
[basis_matrix_opt2,mask_opt,BFP] = basis_matrix_pixel_based_v3_in(Microscopy,pmask);
loc = [10, 0];

imgSz = Microscopy.image_size;
theta =90; phi = 45; gamma=1; 
[muxx,muyy,muzz,muxy,muxz,muyz] = Quickly_rotating_matrix_angleD_gamma_to_M(theta,phi,gamma);
I2 = basis_matrix_opt2*[muxx,muyy,muzz,muxy,muxz,muyz].';
I2 = signal*I2;
I2 = reshape(I2,imgSz,imgSz*2);
I2x = I2(:,1:imgSz); I2y = fliplr(I2(:,imgSz+1:end));
I2 = imtranslate([I2x,I2y],loc);
angle_GT2 = [theta,phi,gamma];
SM2_GT = [loc*58.5,Microscopy.z2*10^9,signal,[muxx,muyy,muzz,muxy,muxz,muyz] ];
% 
% % SM3
% Microscopy.rot=0; Microscopy.z=-600*10^-9; zf_init=0; Microscopy.z2 = 400*10^-9;
% [basis_matrix_opt2,mask_opt,BFP] = basis_matrix_pixel_based_v3_in(Microscopy,pmask);
% loc = [20, 23];
% 
% imgSz = Microscopy.image_size;
% theta =45; phi = 45; gamma=1; 
% [muxx,muyy,muzz,muxy,muxz,muyz] = Quickly_rotating_matrix_angleD_gamma_to_M(theta,phi,gamma);
% I3 = basis_matrix_opt2*[muxx,muyy,muzz,muxy,muxz,muyz].';
% I3 = signal*I3;
% I3 = reshape(I3,imgSz,imgSz*2);
% I3x = I3(:,1:imgSz); I3y = fliplr(I3(:,imgSz+1:end));
% I3 = imtranslate([I3x,I3y],loc);
% angle_GT3 = [theta,phi,gamma];
% SM3_GT = [loc*58.5,Microscopy.z2*10^9,signal,[muxx,muyy,muzz,muxy,muxz,muyz] ];

% I = I1+I2+I3+background;
I = I1+background;
SMLM_img = poissrnd(I);
%SMLM_img = I;
bkg = ones(size(SMLM_img))*background;


%SM_GT = [SM1_GT;SM2_GT;SM3_GT];
%angle_GT = [angle_GT1;angle_GT2;angle_GT3];
SM_GT = [SM1_GT;SM1_GT];
angle_GT = [angle_GT1;angle_GT1];
end
