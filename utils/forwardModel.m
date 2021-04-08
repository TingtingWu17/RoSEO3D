function [I_est,basisImage] = forwardModel(gamma,FimgE,maskz2,imageSize,b,sumNormalizex,sumNormalizey)

    x = reshape(gamma,[],9);
    I_est = zeros(imageSize,imageSize*2);
 for ii = 1:size(x,1)  
     %x = [0*10^-9,0*10^-9,0*10^-9,1/3,1/3,1/3,0,0,0];
     %x = [58.5*10^-9,58.5*10^-9,0,1,0,0,0,0,0];
     [image,basisImage] = forwardModel_1SM(x(ii,:),FimgE,maskz2,imageSize,sumNormalizex,sumNormalizey);
    I_est = I_est+image;   
 %figure();imagesc(image); axis image
 end
 
 I_est = I_est+b;
 

end



function [image,basisImage] = forwardModel_1SM(x,FimgE,maskz2,imageSize,sumNormalizex,sumNormalizey)

pixSize = 58.5; %nm
%x(1:2) = x(1:2)+pixSize; %remove the coordinate difference after fft
FimgExx = FimgE.imgExx;
FimgExy = FimgE.imgExy;
FimgExz = FimgE.imgExz;
FimgEyx = FimgE.imgEyx;
FimgEyy = FimgE.imgEyy;
FimgEyz = FimgE.imgEyz;

N = size(FimgExx,1);
%[Xgrid,Ygrid] = meshgrid(-(N-1)/2:(N-1)/2,-(N-1)/2:(N-1)/2);
[Xgrid,Ygrid] = meshgrid(-1/2:1/(N-1):1/2,-1/2:1/(N-1):1/2);
[phi,rho] = cart2pol(Xgrid*6.8712,Ygrid*6.8712);

theta = asin(rho);
%maskz2 = cos(theta)/1.334;

%maskdx = maskdxy;
%maskdy = rot90(maskdx,3);
maskdx = 2*pi*Xgrid;
maskdy = rot90(maskdx,3);

%postMask = x(1)*maskdx+x(2)*maskdy+x(3)*maskz2;
postMask = x(1)/pixSize*maskdx+x(2)/pixSize*maskdy+x(3)*10^-9*maskz2;
postMask = exp(1i*postMask);

imgExx = rot90(fftshift(fft2(FimgExx.*postMask))',2);
imgExy = rot90(fftshift(fft2(FimgExy.*postMask))',2);
imgExz = rot90(fftshift(fft2(FimgExz.*postMask))',2);
imgEyx = rot90(fftshift(fft2(FimgEyx.*postMask)),2);
imgEyy = rot90(fftshift(fft2(FimgEyy.*postMask)),2);
imgEyz = rot90(fftshift(fft2(FimgEyz.*postMask)),2);

roi = @(img)img(-(imageSize - 1)/2+N/2+1:1:(imageSize - 1)/2+N/2+1, ... .
          -(imageSize - 1)/2+N/2+1:1:(imageSize - 1)/2+N/2+1, :);

basisImagex(:,:,1) = roi(abs(imgExx).^2);
basisImagex(:,:,2) = roi(abs(imgExy).^2);
basisImagex(:,:,3) = roi(abs(imgExz).^2);
basisImagex(:,:,4) = roi(2*real(conj(imgExx).*imgExy));
basisImagex(:,:,5) = roi(2*real(conj(imgExx).*imgExz));
basisImagex(:,:,6) = roi(2*real(conj(imgExy).*imgExz));  


          
basisImagey(:,:,1) = roi(abs(imgEyx).^2);
basisImagey(:,:,2) = roi(abs(imgEyy).^2);
basisImagey(:,:,3) = roi(abs(imgEyz).^2);
basisImagey(:,:,4) = roi(2*real(imgEyx.*conj(imgEyy)));
basisImagey(:,:,5) = roi(2*real(imgEyx.*conj(imgEyz)));
basisImagey(:,:,6) = roi(2*real(imgEyy.*conj(imgEyz)));

      
imgx = basisImagex(:,:,1)*x(4)+basisImagex(:,:,2)*x(5)+basisImagex(:,:,3)*x(6)+basisImagex(:,:,4)*x(7)+...
       basisImagex(:,:,5)*x(8)+basisImagex(:,:,6)*x(9);
imgy = basisImagey(:,:,1)*x(4)+basisImagey(:,:,2)*x(5)+basisImagey(:,:,3)*x(6)+basisImagey(:,:,4)*x(7)+...
       basisImagey(:,:,5)*x(8)+basisImagey(:,:,6)*x(9);

image = cat(2,imgx/sumNormalizex,imgy/sumNormalizey);
basisImage = cat(2, basisImagex/sumNormalizex,basisImagey/sumNormalizey);

end