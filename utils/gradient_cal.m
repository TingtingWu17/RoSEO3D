function fg = gradient_cal(gamma,I_est, I_obs,imageSize,maskz2,basisImage,FimgE,sumNormalizex,sumNormalizey)
dx = 1/58.5;
dz = 1;
g1 = (1-I_est./(I_obs+10^-9));
g1 = reshape(g1,1,[]);

    x = reshape(gamma,[],9);
    fg = zeros(size(x));

 for ii = 1:size(x,1) 
     x_cur = x(ii,:);
    %[basisImagex,basisImagey] = B_cal(x_cur,FimgE,maskz2,imageSize);   
    fg(ii,4:9) = g1*reshape(basisImage,[],6);
    %x derivative
    x_cur_dx = x_cur; x_cur_dx(1) = x_cur(1)+dx;
    I_dx = forwardModel_1SM(x_cur_dx,FimgE,maskz2,imageSize,sumNormalizex,sumNormalizey);
    fg(ii,1) = g1*reshape((I_dx-I_est)/dx,[],1);
    %y derivative
    x_cur_dy = x_cur; x_cur_dy(2) = x_cur(2)+dx;
    I_dy = forwardModel_1SM(x_cur_dy,FimgE,maskz2,imageSize,sumNormalizex,sumNormalizey);
    fg(ii,2) = g1*reshape((I_dy-I_est)/dx,[],1);
    
    %z derivative
    x_cur_dz = x_cur; x_cur_dz(3) = x_cur(3)+dz;
    I_dz = forwardModel_1SM(x_cur_dz,FimgE,maskz2,imageSize,sumNormalizex,sumNormalizey);
    fg(ii,3) = g1*reshape((I_dz-I_est)/dz,[],1);    

 end

end





function image = forwardModel_1SM(x,FimgE,maskz2,imageSize,sumNormalizex,sumNormalizey)

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


end