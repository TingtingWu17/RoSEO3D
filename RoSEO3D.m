function SM_est = RoSEO3D(SMLM_img, b, n1,FPSFx, FPSFy, imgPara,varargin)


s = opt2struct(varargin);

if isfield(s, 'PSF_size_opt')
    PSF_size_opt = s.PSF_size_opt;
else
    PSF_size_opt = 21; % number of sub-frames for sub_frame analysis
end

%check for proper input images
%------------------------------------------------------------
if max(SMLM_img(:)) < 2 * mean(b(:))
    SM_est = [];
    return
end
if any(SMLM_img < 0)
    error('input image must not be background subtracted')
end

img_size = size(SMLM_img, 1); % image side-length size

if 2 * img_size ~= size(SMLM_img, 2)

    error('input image must be a square region')

end

%make sure image size is an odd number
if mode(img_size, 2) == 0
    error('side length (in number of camera pixels) of input image must be an odd value')
end

%% build the imaging parameters for the microscopy

imgPara.PSF_size = PSF_size_opt;

%% Step 1: find the local maximum with RoSEO
[GradMap,loc_data] = RoSEO_step1(n1, SMLM_img, b, FPSFx, FPSFy, varargin{1},varargin{2});
maxNumBeads=5;
PSFSize=5;
thred = 8;
LMaxFinder =vision.LocalMaximaFinder(maxNumBeads,[PSFSize PSFSize],'Threshold',7);
lc_max=single(LMaxFinder(GradMap)).'-[(img_size+1)/2+1,(img_size+1)/2+1].';  %the 3 is from the coordinate difference
if size(lc_max,2)>1
lc_max(:,lc_max(:,2)<-(img_size-1)/2 | lc_max(:,1)<-(img_size-1)/2) = [];
end
lc_max = combine_simimar_points(lc_max,thred);
loc_rec = cat(1,lc_max*pix_sizex,zeros(1,size(lc_max,2)));
if isempty(lc_max)
    SM_est = [];
    NLL = [];
    return
end


%% Step2: estimate orientation for ROI selected region
%initialize the gamma v1
% ROI selected region

gamma_init = [];%zeros(15,size(lc_max,2));
loc_rec_init = [];%zeros(3,size(lc_max,2));
NLL_save = [];%zeros(size(lc_max,1),1);
%SMLM_img_re = reshape((SMLM_img),[],1);
SMLM_remove = SMLM_img;
for ii = 1:size(lc_max,2)  
    
    % select ROI
    lc_cur = loc_rec(:,ii);
    distance_all = round(sqrt(sum((repmat(lc_cur,1,size(loc_rec,2))-loc_rec).^2,1))/58.5);
    distance_all(distance_all==0) = [];
    PSF_size = min(21,max(11,min(distance_all)));
    if mod(PSF_size,2)==0
        PSF_size = PSF_size+1;
    end
    
    if size(lc_max,2)==1
        PSF_size = 21;
    end
    %PSF_size = imgPare.img_size;
    lc_ind_cur = lc_cur./[pix_sizex,pix_sizex,pix_sizez].'+[(img_size+1)/2,(img_size+1)/2,0].';
    indx = max(min(lc_ind_cur(1)+[-(PSF_size-1)/2:(PSF_size-1)/2],img_size),1);
    indy = max(min(lc_ind_cur(2)+[-(PSF_size-1)/2:(PSF_size-1)/2],img_size),1);
    SMLM_img_crop = cat(2,SMLM_remove(indy,indx),SMLM_remove(indy,indx+img_size));
    b_crop = cat(2,b(indy,indx),b(indy,indx+img_size));
    imgPara.img_sizex = length(indx);
    imgPara.img_sizey = length(indy);
    
    % FISTA estimate the orientation and 3D localization
    SM_est_init_cur = [lc_cur.',sum(SMLM_img_crop(:)-b_crop(:)),1/3,1/3,1/3,0,0,0].';
    gamma_cur = [SM_est_init_cur(4)*SM_est_init_cur(5:10);zeros(9,1)];
    
    [gammanew_cur,loc_rec_cur,NLL,G_re] = FIST_optimize_v2(gamma_cur,lc_cur,reshape(SMLM_img_crop,[],1),reshape(b_crop,[],1),MaxIt,Lmax,imgPara);
    %if NLL<0.1
    gamma_init = [gamma_init,gammanew_cur];
    loc_rec_init = [loc_rec_init,loc_rec_cur];
    NLL_save = [NLL_save,NLL];

end

if isempty(gamma_init)
    SM_est=[];
return;
end


gamma_est = gamma_init;
loc_rec_est = loc_rec_init;
%% step 3: conditional gradient to update the whole FoV
% 
imgPara.img_sizex = img_size;
imgPara.img_sizey = img_size;
N_max = size(lc_max,2);

[gamma_est,loc_re_est,G_re,NLL] = FIST_optimize_step2_v2(gamma_est,loc_rec_est,reshape(SMLM_img,[],1),reshape(b,[],1),imgPara);
%figure();
%imagesc(reshape(G_re*gamma_est,imgPare.img_size,imgPare.img_size*2)); axis image;
% 
SM_est = get_loc_data2(gamma_est, loc_rec_est);     
end

