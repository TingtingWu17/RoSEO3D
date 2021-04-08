function SM_est = RoSEO3D(obj, SMLM_img, b, PSFx, PSFy,n1,FPSFx, FPSFy, imgPara,varargin)


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

%% maximum peak 1
% num_iter=5;
% wavelet_level=2;
% img = SMLM_img(:,1:img_size)+SMLM_img(:,(1:img_size)+img_size);
% Energy=Wavelet_backg_est(reshape(img,img_size,img_size,1),...
%     mean(mean(b(:,1:img_size)+b(:,(1:img_size)+img_size))),wavelet_level,'db3',num_iter);
% Energy1 = Energy;
% Energy1(Energy<=.5*max(Energy(:)))=0;
% maxNumBeads=5;
% PSFSize=9;
% LMaxFinder =vision.LocalMaximaFinder(maxNumBeads,[PSFSize PSFSize]);
% lc_max=(single(LMaxFinder(Energy1))-[(img_size-1)/2,(img_size-1)/2]).';
% loc_rec = cat(1,lc_max*pix_sizex,zeros(1,size(lc_max,2)));

%% initialize the gamma v1
%initialize the gamma
%------------------------------------------------------------
%
% forward_model = @(x,G_re) G_re*x+reshape(b,[],1);
% NLL_cal = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest));
% gradient_cal = @(x,Iobs,Iest,G_re) G_re.'*(1-Iobs./Iest);
% 
% 
% Gx = zeros(img_size,img_size,15*size(lc_max,1));
% Gy = zeros(img_size,img_size,15*size(lc_max,1));
% %gamma_init = zeros(27*size(lc_max,1),1);
% 
% 
% Niter = 1;
% jj = 0;
% while jj<Niter
% for ii = 1:size(lc_max,1)    
%     ind_z = min(max(round(loc_rec(ii,3)/pix_sizez),1),10);
%     ind_x = round(loc_rec(ii,1)/pix_sizex);
%     ind_y = round(loc_rec(ii,2)/pix_sizex);    
%     Gx(:,:,15*(ii-1)+1:15*ii) = imtranslate(squeeze(Bx(:,:,ind_z,:)),[ind_x,ind_y]);
%     Gy(:,:,15*(ii-1)+1:15*ii) = imtranslate(squeeze(By(:,:,ind_z,:)),[ind_x,ind_y]);
% end
% G = cat(2,Gx,Gy);
% G_re = reshape(G,[],15*size(lc_max,1));
% gamma_init = pinv(G_re)*reshape((SMLM_img-b),[],1);
% SM_est= get_loc_data2(gamma_init, loc_rec);
% loc_rec_new = SM_est(:,1:3);
% loc_rec = loc_rec+max(min((loc_rec_new-loc_rec),repmat([pix_sizex,pix_sizex,pix_sizez],size(lc_max,1),1)),repmat(-[pix_sizex,pix_sizex,pix_sizez],size(lc_max,1),1));
% 
% SMLM_img_re = reshape((SMLM_img),[],1);
% jj = jj+1;
% %figure();
% %imagesc(reshape(forward_model(gamma_init,G_re),41,82));
% end
% 
% %gamma_init = reshape(gamma_init,[],15);
% 
% % estimate the orientation and localization based on NLL
% 
%% Step2: estimate orientation for patch by patch
% %initialize the gamma v1
% % ROI selected region
% f_forwardModel = @(x,G_re,b) abs((G_re*x+b));
% f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
% f_gradient = @(Iobs,Iest,G_re) (G_re.'*(1-Iobs./(Iest+10^-16)));
% ß
% gamma_init = [];%zeros(15,size(lc_max,2));
% loc_rec_init = [];%zeros(3,size(lc_max,2));
% NLL_save = [];%zeros(size(lc_max,1),1);
% %SMLM_img_re = reshape((SMLM_img),[],1);
% SMLM_remove = SMLM_img;
% 
% 
% % divide SMs into different patches
% [lc_indx_out_center,indx_out] = define_pathes(lc_max,PSF_size_opt);
% loc_rec_center = cat(1,lc_indx_out_center*pix_sizex,zeros(1,size(lc_indx_out_center,2)));
% for ii = 1:size(lc_indx_out_center,2)
%     %crop the whole patch
%     if size(lc_indx_out_center,2)==1
%       lc_center_cur = loc_rec_center;  
%     else
%     lc_center_cur = loc_rec_center(:,ii);
%     end
%     distance = round(sqrt(sum((repmat(lc_center_cur,1,sum(indx_out==ii))-loc_rec(:,indx_out==ii)).^2,1))/pix_sizex);
%     PSF_size_whole = max(distance)*2+PSF_size_opt;
%     if mod(PSF_size_whole,2)==0
%         PSF_size_whole=PSF_size_whole+1;
%     end
%     lc_ind_cur = lc_center_cur./[pix_sizex,pix_sizex,pix_sizez].'+[(img_size+1)/2,(img_size+1)/2,0].';
%     indx = max(min(lc_ind_cur(1)+[-(PSF_size_whole-1)/2:(PSF_size_whole-1)/2],img_size),1);
%     indy = max(min(lc_ind_cur(2)+[-(PSF_size_whole-1)/2:(PSF_size_whole-1)/2],img_size),1);
%     SMLM_img_crop_whole = cat(2,SMLM_remove(indy,indx),SMLM_remove(indy,indx+img_size));
%     b_crop_whole = cat(2,b(indy,indx),b(indy,indx+img_size));
% 
%     
%     gamma_patch_init = [];
%     loc_rec_patch_init = [];
%     for jj = 1:sum(indx_out==ii)
%         [~,ind] = max(indx_out==ii);
%         lc_cur = loc_rec(:,ind+jj-1);
%         distance_all = round(sqrt(sum((repmat(lc_cur,1,sum(indx_out==ii))-loc_rec(:,indx_out==ii)).^2,1))/pix_sizex);
%         distance_all(distance_all==0) = [];
%         PSF_size = min(PSF_size_opt,max(11,min(distance_all)));
%         if mod(PSF_size,2)==0
%             PSF_size = PSF_size+1;
%         end
% 
%         if size(sum(indx_out==ii),2)==1
%             PSF_size = PSF_size_opt;
%         end
%         %PSF_size = imgPare.img_size;
%         lc_ind_cur = lc_cur./[pix_sizex,pix_sizex,pix_sizez].'+[(img_size+1)/2,(img_size+1)/2,0].';
%         indx = max(min(lc_ind_cur(1)+[-(PSF_size-1)/2:(PSF_size-1)/2],img_size),1);
%         indy = max(min(lc_ind_cur(2)+[-(PSF_size-1)/2:(PSF_size-1)/2],img_size),1);
%         SMLM_img_crop = cat(2,SMLM_remove(indy,indx),SMLM_remove(indy,indx+img_size));
%         b_crop = cat(2,b(indy,indx),b(indy,indx+img_size));
%         imgPare.img_sizex = length(indx);
%         imgPare.img_sizey = length(indy);
%         %psudo inverse
%         SM_est_init_cur = [lc_cur.',sum(SMLM_img_crop(:)-b_crop(:)),1/3,1/3,1/3,0,0,0].';
%         gamma_cur = [SM_est_init_cur(4)*SM_est_init_cur(5:10);zeros(9,1)];
%         [loc_re_cur,gamma_out] = psudo_inverse(gamma_cur,lc_cur,reshape(SMLM_img_crop,[],1),reshape(b_crop,[],1),MaxIt,Lmax,imgPare);
%         gamma_patch_init = [gamma_patch_init;gamma_out];
%         loc_rec_patch_init = [loc_rec_patch_init,loc_re_cur];
%     end
%     
%     imgPare.img_sizex = size(SMLM_img_crop_whole,2)/2;
%     imgPare.img_sizey = size(SMLM_img_crop_whole,1);
%     loc_rec_patch_init = loc_rec_patch_init-repmat(lc_center_cur,1,size(loc_rec_patch_init,2));
%     [gammanew_cur,loc_rec_cur,NLL,G_re] = FIST_optimize_v2(gamma_patch_init,loc_rec_patch_init,lc_center_cur,reshape(SMLM_img_crop_whole,[],1),reshape(b_crop_whole,[],1),MaxIt,Lmax,imgPare);
%     
%     gamma_init = [gamma_init;gammanew_cur];
%     loc_rec_init = [loc_rec_init,loc_rec_cur];
% end
% 
% 
% 
% 
% 
% if length(gamma_init)==0
% SM_est=[];
% NLL=[];
% return;
% end

%% Step2: estimate orientation for ROI selected region
%initialize the gamma v1
% ROI selected region
f_forwardModel = @(x,G_re,b) abs((G_re*x+b));
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
f_gradient = @(Iobs,Iest,G_re) (G_re.'*(1-Iobs./(Iest+10^-16)));

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
    %end
    %imgPare.img_size = img_size;
    %[G_re,~,~] = update_basisMatrix(1,gammanew_cur,loc_rec_cur,imgPare);
    %I_est = f_forwardModel(gammanew_cur,G_re,reshape(b,[],1));
    %figure();
    %imagesc(reshape(G_re*gammanew_cur,imgPare.img_size,imgPare.img_size*2)); axis image;
    %SMLM_remove = reshape(reshape(SMLM_remove,[],1)-I_est,img_size,img_size*2);
    %SMLM_remove(SMLM_remove<0) = mean(mean(b));
end

if length(gamma_init)==0
SM_est=[];
NLL=[];
return;
end

%% step 3: filter SMs with close signal, location and orientation
% SM_est = get_loc_data2(gamma_init, loc_rec_init);
% SM_similar_box = SM_est(1,:);
% gamma_similar_box = gamma_init(:,1);
% loc_rec_similar_box = loc_rec_init(:,1);
% NLL_similar_box = NLL_save(1);
% thred_distance = norm([58.5*5,58.5*5]); %in unit of nm
% indx_similar_box = zeros(size(SM_est,1),1);
% indx_similar_box(1) = 1;
% if size(SM_est,1)>1
% for ii = 2:size(SM_est,1)  
%     add=1;
%     SM_est_cur = SM_est(ii,:);
%     gamma_cur = gamma_init(:,ii);
%     loc_rec_cur = loc_rec_init(:,ii);
%     NLL_cur = NLL_save(ii);
%     % compare with the SMs in similar box
%     for jj = 1:size(SM_similar_box,1)  
%         distance = norm(SM_est_cur(1:2)-SM_similar_box(jj,1:2));
%         
%         if distance<thred_distance
%             indx_similar_box(ii) = jj;
%             add = 0;
%            if NLL_cur<NLL_similar_box(jj)
%                 SM_similar_box(jj,:) = SM_est_cur;
%                 gamma_similar_box(:,jj) = gamma_cur;
%                 loc_rec_similar_box(:,jj) = loc_rec_cur;
%                 NLL_similar_box(jj) = NLL_cur;
%            end
%         end
%     end
%     if add==1
%             indx_similar_box(ii) = size(SM_similar_box,1)+1;
%             SM_similar_box = [SM_similar_box;SM_est_cur];  
%             gamma_similar_box = [gamma_similar_box,gamma_cur];
%             loc_rec_similar_box = [loc_rec_similar_box,loc_rec_cur];
%             NLL_similar_box = [NLL_similar_box,NLL_cur];
%     end
%     
% end
% end

%gamma_est = gamma_similar_box;
%loc_rec_est = loc_rec_similar_box;

gamma_est = gamma_init;
loc_rec_est = loc_rec_init;
%% step 4: conditional gradient to update the whole FoV
% 
imgPara.img_sizex = img_size;
imgPara.img_sizey = img_size;
N_max = size(lc_max,2);

[gamma_est,loc_re_est,G_re,NLL] = FIST_optimize_step2_v2(gamma_est,loc_rec_est,reshape(SMLM_img,[],1),reshape(b,[],1),imgPara);
%figure();
%imagesc(reshape(G_re*gamma_est,imgPare.img_size,imgPare.img_size*2)); axis image;
% 
SM_est = get_loc_data2(gamma_est, loc_rec_est);     

%% step 4: global grid flow to update the whole FoV

% 
% imgPare.img_sizex = img_size;
% imgPare.img_sizey = img_size;
% % at each estimation put 26 molecules
% [loc_grid,gamma_grid]=build_grid_point(loc_rec_init,gamma_init,imgPare);
% [gamma_est,loc_rec_est] = FIST_optimize_step2_gradient(gamma_grid,loc_grid,reshape(SMLM_img,[],1),reshape(b,[],1),imgPare);
% SM_est = get_loc_data2(gamma_est, loc_rec_est);  
%      
% % temp2 = get_loc_data2(gammanew_cur, loc_grid);
% % loc_cur = temp2(:,1:3);
% % 
% 
% figure();
% scatter3(loc_grid(1,:),loc_grid(2,:),loc_grid(3,:)); hold on;
% scatter3(SM_est(:,1),SM_est(:,2),SM_est(:,3));
% plot3([loc_grid(1,:).',SM_est(:,1)].',[loc_grid(2,:).',SM_est(:,2)].',[loc_grid(3,:).',SM_est(:,3)].');
% %     
end



function [loc_grid,gamma_grid]=build_grid_point(loc_rec_init,gamma_init,imgPare)
gamma_init = reshape(gamma_init,15,[]);
pix_sizex = imgPare.pix_sizex;
pix_sizez = imgPare.pix_sizez;
loc_grid = zeros(3,27*size(gamma_init,2));
gamma_grid = zeros(15,27*size(gamma_init,2));
for ii = 1:size(gamma_init,2)
    
S_scdM = gamma_init(1:6,:);

loc_grid(:,27*(ii-1)+1:27*ii) = repmat(loc_rec_init(:,ii),1,27)+generate_loc_grid;
gamma_grid(:,27*(ii-1)+1:27*ii) = [repmat(S_scdM(:,ii),1,27)/27;zeros(9,27)];

end

    function loc_grid_cur = generate_loc_grid()
        loc_grid_cur = zeros(3,27);
        indx = 0;
        for i = [-2,0,2]*pix_sizex
            for j = [-2,0,2]*pix_sizex
                for k = [-2,0,2]*pix_sizez
                    indx = indx+1;
                    loc_grid_cur(:,indx) = [i,j,k].';
                end
            end
        end   
        
    end


gamma_grid(:,27*(ii-1)+14)=0.001*gamma_grid(:,27*(ii-1)+14);
end








