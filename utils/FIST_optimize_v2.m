function [gammanew_out,loc_re_new_out,NLL_output] = FIST_optimize_v2(loc_re_init,SMLM_img_re,b,imgPara)
N = length(loc_re_init(:))/3;

%% step 1: psudo + small interation of FISTA to decide good initial points

%---------------------------FISTA optimize parameters---------------------------
MaxIt = 40;
Lmax = 0.0001;
r = 58.5/100/2*2;


gammaold_save = zeros(15,size(imgPara.Bx,3));
NLL_cur_save = zeros(1,size(imgPara.Bx,3));
loc_re_old_save = zeros(3,size(imgPara.Bx,3));

for  ii= 1:length(imgPara.axial_grid_points)

loc_re_old = [0,0,imgPara.axial_grid_points(ii)].';
[G_re,~,~] = update_basisMatrix(N,zeros(15,1),loc_re_old,imgPara);
gammaold = (pinv(G_re)*(SMLM_img_re-b));
[G_re,loc_re_old,gammaold] = update_basisMatrix(N,gammaold,loc_re_old,imgPara);
gammaold = f_projection(gammaold,r);


[gammanew_out,loc_re_new_out,NLL_cur_out] = FISTA(SMLM_img_re,b,gammaold,loc_re_old,r,G_re,MaxIt,Lmax,imgPara); 
gammaold_save(:,ii) = gammanew_out;
loc_re_old_save(:,ii) = loc_re_new_out;
NLL_cur_save(ii) = NLL_cur_out;
end
%---------------------------choose the best initial points from candidates---------------------------

[~,indx] = min(NLL_cur_save);
gammaold = gammaold_save(:,indx);
loc_re_old = loc_re_old_save(:,indx);


%% step 2: FISTA for ROI 



%---------------------------FISTA optimize parameters---------------------------

MaxIt = 40;
Lmax = 0.0001;
r = 58.5/100/2*4;
r_R = 58.5/100/2*2;
gammaold = reshape(gammaold,15,[]);

[G_re,~,~] = update_basisMatrix(N,gammaold,loc_re_old,imgPara);

[gammanew_out,loc_re_new_out,NLL_cur_out] = FISTA(SMLM_img_re,b,gammaold,loc_re_old,r,G_re,MaxIt,Lmax,imgPara);      

%---------------------------output---------------------------
loc_re_new_out = loc_re_new_out+loc_re_init;
NLL_output = NLL_cur_out;
end








