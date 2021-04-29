function [gammanew_out,loc_re_new_out] = FIST_optimize_step2_v2(gamma_init,loc_re_init,SMLM_img_re,b,imgPara)


gamma_init = reshape(gamma_init,[],1);
N = length(gamma_init)/15;

%---------------------------FISTA optimize parameters---------------------------
MaxIt = 100;
gammaold = gamma_init;
Lmax = 0.0001;
loc_re_old = loc_re_init;
r = 58.5/100/2*3;
[G_re,~,~] = update_basisMatrix(N,gammaold,loc_re_old,imgPara);

%--------------------------
[gammanew_out,loc_re_new_out,~] = FISTA(SMLM_img_re,b,gammaold,loc_re_old,r,G_re,MaxIt,Lmax,imgPara);  
                              
 
end




