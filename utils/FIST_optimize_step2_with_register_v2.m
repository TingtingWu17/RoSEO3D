function [gammanew_out,loc_re_new_out] = FIST_optimize_step2_with_register_v2(gamma_init,loc_re_init,SMLM_img_re,b,imgPara)


gamma_init = reshape(gamma_init,[],1);
N = length(gamma_init)/21;

%---------------------------FISTA optimize parameters---------------------------
MaxIt = 100;
gammaold = gamma_init;
Lmax = 0.0001;
loc_re_old = loc_re_init;
r = 58.5/100/2*3;
r_R = 58.5/100/2*2;
[G_re,~,~] = update_basisMatrix_w_register(N,gammaold,loc_re_old,imgPara);

%--------------------------
[gammanew_out,loc_re_new_out,~] = FISTA_w_register(SMLM_img_re,b,gammaold,loc_re_old,r,r_R,G_re,MaxIt,Lmax,imgPara);  
                              
 
end