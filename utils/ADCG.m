function [gamma_est,loc_re_est] = ADCG(gamma_est,loc_rec_est,N_max,SMLM_img_re,b_re,imgPare)

f_forwardModel = @(x,G_re) abs((G_re*x+b));
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
f_gradient = @(Iobs,Iest,G_re) (G_re.'*(1-Iobs./(Iest+10^-16)));



[gamma_est,loc_re_est,G_re] = FIST_optimize_step2_v2(gamma_est,loc_rec_est,SMLM_img_re,b_re,imgPare);
I_est = f_forwardModel(gamma_est,G_re);
NLL_cur = f_loss(SMLM_img_re,I_est);
gradient_cur = f_gradient(SMLM_img_re,I_est,G_re);
     
while size(gamma_est,2)<N_max+1
    % early terminate condition 
    
    
    % add next SM
    gamma_add
    
end


end