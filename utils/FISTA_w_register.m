function [gammanew_out,loc_re_new_out,NLL_cur_out] = FISTA_w_register(SMLM_img_re,b,gammaold,loc_re_old,r,r_R,G_re,MaxIt,Lmax,imgPara)

f_forwardModel = @(x,G_re) abs((G_re*x+b));
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
f_gradient = @(Iobs,Iest,G_re) (G_re.'*(1-Iobs./(Iest+10^-16)));



N = length(gammaold)/21;

z = gammaold;
t = 1;
i=1;
NLL_cur_save = [10^10];
 while (i < MaxIt)

     
     I_est = f_forwardModel(z,G_re);
     NLL_cur = f_loss(SMLM_img_re,I_est);
     gradient_cur = f_gradient(SMLM_img_re,I_est,G_re);

     
     %backtracking
     l_it = Lmax;
     ik = 1;
     eta = 1.1;
     if NLL_cur<min(NLL_cur_save)
        gammanew_out = gammaold;
        loc_re_new_out = loc_re_old;
        NLL_cur_out = NLL_cur;
     end
     NLL_cur_save = [NLL_cur_save,NLL_cur];

     
     % check ealy stop
     if i>=20     
         decrease = abs((NLL_cur_save(1:end-1)-NLL_cur_save(2:end)));
         if mean(decrease(end-9:end))<5*10^-3 
            break;
         end
     end
     
     z_update = f_projection_w_register(z-1/l_it.*gradient_cur,r,r_R);
     I_est_update = f_forwardModel(z_update,G_re);
     f_update = f_loss(SMLM_img_re,I_est_update);
     comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_update-z))+l_it/2*sum((z_update-z).^2);
     while comp1
         l_it = (eta.^ik)* l_it;
         z_update = f_projection_w_register(z-1/l_it.*gradient_cur,r,r_R);
         I_est_update = f_forwardModel(z_update,G_re);
         f_update = f_loss(SMLM_img_re,I_est_update);
         comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_update-z))+l_it/2*sum((z_update-z).^2);
         ik = 1+ik;             
     end

     l_it = (eta.^ik)* l_it;
     gammanew = f_projection_w_register(z-1/l_it.*gradient_cur,r,r_R);
     
     %update grid points and basis  
    [G_re,loc_re_new,gammanew] = update_basisMatrix_w_register(N,gammanew,loc_re_old,imgPara);
    % update gammaold with respect to the new grid point        
    for ii = 1:N
        gammaold = reshape(gammaold,[],N);
        S_scdM = gammaold(1:6,ii);
        M_dx_new = S_scdM(1:3)*(loc_re_old(1,ii)-loc_re_new(1,ii))/100+gammaold(7:9,ii);
        M_dy_new = S_scdM(1:3)*(loc_re_old(2,ii)-loc_re_new(2,ii))/100+gammaold(10:12,ii);
        M_dz_new = S_scdM(1:3)*(loc_re_old(3,ii)-loc_re_new(3,ii))/100+gammaold(13:15,ii);
        gammaold(:,ii) = [S_scdM;M_dx_new;M_dy_new;M_dz_new;gammaold(16:21,ii)];
        loc_re_old = loc_re_new;
    end
    gammaold = reshape(gammaold,[],1);

     t_new = 1 / 2 + (sqrt(1 + 4 * t^2)) / 2;
     z = gammanew + ((t - 1) / t_new) * (gammanew - gammaold);
     gammaold = gammanew;
     t = t_new;
     i = i + 1;   
     
    
     
 end
        
end
