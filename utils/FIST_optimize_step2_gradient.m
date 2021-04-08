function [gammanew,loc_re_init] = FIST_optimize_step2_gradient(gamma_init,loc_re_init,SMLM_img_re,b,imgPare)
%r = (sqrt((imgPare.pix_sizex/2)^2+(imgPare.pix_sizex/2)^2+(imgPare.pix_sizez/2)^2))/100;
r = 58.5/100/2/40;
rec_val = 0.22;

f_forwardModel = @(x,G_re) abs((G_re*x+b));
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
f_gradient = @(Iobs,Iest,G_re) (G_re.'*(1-Iobs./(Iest+10^-16)));
f_progect = @(gamma) projection_cone3(gamma, r);
%f_progect = @(gamma)gamma;
gammanew_update = @(x, grad, l_it) f_progect(x-1/l_it*(grad + gradhmu3(x, rec_val)));

gamma_init = reshape(gamma_init,[],1);
N = length(gamma_init)/15;


gammaold = gamma_init;
loc_re_old = loc_re_init;
z = gamma_init;

MaxIt = 100;
t = 1;
i=1;
Lmax = 0.0001;
l_it = Lmax;



[G_re,~,~] = update_basisMatrix(N,gammaold,loc_re_init,imgPare);

NLL_cur_save = [10^9];
 while (i < MaxIt)

     
     if (mod(i, 20) == 0 && i < 250)
            l_it = l_it / 5;
     end
     I_est = f_forwardModel(z,G_re);
     I_est(I_est<0) = mean(b)/2;
     
     NLL_cur = f_loss(SMLM_img_re,I_est);
     gradient_cur = f_gradient(SMLM_img_re,I_est,G_re);
     
%      if NLL_cur<min(NLL_cur_save)
%         gammanew_out = gammaold;
%         G_re_out = G_re;
%         loc_re_new_out = loc_re_new;
%         NLL_cur_out = NLL_cur;
%      end
     
     %backtracking
     ik = 1;
     eta = 1.1;
     NLL_cur_save = [NLL_cur_save,NLL_cur];
     z_pjt_update = f_progect(z-1/l_it.*gradient_cur);
     I_est_update = f_forwardModel(z_pjt_update,G_re); 
     %I_est_update(I_est_update<0) = mean(b)/2;
     f_update = f_loss(SMLM_img_re,I_est_update);
     
      % check ealy stop
     if i>=20     
         %convergece = abs(NLL_cur_save(2:end-1)-NLL_cur_save(3:end))./abs(NLL_cur_save(1:end-2)-NLL_cur_save(3:end));
         decrease = abs((NLL_cur_save(1:end-1)-NLL_cur_save(2:end)));
         if mean(decrease(end-9:end))<5*10^-3 
            break;
         end
     end
     

     comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_pjt_update-z))+l_it/2*sum((z_pjt_update-z).^2);%(f-theta*l_it*norm(fg)^2);
     %temp1 = NLL_cur+sum(gradient_cur.*(1/l_it.*gradient_cur))+l_it/2*sum((1/l_it.*gradient_cur).^2);
     %temp2 = f_update;
     while comp1
         l_it = (eta.^ik)* l_it;
         z_pjt_update = f_progect(z-1/l_it.*gradient_cur);
         I_est_update = f_forwardModel(z_pjt_update,G_re);
         
         %I_est_update(I_est_update<0) = mean(b)/2;
         %I_est_update = max(0,I_est_update);
         f_update = f_loss(SMLM_img_re,I_est_update);
         comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_pjt_update-z))+l_it/2*sum((z_pjt_update-z).^2);%(f-theta*l_it*norm(fg)^2);
         ik = 1+ik;             
     end

     l_it = (eta.^ik)* l_it+1;
     gammanew = gammanew_update(z, gradient_cur, l_it);
     
     %update grid points and basis
     
         
%         [G_re,loc_re_new,gammanew] = update_basisMatrix(N,gammanew,loc_re_old,imgPare);
%         % update gammaold with respect to the new grid point  
%         for ii = 1:N
%             gammaold = reshape(gammaold,[],N);
%             S_scdM = gammaold(1:6,ii);
%             M_dx_new = S_scdM(1:3)*(loc_re_old(1,ii)-loc_re_new(1,ii))/100+gammaold(7:9,ii);
%             M_dy_new = S_scdM(1:3)*(loc_re_old(2,ii)-loc_re_new(2,ii))/100+gammaold(10:12,ii);
%             M_dz_new = S_scdM(1:3)*(loc_re_old(3,ii)-loc_re_new(3,ii))/100+gammaold(13:15,ii);
%             gammaold(:,ii) = [S_scdM;M_dx_new;M_dy_new;M_dz_new];
%             loc_re_old = loc_re_new;
%         end
%         gammaold = reshape(gammaold,[],1);
%  
     
     
     
     
     t_new = 1 / 2 + (sqrt(1 + 4 * t^2)) / 2;
     z = gammanew + ((t - 1) / t_new) * (gammanew - gammaold);
     gammaold = gammanew;
     t = t_new;
     i = i + 1;   
     
     
 end

end