function [gammanew_out,loc_re_new_out,NLL_output,G_re_out] = FIST_optimize_v2(gamma_init,loc_re_init,SMLM_img_re,b,MaxIt,Lmax,imgPare)
r = 58.5/100/2*2;

f_forwardModel = @(x,G_re) abs((G_re*x+b));
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
f_gradient = @(Iobs,Iest,G_re) (G_re.'*(1-Iobs./(Iest+10^-16)));
N = length(gamma_init)/15;
temp = ones(3,1)*100;


MaxIt = 40;
t = 1;
i=1;
Lmax = 0.0001;
NLL_cur_saveFista = [10^10];
%l_it = Lmax;
%%
gammaold_save = zeros(15,size(imgPare.Bx,3));
NLL_cur_save = zeros(1,size(imgPare.Bx,3));
loc_re_old_save = zeros(3,size(imgPare.Bx,3));
for ii = 0:size(imgPare.Bx,3)-1
%imgPare.steps = 'psudo';
loc_re_old = [0,0,(ii-4)*imgPare.pix_sizez].';
[G_re,~,~] = update_basisMatrix(N,zeros(15,1),loc_re_old,imgPare);
gammaold = (pinv(G_re)*(SMLM_img_re-b));
[G_re,loc_re_old,gammaold] = update_basisMatrix(N,gammaold,loc_re_old,imgPare);
gammaold = f_projection(gammaold,r);
z = gammaold;
t = 1;
NLL_cur_saveFista = [10^10];
i=1;
 while (i < MaxIt)

     
     I_est = f_forwardModel(z,G_re);
     I_est(I_est<0) = mean(b)/2;
     %I_est(I_est<0) = SMLM_img_re(I_est<0);
     %I_est = max(0,I_est);
     NLL_cur = f_loss(SMLM_img_re,I_est);
     gradient_cur = f_gradient(SMLM_img_re,I_est,G_re);

     if NLL_cur<min(NLL_cur_saveFista)
        gammanew_out = gammaold;
        G_re_out = G_re;
        loc_re_new_out = loc_re_old;
        NLL_cur_out = NLL_cur;
     end
     
     %backtracking
     l_it = Lmax;
     ik = 1;
     eta = 1.1;
     
     NLL_cur_saveFista = [NLL_cur_saveFista,NLL_cur];
     
     % check ealy stop
     if i>=20     
         %convergece = abs(NLL_cur_save(2:end-1)-NLL_cur_save(3:end))./abs(NLL_cur_save(1:end-2)-NLL_cur_save(3:end));
         decrease = abs((NLL_cur_save(1:end-1)-NLL_cur_save(2:end)));
         if mean(decrease(end-9:end))<5*10^-3 
            break;
         end
     end
     
     z_update = f_projection(z-1/l_it.*gradient_cur,r);
     I_est_update = f_forwardModel(z_update,G_re);
     I_est_update(I_est_update<0) = mean(b)/2;
     %I_est_update(I_est_update<0) = SMLM_img_re(I_est_update<0);
     %I_est_update = max(0,I_est_update);
     f_update = f_loss(SMLM_img_re,I_est_update);
     %comp1 = f_update>NLL_cur-l_it/2*sum((gradient_cur).^2);
     comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_update-z))+l_it/2*sum((z_update-z).^2);%(f-theta*l_it*norm(fg)^2);
     %temp1 = NLL_cur+sum(gradient_cur.*(-1/l_it.*gradient_cur))+l_it/2*sum((1/l_it.*gradient_cur).^2);
     %temp2 = f_update;
     while comp1
         l_it = (eta.^ik)* l_it;
         z_update = f_projection(z-1/l_it.*gradient_cur,r);
         I_est_update = f_forwardModel(z_update,G_re);
         I_est_update(I_est_update<0) = mean(b);
         f_update = f_loss(SMLM_img_re,I_est_update);
         comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_update-z))+l_it/2*sum((z_update-z).^2);%(f-theta*l_it*norm(fg)^2);
         ik = 1+ik;             
     end

     l_it = (eta.^ik)* l_it;
     gammanew = f_projection(z-1/l_it.*gradient_cur,r);
     
     %update grid points and basis  
    [G_re,loc_re_new,gammanew] = update_basisMatrix(N,gammanew,loc_re_old,imgPare);
    % update gammaold with respect to the new grid point        
    S_scdM = gammaold(1:6);
    M_dx_new = S_scdM(1:3)*(loc_re_old(1)-loc_re_new(1))/100+gammaold(7:9);
    M_dy_new = S_scdM(1:3)*(loc_re_old(2)-loc_re_new(2))/100+gammaold(10:12);
    M_dz_new = S_scdM(1:3)*(loc_re_old(3)-loc_re_new(3))/100+gammaold(13:15);
    gammaold = [S_scdM;M_dx_new;M_dy_new;M_dz_new];
    loc_re_old = loc_re_new;


     t_new = 1 / 2 + (sqrt(1 + 4 * t^2)) / 2;
     z = gammanew + ((t - 1) / t_new) * (gammanew - gammaold);
%      z2 = f_projection(z,r);
%      if sum(z2-z)~=0
%          z = gammanew; t=1;
%      end
     gammaold = gammanew;
     t = t_new;
     i = i + 1;   
     
    
     
 end

gammaold_save(:,ii+1) = gammanew_out;
loc_re_old_save(:,ii+1) = loc_re_new_out;
NLL_cur_save(ii+1) = NLL_cur_out;
end
%%
[~,~,C]=mode(loc_re_old_save(3,:));
C = C{1};
loc_re_old2 = [];
gammaold_save2 = [];
NLL_cur_save2 = [];
% for ii = 1:length(C)
%     NLL_cur_save2 = [NLL_cur_save2,NLL_cur_save(loc_re_old_save(3,:)==C(ii))];
%     gammaold_save2 = [gammaold_save2,gammaold_save(:,loc_re_old_save(3,:)==C(ii))];
%     loc_re_old2 = [loc_re_old2,loc_re_old_save(:,loc_re_old_save(3,:)==C(ii))];
% end
% [~,indx] = min(NLL_cur_save2);
% gammaold = gammaold_save2(:,indx);
%loc_re_old = loc_re_old2(:,indx);
[~,indx] = min(NLL_cur_save);
gammaold = gammaold_save(:,indx);
loc_re_old = loc_re_old_save(:,indx);

%loc_re_old = loc_re_init;
%gammaold = gamma_init;
if gammaold(3)<0
    aa = 1
end
z = gammaold;
NLL_cur_saveFista = [10^10];
MaxIt = 20;
t = 1;
i=1;
Lmax = 0.0001;

[G_re,~,~] = update_basisMatrix(N,gammaold,loc_re_old,imgPare);
%imgPare.steps = 'others';

%% constrain
% A = zeros(15,15); b = zeros(15,1);
% A(1,1)=1;
% A(2,2)=1;
% A(3,3)=1;
% A(4,1:3)=-1/2;A(4,4)=1;
% A(5,1:3)=-1/2;A(5,4)=-1;
% A(6,1:3)=-1/2;A(6,5)=1;
% A(7,1:3)=-1/2;A(7,5)=-1;
% A(8,1:3)=-1/2;A(8,6)=1;
% A(9,1:3)=-1/2;A(9,6)=-1;

r = 58.5/100/2*8;
%%
 while (i < MaxIt)
%      if (mod(i, 20) == 0 && i < 250)
%             l_it = l_it / 5;
%      end
     
     if (mod(i, 30) == 0 && i < 250)
            r = r / 2;
     end
     
     I_est = f_forwardModel(z,G_re);
     I_est(I_est<0) = mean(b)/2;
     %I_est(I_est<0) = SMLM_img_re(I_est<0);
     %I_est = max(0,I_est);
     NLL_cur = f_loss(SMLM_img_re,I_est);
     gradient_cur = f_gradient(SMLM_img_re,I_est,G_re);

     if NLL_cur<min(NLL_cur_saveFista)
        gammanew_out = gammaold;
        G_re_out = G_re;
        loc_re_new_out = loc_re_old;
        NLL_cur_out = NLL_cur;
     end
     
     %backtracking
     l_it = Lmax;
     ik = 1;
     eta = 1.1;
     
     NLL_cur_saveFista = [NLL_cur_saveFista,NLL_cur];
     
     % check ealy stop
     if i>=20     
         %convergece = abs(NLL_cur_save(2:end-1)-NLL_cur_save(3:end))./abs(NLL_cur_save(1:end-2)-NLL_cur_save(3:end));
         decrease = abs((NLL_cur_save(1:end-1)-NLL_cur_save(2:end)));
         if mean(decrease(end-9:end))<5*10^-3 
            break;
         end
     end
     
     z_update = f_projection(z-1/l_it.*gradient_cur,r);
     I_est_update = f_forwardModel(z_update,G_re);
     I_est_update(I_est_update<0) = mean(b)/2;
     %I_est_update(I_est_update<0) = SMLM_img_re(I_est_update<0);
     %I_est_update = max(0,I_est_update);
     f_update = f_loss(SMLM_img_re,I_est_update);
     %comp1 = f_update>NLL_cur-l_it/2*sum((gradient_cur).^2);
     comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_update-z))+l_it/2*sum((z_update-z).^2);%(f-theta*l_it*norm(fg)^2);
     %temp1 = NLL_cur+sum(gradient_cur.*(-1/l_it.*gradient_cur))+l_it/2*sum((1/l_it.*gradient_cur).^2);
     %temp2 = f_update;
     while comp1
         l_it = (eta.^ik)* l_it;
         z_update = f_projection(z-1/l_it.*gradient_cur,r);
         I_est_update = f_forwardModel(z_update,G_re);
         I_est_update(I_est_update<0) = mean(b);
         f_update = f_loss(SMLM_img_re,I_est_update);
         comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_update-z))+l_it/2*sum((z_update-z).^2);%(f-theta*l_it*norm(fg)^2);
         ik = 1+ik;             
     end

     l_it = (eta.^ik)* l_it;
     gammanew = f_projection(z-1/l_it.*gradient_cur,r);
     
     %update grid points and basis  
    [G_re,loc_re_new,gammanew] = update_basisMatrix(N,gammanew,loc_re_old,imgPare);
    % update gammaold with respect to the new grid point        
    S_scdM = gammaold(1:6);
    M_dx_new = S_scdM(1:3)*(loc_re_old(1)-loc_re_new(1))/100+gammaold(7:9);
    M_dy_new = S_scdM(1:3)*(loc_re_old(2)-loc_re_new(2))/100+gammaold(10:12);
    M_dz_new = S_scdM(1:3)*(loc_re_old(3)-loc_re_new(3))/100+gammaold(13:15);
    gammaold = [S_scdM;M_dx_new;M_dy_new;M_dz_new];
    loc_re_old = loc_re_new;


     
     
     
     t_new = 1 / 2 + (sqrt(1 + 4 * t^2)) / 2;
     z = gammanew + ((t - 1) / t_new) * (gammanew - gammaold);
%      z2 = f_projection(z,r);
%      if sum(z2-z)~=0
%          z = gammanew; t=1;
%      end
     gammaold = gammanew;
     t = t_new;
     i = i + 1;   
     
    
     
 end
     

loc_re_new_out = loc_re_new_out+loc_re_init;
NLL_output = NLL_cur_out;
end








