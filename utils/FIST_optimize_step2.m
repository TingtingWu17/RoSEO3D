function gammanew = FIST_optimize_step2(gamma_init,gamma_grid,loc_grid,loc_re_init,SMLM_img_re,b,imgPare,step)
r = (sqrt((imgPare.pix_sizex/2)^2+(imgPare.pix_sizex/2)^2+(imgPare.pix_sizez/2)^2))/100;
rec_val = 0.22;

f_forwardModel = @(x,G_re) abs((G_re*x.'+b));
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
f_gradient = @(Iobs,Iest,G_re) ((G_re.'*(1-Iobs./(Iest+10^-16))).');
f_progect = @(gamma) projection_cone2(gamma, r);
gammanew_update = @(x, grad, l_it) f_progect(x-1/l_it*(grad + gradhmu3(x, rec_val)));
%gammanew_update = @(x, grad, l_it) f_progect(x-1/l_it*(grad));
 
gamma_init = reshape(gamma_init,1,[]);
gamma_grid = reshape(gamma_grid,1,[]);


N = length(gamma_init)/15;



MaxIt = 100;
gammaold = gamma_grid;
z = gamma_grid;
t = 1;
i=1;
Lmax = 0.001;
l_it = Lmax;


%loc_re_old = [0,0,500];

G_re = update_basisMatrix_step2(loc_grid,imgPare);

% psudo inverse to get initialize value
%gammaold = (pinv(G_re)*(SMLM_img_re-b))';
z = gammaold;


NLL_cur_save = [];
 while (i < MaxIt)
     if i==90
         aaa=1;
     end
     
     if i==50
         aaa=1;
     end
     
     
     I_est = f_forwardModel(z,G_re);
     %I_est(I_est<0) = SMLM_img_re(I_est<0);
     %I_est = max(0,I_est);
     NLL_cur = f_loss(SMLM_img_re,I_est);
     gradient_cur = f_gradient(SMLM_img_re,I_est,G_re);

     %backtracking
     ik = 1;
     eta = 1.1;
     %l_it = Lmax;
     NLL_cur_save = [NLL_cur_save,NLL_cur];
     z_pjt_update = f_progect(z-1/l_it.*gradient_cur);
     I_est_update = f_forwardModel(z_pjt_update,G_re);
     %I_est_update(I_est_update<0) = SMLM_img_re(I_est_update<0);
     %I_est_update = max(0,I_est_update);
     f_update = f_loss(SMLM_img_re,I_est_update);
     %comp1 = f_update>NLL_cur-l_it/2*sum((gradient_cur).^2);
     comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_pjt_update-z))+l_it/2*sum((z_pjt_update-z).^2);%(f-theta*l_it*norm(fg)^2);
     temp1 = NLL_cur+sum(gradient_cur.*(z_pjt_update-z))+l_it/2*sum((z_pjt_update-z).^2);
     temp2 = f_update;
     while comp1
         l_it = (eta.^ik)* l_it;
         z_pjt_update = f_progect(z-1/l_it.*gradient_cur);
         I_est_update = f_forwardModel(z_pjt_update,G_re);
         %I_est_update = max(0,I_est_update);
         f_update = f_loss(SMLM_img_re,I_est_update);
         comp1 = f_update>NLL_cur+sum(gradient_cur.*(z_pjt_update-z))+l_it/2*sum((z_pjt_update-z).^2);%(f-theta*l_it*norm(fg)^2);
         ik = 1+ik;             
     end

     l_it = (eta.^ik)* l_it+1;
     %gammanew = z-1/l_it.*gradient_cur;
     gammanew = gammanew_update(z, gradient_cur, l_it);
     
     %update grid points and basis
     t_new = 1 / 2 + (sqrt(1 + 4 * t^2)) / 2;
     z = gammanew + ((t - 1) / t_new) * (gammanew - gammaold);
     gammaold = gammanew;
     t = t_new;
     i = i + 1;   
     
     
 end
     

end



function G_re = update_basisMatrix_step2(loc_rec,imgPare)


N = size(loc_rec,1);
img_size = imgPare.img_size;
pix_sizex = imgPare.pix_sizex;
pix_sizez = imgPare.pix_sizez;
Bx = imgPare.Bx;
By = imgPare.By;
B_size = size(Bx,1);
ind = (B_size-1)/2+1+[-(img_size-1)/2:(img_size-1)/2];
Bx = Bx(ind,ind,:,:);
By = By(ind,ind,:,:);

Gx = zeros(img_size,img_size,15*N);
Gy = zeros(img_size,img_size,15*N);



for ii = 1:N 
    ind_z = min(max(loc_rec(ii,3)/pix_sizez,1),10);
    ind_x = loc_rec(ii,1)/pix_sizex;
    ind_y = loc_rec(ii,2)/pix_sizex;      
    
    Gx(:,:,15*(ii-1)+1:15*ii) = imtranslate(squeeze(Bx(:,:,ind_z,:)),[ind_x,ind_y]);
    Gy(:,:,15*(ii-1)+1:15*ii) = imtranslate(squeeze(By(:,:,ind_z,:)),[ind_x,ind_y]);
end

G = cat(2,Gx,Gy);
G_re = reshape(G,[],15*N);

end