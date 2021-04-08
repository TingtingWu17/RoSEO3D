function [loc_re_new_out,gamma_out] = psudo_inverse(gamma_init,loc_re_init,SMLM_img_re,b,MaxIt,Lmax,imgPare)

f_forwardModel = @(x,G_re) abs((G_re*x+b));
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
f_gradient = @(Iobs,Iest,G_re) (G_re.'*(1-Iobs./(Iest+10^-16)));
N = length(gamma_init)/15;
temp = ones(3,1)*100;


MaxIt = 100;
t = 1;
i=1;
Lmax = 0.0001;
l_it = Lmax;

%l_it = Lmax;
gammaold_save = zeros(15,size(imgPare.Bx,3));
NLL_cur_save = zeros(1,size(imgPare.Bx,3));
loc_re_old_save = zeros(3,size(imgPare.Bx,3));
for ii = 0:size(imgPare.Bx,3)-1
%imgPare.steps = 'psudo';
loc_re_old = [0,0,ii*imgPare.pix_sizez].';
[G_re,~,~] = update_basisMatrix(N,zeros(15,1),loc_re_old,imgPare);
gammaold = (pinv(G_re)*(SMLM_img_re-b));
count = 0;
while sum(abs(loc_re_old-temp))~=0 && count<10
count = count+1;
[G_re,loc_re_new,gammanew] = update_basisMatrix(N,gammaold,loc_re_old,imgPare);
gammaold = (pinv(G_re)*(SMLM_img_re-b));
temp = loc_re_old;
loc_re_old = loc_re_new;
end


gammaold_save(:,ii+1) = gammaold;
loc_re_old_save(:,ii+1) = loc_re_old;
I_est = f_forwardModel(gammaold,G_re);
NLL_cur_save(ii+1) = f_loss(SMLM_img_re,I_est);


end
[~,indx] = min(NLL_cur_save);
gammaold = gammaold_save(:,indx);
loc_re_old = loc_re_old_save(:,indx);


     

loc_re_new_out = loc_re_old+loc_re_init;
gamma_out=gammaold;
end




