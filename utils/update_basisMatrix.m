function [G_re,loc_rec,gamma_new] = update_basisMatrix(N,gamma, loc_rec,imgPare)

SM_est= get_loc_data2(gamma, loc_rec);
loc_rec_new_all = (SM_est(:,1:3)).';
loc_rec_adds = loc_rec_new_all-loc_rec;

img_sizex = imgPare.img_sizex;
img_sizey = imgPare.img_sizey;
pix_sizex = imgPare.pix_sizex;
pix_sizez = imgPare.pix_sizez;
Bx = imgPare.Bx;
By = imgPare.By;
B_size = size(Bx,1);


Gx = zeros(img_sizey,img_sizex,15*N);
Gy = zeros(img_sizey,img_sizex,15*N);
%gamma_init = zeros(27*size(lc_max,1),1);
gamma = reshape(gamma,15,[]);
gamma_new = zeros(15,N);
z_slice = size(Bx,3);

for ii = 1:N 
    ind_z = min(max(round(loc_rec_adds(3,ii)/pix_sizez)+loc_rec(3,ii)/pix_sizez,-4),z_slice-5);
    ind_x = max(min(round(loc_rec_adds(1,ii)/pix_sizex)+loc_rec(1,ii)/pix_sizex,img_sizex),-img_sizex);
    ind_y = max(min(round(loc_rec_adds(2,ii)/pix_sizex)+loc_rec(2,ii)/pix_sizex,img_sizey),-img_sizey);  
    
    ind_z_add = ind_z-loc_rec(3,ii)/pix_sizez;
    ind_x_add = ind_x-loc_rec(1,ii)/pix_sizex;
    ind_y_add = ind_y-loc_rec(2,ii)/pix_sizex;

    loc_rec(1:3,ii) = [ind_x_add,ind_y_add,ind_z_add].'.*[pix_sizex,pix_sizex,pix_sizez].'+loc_rec(:,ii);
    gamma_new(:,ii) = gamma_update_2_loc_rec([ind_x_add,ind_y_add,ind_z_add].*[pix_sizex,pix_sizex,pix_sizez],gamma(:,ii));
    
    indx_trans = (B_size-1)/2+1+[-(img_sizex-1)/2:(img_sizex-1)/2]-ind_x;
    indy_trans = (B_size-1)/2+1+[-(img_sizey-1)/2:(img_sizey-1)/2]-ind_y;
    
        
    Gx(:,:,15*(ii-1)+1:15*ii) = squeeze(Bx(indy_trans,indx_trans,ind_z+5,:));
    Gy(:,:,15*(ii-1)+1:15*ii) = squeeze(By(indy_trans,indx_trans,ind_z+5,:));
    %Gx(:,:,15*(ii-1)+1:15*ii) = imtranslate(squeeze(Bx(:,:,ind_z+1,:)),[ind_x,ind_y]);
    %Gy(:,:,15*(ii-1)+1:15*ii) = imtranslate(squeeze(By(:,:,ind_z+1,:)),[ind_x,ind_y]);
end

G = cat(2,Gx,Gy);
G_re = reshape(G,[],15*N);
gamma_new = reshape(gamma_new,[],1);

end