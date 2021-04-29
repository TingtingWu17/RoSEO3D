function [G_re,loc_rec,gamma_new] = update_basisMatrix_w_register(N,gamma, loc_rec,imgPara)

SM_est= get_loc_data2_w_register1(gamma, loc_rec);
loc_rec_new_all = (SM_est(:,1:3)).';
loc_rec_adds = loc_rec_new_all-loc_rec;

img_sizex = imgPara.img_sizex;
img_sizey = imgPara.img_sizey;
pix_sizex = imgPara.pix_sizex;
pix_sizez = imgPara.pix_sizez;
Bx = imgPara.Bx;
By = imgPara.By;
B_size = size(Bx,1);


Gx = zeros(img_sizey,img_sizex,15*N);
Gy = zeros(img_sizey,img_sizex,15*N);
%gamma_init = zeros(27*size(lc_max,1),1);
gamma = reshape(gamma,21,[]);
gamma_new = zeros(21,N);
z_slice = size(Bx,3);
z_min = imgPara.axial_grid_points(1);
z_max = imgPara.axial_grid_points(end);
z_min_indx = z_min/pix_sizez;

for ii = 1:N 
    ind_z = round(min(max(loc_rec_new_all(3,ii),z_min),z_max)/pix_sizez);
    ind_x = max(min(round(loc_rec_new_all(1,ii)/pix_sizex),img_sizex),-img_sizex);
    ind_y = max(min(round(loc_rec_new_all(2,ii)/pix_sizex),img_sizey),-img_sizey);  
    
    ind_z_add = ind_z-loc_rec(3,ii)/pix_sizez;
    ind_x_add = ind_x-loc_rec(1,ii)/pix_sizex;
    ind_y_add = ind_y-loc_rec(2,ii)/pix_sizex;

    loc_rec(1:3,ii) = [ind_x,ind_y,ind_z].'.*[pix_sizex,pix_sizex,pix_sizez].';
    gamma_new(:,ii) = gamma_update_2_loc_rec_w_register([ind_x_add,ind_y_add,ind_z_add].*[pix_sizex,pix_sizex,pix_sizez],gamma(:,ii));
    
    indx_trans = (B_size-1)/2+1+[-(img_sizex-1)/2:(img_sizex-1)/2]-ind_x;
    indy_trans = (B_size-1)/2+1+[-(img_sizey-1)/2:(img_sizey-1)/2]-ind_y;
    
        
    Gx(:,:,21*(ii-1)+1:(21*ii-6)) = squeeze(Bx(indy_trans,indx_trans,ind_z-z_min_indx+1,:));
    Gy(:,:,21*(ii-1)+1:(21*ii-6)) = squeeze(By(indy_trans,indx_trans,ind_z-z_min_indx+1,:));
    Gx(:,:,21*(ii-1)+16:21*ii) = Gx(:,:,21*(ii-1)+7:(21*(ii-1)+12));
    Gy(:,:,21*(ii-1)+16:21*ii) = -Gy(:,:,21*(ii-1)+7:(21*(ii-1)+12));
end

G = cat(2,Gx,Gy);
G_re = reshape(G,[],21*N);
gamma_new = reshape(gamma_new,[],1);

end



function gamma_new = gamma_update_2_loc_rec_w_register(loc_remove,gamma1)


S_scdM = gamma1(1:6);
M_dx = gamma1(7:9);
M_dy = gamma1(10:12);
M_dz = gamma1(13:15);
M_dRxy = gamma1(16:21);

M_dx_new = M_dx-S_scdM(1:3)*loc_remove(1)/100;
M_dy_new = M_dy-S_scdM(1:3)*loc_remove(2)/100;
M_dz_new = M_dz-S_scdM(1:3)*loc_remove(3)/100;

gamma_new = [S_scdM;M_dx_new;M_dy_new;M_dz_new;M_dRxy];

end