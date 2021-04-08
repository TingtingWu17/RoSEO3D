PSFx_isotropic = (1/3*PSFx.XXx(:,:,3)+1/3*PSFx.YYx(:,:,3)+1/3*PSFx.ZZx(:,:,3));
PSFy_isotropic = (1/3*PSFy.XXy(:,:,3)+1/3*PSFy.YYy(:,:,3)+1/3*PSFy.ZZy(:,:,3));
SM_est = RoSEO(n1, SMLM_img, backg, PSFx, PSFy);
 
