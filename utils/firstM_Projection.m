function gammaNew = firstM_Projection(gammaold)

gamma = reshape(gammaold,15,[]);
sxx = gamma(1,:)+10^-9; 
syy = gamma(2,:)+10^-9; 
szz = gamma(3,:)+10^-9; 
sxy = gamma(4,:); 
sxz = gamma(5,:); 
syz = gamma(6,:); 
s = sxx+syy+szz; 
sxxdx = gamma(7,:);
syydx = gamma(8,:);
szzdx = gamma(9,:);
sxxdy = gamma(10,:);
syydy = gamma(11,:);
szzdy = gamma(12,:);
sxxdz = gamma(13,:);
syydz = gamma(14,:);
szzdz = gamma(15,:);
if s<0
    gammaNew = zeros(size(gammaold));
    return
end

xx = sxx./s;  xy = sxy./s;
yy = syy./s;  xz = sxz./s;
zz = szz./s;  yz = syz./s;
dx_xx = sxxdx./sxx; dx_yy = syydx./syy; dx_zz = szzdx./szz;
dy_xx = sxxdy./sxx; dy_yy = syydy./syy; dy_zz = szzdy./szz;
dz_xx = sxxdz./sxx; dz_yy = syydz./syy; dz_zz = szzdz./szz;


for ii = 1:length(xx)
    M = [xx(ii), xy(ii), xz(ii); ... .
         xy(ii), yy(ii), yz(ii); ...
         xz(ii), yz(ii), zz(ii)];
[V, D] = eig(M);
mux = real(V(1, 3));
muy = real(V(2, 3));
muz = real(V(3, 3));

pjt_factor = max(sqrt(mux^2+muy^2+muz^2),1);
mux_pjc(ii) = mux/pjt_factor;
muy_pjc(ii) = muy/pjt_factor;
muz_pjc(ii) = muz/pjt_factor;

lambda_pjc(ii) = max(min(1.5 * real(D(3, 3)) - .5,1),0);

end

sxx_pjt = s*(lambda_pjc .* mux_pjc.^2 + (1 - lambda_pjc) / 3);
syy_pjt = s*(lambda_pjc .* muy_pjc.^2 + (1 - lambda_pjc) / 3);
szz_pjt = s*(lambda_pjc .* muz_pjc.^2 + (1 - lambda_pjc) / 3);
sxy_pjt = s*lambda_pjc .* mux_pjc .* muy_pjc;
sxz_pjt = s*lambda_pjc .* mux_pjc .* muz_pjc;
syz_pjt = s*lambda_pjc .* muy_pjc .* muz_pjc;

sxxdx_pjt = sxx_pjt*dx_xx;
syydx_pjt = syy_pjt*dx_yy;
szzdx_pjt = szz_pjt*dx_zz;
sxxdy_pjt = sxx_pjt*dy_xx;
syydy_pjt = syy_pjt*dy_yy;
szzdy_pjt = szz_pjt*dy_zz;
sxxdz_pjt = sxx_pjt*dz_xx;
syydz_pjt = syy_pjt*dz_yy;
szzdz_pjt = szz_pjt*dz_zz;

gammaNew = [sxx_pjt;syy_pjt;szz_pjt;sxy_pjt;sxz_pjt;syz_pjt;...
            sxxdx_pjt;syydx_pjt;szzdx_pjt;...
            sxxdy_pjt;syydy_pjt;szzdy_pjt;...
            sxxdz_pjt;syydz_pjt;szzdz_pjt];

gammaNew = reshape(gammaNew,[],1);
end