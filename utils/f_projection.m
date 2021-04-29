function gammaNew =f_projection(gamma,r)

%% projection step-1: limit the jumping distance
gamma = reshape(gamma,15,[]);
N = size(gamma,2);
projected_gamma = zeros(size(gamma));

%rearrange molecular parameters
%------------------------------
xx = gamma(1,:);
yy = gamma(2,:);
zz = gamma(3,:);

xxdx = gamma(7,:);
yydx = gamma(8,:);
zzdx = gamma(9,:);
xxdy = gamma(10,:);
yydy = gamma(11,:);
zzdy = gamma(12,:);
xxdz = gamma(13,:);
yydz = gamma(14,:);
zzdz = gamma(15,:);

%zzdyR = gamma(17*N+1:18*N, :);
projected_gammaxy_t = gamma(4,:);
projected_gammaxz_t = gamma(5,:);
projected_gammayz_t = gamma(6,:);

% x1=gamma(1:N,:);
% x2=gamma(N+1:2*N,:);
% x3=gamma(2*N+1:3*N,:);
% x23=sqrt(x2.^2+x3.^2);

xxjoint = sqrt(xxdx.^2+xxdy.^2+xxdz.^2);
yyjoint = sqrt(yydx.^2+yydy.^2+yydz.^2);
zzjoint = sqrt(zzdx.^2+zzdy.^2+zzdz.^2);


%compute the conditions of the projection operator
%-------------------------------------------------
%**************
%XX
cxx = (xx + r * xxjoint) / (1 + r^2);
cxx1 = xxjoint <= xx * r;
cxx2 = xxjoint <= -xx / r;
cxx3 = xxjoint > abs(xx*r);
cxx_t1 = (1 - cxx2) .* cxx1;
cxx_t2 = (1 - cxx2) .* cxx3 .* cxx;

projected_gammaxx_t = cxx_t1 .* xx + cxx_t2;
projected_gammaxx_tx = cxx_t1 .* xxdx + cxx_t2 .* r .* (xxdx) ./ (eps + xxjoint);
projected_gammaxx_ty = cxx_t1 .* xxdy + cxx_t2 .* r .* xxdy ./ (eps + xxjoint);
projected_gammaxx_tz = cxx_t1 .* xxdz + cxx_t2 .* r .* xxdz ./ (eps + xxjoint);
projected_gammaxx_t = projected_gammaxx_t .* (projected_gammaxx_t > 0);


%YY
cyy = (yy + r * yyjoint) / (1 + r^2);
cyy1 = yyjoint <= yy * r;
cyy2 = yyjoint <= -yy / r;
cyy3 = yyjoint > abs(yy*r);
cyy_t1 = (1 - cyy2) .* cyy1;
cyy_t2 = (1 - cyy2) .* cyy3 .* cyy;

projected_gammayy_t = cyy_t1 .* yy + cyy_t2;
projected_gammayy_tx = cyy_t1 .* yydx + cyy_t2 .* r .* (yydx) ./ (eps + yyjoint);
projected_gammayy_ty = cyy_t1 .* yydy + cyy_t2 .* r .* yydy ./ (eps + yyjoint);
projected_gammayy_tz = cyy_t1 .* yydz + cyy_t2 .* r .* yydz ./ (eps + yyjoint);
projected_gammayy_t = projected_gammayy_t .* (projected_gammayy_t > 0);



%ZZ
czz = (zz + r * zzjoint) / (1 + r^2);
czz1 = zzjoint <= zz * r;
czz2 = zzjoint <= -zz / r;
czz3 = zzjoint > abs(zz*r);
czz_t1 = (1 - czz2) .* czz1;
czz_t2 = (1 - czz2) .* czz3 .* czz;

projected_gammazz_t = czz_t1 .* zz + czz_t2;
projected_gammazz_tx = czz_t1 .* zzdx + czz_t2 .* r .* (zzdx) ./ (eps + zzjoint);
projected_gammazz_ty = czz_t1 .* zzdy + czz_t2 .* r .* zzdy ./ (eps + zzjoint);
projected_gammazz_tz = czz_t1 .* zzdz + czz_t2 .* r .* zzdz ./ (eps + zzjoint);
projected_gammazz_t = projected_gammazz_t .* (projected_gammazz_t > 0);


%rearrange the molecular estimates
%---------------------------------
% projected_gamma=vertcat(projected_gamma_t,....
%     projected_gamma_tx,projected_gamma_ty);


projected_gamma(1,:) = projected_gammaxx_t;
projected_gamma(2,:) = projected_gammayy_t;
projected_gamma(3,:) = projected_gammazz_t;
projected_gamma(4,:) = projected_gammaxy_t;
projected_gamma(5,:) = projected_gammaxz_t;
projected_gamma(6,:) = projected_gammayz_t;
projected_gamma(7,:) = projected_gammaxx_tx;
projected_gamma(8,:) = projected_gammayy_tx;
projected_gamma(9,:) = projected_gammazz_tx;
projected_gamma(10,:) = projected_gammaxx_ty;
projected_gamma(11,:) = projected_gammayy_ty;
projected_gamma(12,:) = projected_gammazz_ty;
projected_gamma(13,:) = projected_gammaxx_tz;
projected_gamma(14,:) = projected_gammayy_tz;
projected_gamma(15,:) = projected_gammazz_tz;

projected_gamma = reshape(projected_gamma,[],1);
%projected_gamma = reshape(gamma,[],1);

%% projection step-2: first moments constraint

gammaold = projected_gamma;
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

sxx_pjt = s.*(lambda_pjc .* mux_pjc.^2 + (1 - lambda_pjc) / 3);
syy_pjt = s.*(lambda_pjc .* muy_pjc.^2 + (1 - lambda_pjc) / 3);
szz_pjt = s.*(lambda_pjc .* muz_pjc.^2 + (1 - lambda_pjc) / 3);
sxy_pjt = s.*lambda_pjc .* mux_pjc .* muy_pjc;
sxz_pjt = s.*lambda_pjc .* mux_pjc .* muz_pjc;
syz_pjt = s.*lambda_pjc .* muy_pjc .* muz_pjc;

sxxdx_pjt = sxx_pjt.*dx_xx;
syydx_pjt = syy_pjt.*dx_yy;
szzdx_pjt = szz_pjt.*dx_zz;
sxxdy_pjt = sxx_pjt.*dy_xx;
syydy_pjt = syy_pjt.*dy_yy;
szzdy_pjt = szz_pjt.*dy_zz;
sxxdz_pjt = sxx_pjt.*dz_xx;
syydz_pjt = syy_pjt.*dz_yy;
szzdz_pjt = szz_pjt.*dz_zz;

gammaNew = [sxx_pjt;syy_pjt;szz_pjt;sxy_pjt;sxz_pjt;syz_pjt;...
            sxxdx_pjt;syydx_pjt;szzdx_pjt;...
            sxxdy_pjt;syydy_pjt;szzdy_pjt;...
            sxxdz_pjt;syydz_pjt;szzdz_pjt];

gammaNew = reshape(gammaNew,[],1);

end