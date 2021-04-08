function projected_gamma = f_projection2(gamma,r)

%projected_gamma = gamma;

%projection_cone projects a vector onto a second-order cone
%constraint as well as a support constraint captured by the indx variable.
%->---
%input
%->---
%gamma:             array(3N,n_f) - molecular parameter estimates
%N:             scalar        - number of grid points
%r:             scalar        - constraint parameter associated with the secon-order cone
%n_f:           scalar        - number of frames
%optinal input parameters:
%indx:          logical array(3N,n_f) - indices of grid points with no
%molecules associated with them
%---->-
%output
%---->-
%projected_gamma:  array(3N,n_f) - projected molecular parameter estimates

% gamma = reshape(gamma,15,[]);
% N = size(gamma,2);
% projected_gamma = zeros(size(gamma));
% 
% xx = gamma(1,:); xx(xx<=0)=0;
% yy = gamma(2,:); yy(yy<=0)=0;
% zz = gamma(3,:); zz(zz<=0)=0; s = xx+yy+zz;
% xy = gamma(4,:); xy(xy>s/2)=s(xy>s/2)/2; xy(xy<-s/2)=-s(xy<s/2)/2;
% xz = gamma(5,:); xz(xz>s/2)=s(xz>s/2)/2; xz(xz<-s/2)=s(xz<-s/2)/2;
% yz = gamma(6,:); yz(yz>s/2)=s(yz>s/2)/2; yz(yz<-s/2)=s(yz<s/2)/2;
% xxdx = gamma(7,:); xxdx(xx<=0)=0;
% yydx = gamma(8,:); yydx(yy<=0)=0;
% zzdx = gamma(9,:); zzdx(zz<=0)=0;
% xxdy = gamma(10,:); xxdy(xx<=0)=0;
% yydy = gamma(11,:); yydy(yy<=0)=0;
% zzdy = gamma(12,:); zzdy(zz<=0)=0;
% xxdz = gamma(13,:); xxdz(xx<=0)=0;
% yydz = gamma(14,:); yydz(yy<=0)=0;
% zzdz = gamma(15,:); zzdz(zz<=0)=0;
% 
% projected_gamma(1,:) = xx;
% projected_gamma(2,:) = yy;
% projected_gamma(3,:) = zz;
% projected_gamma(4,:) = xy;
% projected_gamma(5,:) = xz;
% projected_gamma(6,:) = yz;
% projected_gamma(7,:) = xxdx;
% projected_gamma(8,:) = yydx;
% projected_gamma(9,:) = zzdx;
% projected_gamma(10,:) = xxdy;
% projected_gamma(11,:) = yydy;
% projected_gamma(12,:) = zzdy;
% projected_gamma(13,:) = xxdz;
% projected_gamma(14,:) = yydz;
% projected_gamma(15,:) = zzdz;


% %%
% gamma = reshape(gamma,15,[]);
% N = size(gamma,2);
% projected_gamma = zeros(size(gamma));
% 
% 
% 
% %rearrange molecular parameters
% %------------------------------
% xx = gamma(1,:);
% yy = gamma(2,:);
% zz = gamma(3,:);
% 
% xxdx = gamma(7,:);
% yydx = gamma(8,:);
% zzdx = gamma(9,:);
% xxdy = gamma(10,:);
% yydy = gamma(11,:);
% zzdy = gamma(12,:);
% xxdz = gamma(13,:);
% yydz = gamma(14,:);
% zzdz = gamma(15,:);
% 
% %zzdyR = gamma(17*N+1:18*N, :);
% projected_gammaxy_t = gamma(4,:);
% projected_gammaxz_t = gamma(5,:);
% projected_gammayz_t = gamma(6,:);
% 
% % x1=gamma(1:N,:);
% % x2=gamma(N+1:2*N,:);
% % x3=gamma(2*N+1:3*N,:);
% % x23=sqrt(x2.^2+x3.^2);
% 
% xxjoint = sqrt(xxdx.^2+xxdy.^2+xxdz.^2);
% yyjoint = sqrt(yydx.^2+yydy.^2+yydz.^2);
% zzjoint = sqrt(zzdx.^2+zzdy.^2+zzdz.^2);
% 
% 
% %compute the conditions of the projection operator
% %-------------------------------------------------
% %**************
% %XX
% cxx = (xx + r * xxjoint) / (1 + r^2);
% cxx1 = xxjoint <= xx * r;
% cxx2 = xxjoint <= -xx / r;
% cxx3 = xxjoint > abs(xx*r);
% cxx_t1 = (1 - cxx2) .* cxx1;
% cxx_t2 = (1 - cxx2) .* cxx3 .* cxx;
% 
% projected_gammaxx_t = cxx_t1 .* xx + cxx_t2;
% projected_gammaxx_tx = cxx_t1 .* xxdx + cxx_t2 .* r .* (xxdx) ./ (eps + xxjoint);
% projected_gammaxx_ty = cxx_t1 .* xxdy + cxx_t2 .* r .* xxdy ./ (eps + xxjoint);
% projected_gammaxx_tz = cxx_t1 .* xxdz + cxx_t2 .* r .* xxdz ./ (eps + xxjoint);
% projected_gammaxx_t = projected_gammaxx_t .* (projected_gammaxx_t > 0);
% 
% 
% %YY
% cyy = (yy + r * yyjoint) / (1 + r^2);
% cyy1 = yyjoint <= yy * r;
% cyy2 = yyjoint <= -yy / r;
% cyy3 = yyjoint > abs(yy*r);
% cyy_t1 = (1 - cyy2) .* cyy1;
% cyy_t2 = (1 - cyy2) .* cyy3 .* cyy;
% 
% projected_gammayy_t = cyy_t1 .* yy + cyy_t2;
% projected_gammayy_tx = cyy_t1 .* yydx + cyy_t2 .* r .* (yydx) ./ (eps + yyjoint);
% projected_gammayy_ty = cyy_t1 .* yydy + cyy_t2 .* r .* yydy ./ (eps + yyjoint);
% projected_gammayy_tz = cyy_t1 .* yydz + cyy_t2 .* r .* yydz ./ (eps + yyjoint);
% projected_gammayy_t = projected_gammayy_t .* (projected_gammayy_t > 0);
% 
% 
% 
% %ZZ
% czz = (zz + r * zzjoint) / (1 + r^2);
% czz1 = zzjoint <= zz * r;
% czz2 = zzjoint <= -zz / r;
% czz3 = zzjoint > abs(zz*r);
% czz_t1 = (1 - czz2) .* czz1;
% czz_t2 = (1 - czz2) .* czz3 .* czz;
% 
% projected_gammazz_t = czz_t1 .* zz + czz_t2;
% projected_gammazz_tx = czz_t1 .* zzdx + czz_t2 .* r .* (zzdx) ./ (eps + zzjoint);
% projected_gammazz_ty = czz_t1 .* zzdy + czz_t2 .* r .* zzdy ./ (eps + zzjoint);
% projected_gammazz_tz = czz_t1 .* zzdz + czz_t2 .* r .* zzdz ./ (eps + zzjoint);
% projected_gammazz_t = projected_gammazz_t .* (projected_gammazz_t > 0);
% 
% 
% %rearrange the molecular estimates
% %---------------------------------
% % projected_gamma=vertcat(projected_gamma_t,....
% %     projected_gamma_tx,projected_gamma_ty);
% 
% 
% projected_gamma(1,:) = projected_gammaxx_t;
% projected_gamma(2,:) = projected_gammayy_t;
% projected_gamma(3,:) = projected_gammazz_t;
% projected_gamma(4,:) = projected_gammaxy_t;
% projected_gamma(5,:) = projected_gammaxz_t;
% projected_gamma(6,:) = projected_gammayz_t;
% projected_gamma(7,:) = projected_gammaxx_tx;
% projected_gamma(8,:) = projected_gammayy_tx;
% projected_gamma(9,:) = projected_gammazz_tx;
% projected_gamma(10,:) = projected_gammaxx_ty;
% projected_gamma(11,:) = projected_gammayy_ty;
% projected_gamma(12,:) = projected_gammazz_ty;
% projected_gamma(13,:) = projected_gammaxx_tz;
% projected_gamma(14,:) = projected_gammayy_tz;
% projected_gamma(15,:) = projected_gammazz_tz;
% 
% projected_gamma = reshape(projected_gamma,[],1);
% %projected_gamma = reshape(gamma,[],1);
% 
 %%
% projected_gamma = firstM_Projection(projected_gamma);
% %%
% gamma = reshape(projected_gamma,15,[]);
% N = size(gamma,2);
% projected_gamma = zeros(size(gamma));
% 
% % for xx,yy,zz
% xx = gamma(1,:); xx(xx<=0)=0;
% yy = gamma(2,:); yy(yy<=0)=0;
% zz = gamma(3,:); zz(zz<=0)=0; s = xx+yy+zz;
% 
% xxdx = gamma(7,:); xxdx(xx<=0)=0;
% yydx = gamma(8,:); yydx(yy<=0)=0;
% zzdx = gamma(9,:); zzdx(zz<=0)=0;
% xxdy = gamma(10,:); xxdy(xx<=0)=0;
% yydy = gamma(11,:); yydy(yy<=0)=0;
% zzdy = gamma(12,:); zzdy(zz<=0)=0;
% xxdz = gamma(13,:); xxdz(xx<=0)=0;
% yydz = gamma(14,:); yydz(yy<=0)=0;
% zzdz = gamma(15,:); zzdz(zz<=0)=0;
% 
% % for xy
% r = 1/2;
% xy = gamma(4,:); xz = gamma(5,:); yz = gamma(6,:);
% norm = sqrt(xy.^2+xz.^2+yz.^2); 
% ind = norm>s*r;
% xy(ind) = xy(ind)./norm(ind).*r;
% xz(ind) = xz(ind)./norm(ind).*r;
% yz(ind) = yz(ind)./norm(ind).*r;
% 
% factor = (s(ind)+r*norm(ind))/(1+r^2)/(s(ind)+10^-6);
% xx(ind) = xx(ind)*factor;
% yy(ind) = yy(ind)*factor;
% zz(ind) = zz(ind)*factor;
% xxdx(ind) = xxdx(ind)*factor;
% yydx(ind) = yydx(ind)*factor;
% zzdx(ind) = zzdx(ind)*factor;
% xxdy(ind) = xxdy(ind)*factor;
% yydy(ind) = yydy(ind)*factor;
% zzdy(ind) = zzdy(ind)*factor;
% xxdz(ind) = xxdz(ind)*factor;
% yydz(ind) = yydz(ind)*factor;
% zzdz(ind) = zzdz(ind)*factor;
% 
% projected_gamma(1,:) = xx;
% projected_gamma(2,:) = yy;
% projected_gamma(3,:) = zz;
% projected_gamma(4,:) = xy;
% projected_gamma(5,:) = xz;
% projected_gamma(6,:) = yz;
% projected_gamma(7,:) = xxdx;
% projected_gamma(8,:) = yydx;
% projected_gamma(9,:) = zzdx;
% projected_gamma(10,:) = xxdy;
% projected_gamma(11,:) = yydy;
% projected_gamma(12,:) = zzdy;
% projected_gamma(13,:) = xxdz;
% projected_gamma(14,:) = yydz;
% projected_gamma(15,:) = zzdz;
gamma = reshape(projected_gamma,15,[]);
xx = gamma(1,:); xx(xx<=0)=0;
yy = gamma(2,:); yy(yy<=0)=0;
zz = gamma(3,:); zz(zz<=0)=0; s = xx+yy+zz;
xy = gamma(4,:); xy(xy>s/2)=s(xy>s/2)/2; xy(xy<-s/2)=-s(xy<s/2)/2;
xz = gamma(5,:); xz(xz>s/2)=s(xz>s/2)/2; xz(xz<-s/2)=s(xz<-s/2)/2;
yz = gamma(6,:); yz(yz>s/2)=s(yz>s/2)/2; yz(yz<-s/2)=s(yz<-s/2)/2;
xxdx = gamma(7,:); xxdx(xx<=0)=0;
yydx = gamma(8,:); yydx(yy<=0)=0;
zzdx = gamma(9,:); zzdx(zz<=0)=0;
xxdy = gamma(10,:); xxdy(xx<=0)=0;
yydy = gamma(11,:); yydy(yy<=0)=0;
zzdy = gamma(12,:); zzdy(zz<=0)=0;
xxdz = gamma(13,:); xxdz(xx<=0)=0;
yydz = gamma(14,:); yydz(yy<=0)=0;
zzdz = gamma(15,:); zzdz(zz<=0)=0;

projected_gamma(1,:) = xx;
projected_gamma(2,:) = yy;
projected_gamma(3,:) = zz;
projected_gamma(4,:) = xy;
projected_gamma(5,:) = xz;
projected_gamma(6,:) = yz;
projected_gamma(7,:) = xxdx;
projected_gamma(8,:) = yydx;
projected_gamma(9,:) = zzdx;
projected_gamma(10,:) = xxdy;
projected_gamma(11,:) = yydy;
projected_gamma(12,:) = zzdy;
projected_gamma(13,:) = xxdz;
projected_gamma(14,:) = yydz;
projected_gamma(15,:) = zzdz;
end