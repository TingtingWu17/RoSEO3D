function loc_data = get_loc_data3(gammaf, recovStruct,FPSFx, FPSFy,)
%get_loc_data computes molecular parameters (localization data) consist of
%brightness, position and confidence level
%->---
%input
%->---
%gammaf:    array(3N,n_f) -molecular parameter estimates
%grid_p:    array(N,1)    -lateral grid points along x axis
%T:         array(N,N)    -point-wise confidence map
%thr:       scalar [0,1]  -point-wise confidence threshold
%---->-
%output
%---->-
%loc_data:   array with columns [frame _number  brightness, x,y, confidence]

loc_data = [];

%grid points passing through origin along x
grid_p = recovStruct.lateral_grid_p;

%number of grid points along x axis
grid_length = length(grid_p);

%total number of grid points
n_grid_p = recovStruct.n_grid_p;


%number of frames
n_frames = recovStruct.num_frames;

%loop over frames
for i = 1:n_frames

    %current molecular estimates
    z = gammaf(:, i);

    %molecular parameters related to XX basis
    z_xx = [z(1:n_grid_p); z(6 * n_grid_p + 1:7 * n_grid_p); z(7 * n_grid_p + 1:8 * n_grid_p); z(12 * n_grid_p + 1:13 * n_grid_p); z(13 * n_grid_p + 1:14 * n_grid_p)];

    %molecular parameters related to YY basis
    z_yy = [z(n_grid_p + 1:2 * n_grid_p); z(8 * n_grid_p + 1:9 * n_grid_p); z(9 * n_grid_p + 1:10 * n_grid_p); z(14 * n_grid_p + 1:15 * n_grid_p); z(15 * n_grid_p + 1:16 * n_grid_p)];

    %molecular parameters related to ZZ basis
    z_zz = [z(2 * n_grid_p + 1:3 * n_grid_p); z(10 * n_grid_p + 1:11 * n_grid_p); z(11 * n_grid_p + 1:12 * n_grid_p); z(16 * n_grid_p + 1:17 * n_grid_p); z(17 * n_grid_p + 1:18 * n_grid_p)];

    %molecular parameters related toXY basis
    z_xy = [z(3 * n_grid_p + 1:4 * n_grid_p); zeros(n_grid_p, 1); zeros(n_grid_p, 1)];

    %molecular parameters related to XZ basis
    z_xz = [z(4 * n_grid_p + 1:5 * n_grid_p); zeros(n_grid_p, 1); zeros(n_grid_p, 1)];

    %molecular parameters related to YZ basis
    z_yz = [z(5 * n_grid_p + 1:6 * n_grid_p); zeros(n_grid_p, 1); zeros(n_grid_p, 1)];

    %extract brightness-scaled second moments estimates (XX, YY and ZZ)
    z_secondMoments = reshape(z(1:3 * n_grid_p), sqrt(n_grid_p), sqrt(n_grid_p), 3);

    %checking if current frame contains molecular parameters
    if any(z_secondMoments(:) > 0)


        % get the index position of molecules across XX, YY and ZZ
        [I, J] = find(sum(z_secondMoments, 3) > 0);
        loc_t = [J, I];


        k = 1;
        while k <= size(loc_t, 1)


            pm = [grid_p(loc_t(k, 1)), grid_p(loc_t(k, 2))]; %continuous location
            im = (loc_t(k, 1) - 1) * grid_length + loc_t(k, 2); %index


            % XX basis
            %-----------------------------------------------------


            imx_xx = pm(1, 1) + (z_xx(im + n_grid_p) / (eps + z_xx(im))) * 10^2; %x position
            imy_xx = pm(1, 2) + (z_xx(im + 2 * n_grid_p) / (eps + z_xx(im))) * 10^2; %y position
            br_m_xx = z_xx(im); %brightness scales

            imx_xxR = (z_xx(im + 3 * n_grid_p) / (eps + z_xx(im))) * 10^2; %x shift position
            imy_xxR = (z_xx(im + 4 * n_grid_p) / (eps + z_xx(im))) * 10^2; %y shift position  
            % YY basis
            %-----------------------------------------------------

            imx_yy = pm(1, 1) + (z_yy(im + n_grid_p) / (eps + z_yy(im))) * 10^2; %x position
            imy_yy = pm(1, 2) + (z_yy(im + 2 * n_grid_p) / (eps + z_yy(im))) * 10^2; %y position
            br_m_yy = z_yy(im); %brightness scales
            
            imx_yyR = (z_yy(im + 3 * n_grid_p) / (eps + z_yy(im))) * 10^2; %x shift position
            imy_yyR = (z_yy(im + 4 * n_grid_p) / (eps + z_yy(im))) * 10^2; %y shift position

            % ZZ basis
            %-----------------------------------------------------

            imx_zz = pm(1, 1) + (z_zz(im + n_grid_p) / (eps + z_zz(im))) * 10^2; %x position
            imy_zz = pm(1, 2) + (z_zz(im + 2 * n_grid_p) / (eps + z_zz(im))) * 10^2; %y position
            br_m_zz = z_zz(im); %brightness scales
            
            imx_zzR = (z_zz(im + 3 * n_grid_p) / (eps + z_zz(im))) * 10^2; %x shift position
            imy_zzR = (z_zz(im + 4 * n_grid_p) / (eps + z_zz(im))) * 10^2; %y shift position

            % XY
            %-----------------------------------------------------
            br_m_xy = z_xy(im); %brightness scales

            % XZ
            %-----------------------------------------------------
            br_m_xz = z_xz(im); %brightness scales

            % YZ
            %-----------------------------------------------------
            br_m_yz = z_yz(im); %brightness scales

            % combine position estimates
            %-----------------------------------------------------

            br_m = br_m_xx + br_m_yy + br_m_zz; % molecule brightness is the sum across XX, YY and ZZ


            imx = ((imx_xx * br_m_xx) + (imx_yy * br_m_yy) + (imx_zz * br_m_zz)) / (eps + br_m);
            imy = ((imy_xx * br_m_xx) + (imy_yy * br_m_yy) + (imy_zz * br_m_zz)) / (eps + br_m);
            
            imxR = ((imx_yyR * br_m_yy) + (imx_zzR * br_m_zz)) / (eps + br_m_yy+br_m_zz);
            imyR = ((imy_yyR * br_m_yy) + (imy_zzR * br_m_zz)) / (eps + br_m_yy+br_m_zz);

            % map to sencond moments
            %-----------------------------------------------------

            secondM = [br_m_xx, br_m_yy, br_m_zz, br_m_xy, br_m_xz, br_m_yz] / (br_m + eps);


            %update localizaton data
            %----------------------------------------------------
            loc_data = [loc_data; i, imx, imy, br_m, secondM,imxR,imyR];


            k = k + 1;
        end

    end

    
    loc_data(loc_data(4)<0.2*max(loc_data(4))) = [];
    imxR_mean = mean(loc_data(end-1));
    imyR_mean = mean(loc_data(end));
    
        %y_channel
    FXXy = FPSFy.FXXy;
    FYYy = FPSFy.FYYy;
    FZZy = FPSFy.FZZy;
    FXYy = FPSFy.FXYy;
    FXZy = FPSFy.FXZy;
    FYZy = FPSFy.FYZy;
    %gradients
    FXXydx = FPSFy.FXXydx;
    FXXydy = FPSFy.FXXydy;
    FYYydx = FPSFy.FYYydx;
    FYYydy = FPSFy.FYYydy;
    FZZydx = FPSFy.FZZydx;
    FZZydy = FPSFy.FZZydy;
    
    FXXy = imtranslate(FXXy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FYYy = imtranslate(FYYy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FZZy = imtranslate(FZZy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FXYy = imtranslate(FXYy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FXZy = imtranslate(FXZy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FYZy = imtranslate(FYZy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    
    FXXydx = imtranslate(FXXydx,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FXXydy = imtranslate(FXXydy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FYYydx = imtranslate(FYYydx,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FYYydy = imtranslate(FYYydy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FZZydx = imtranslate(FZZydx,[imxR_mean,imyR_mean]/58.5,'bicubic');
    FZZydy = imtranslate(FZZydy,[imxR_mean,imyR_mean]/58.5,'bicubic');
    
    
    
 
            %y_channel
    FXXx = FPSFx.FXXx;
    FYYx = FPSFx.FYYx;
    FZZx = FPSFx.FZZx;
    FXYx = FPSFx.FXYx;
    FXZx = FPSFx.FXZx;
    FYZx = FPSFx.FYZx;
    %gradients
    FXXydx = FPSFy.FXXydx;
    FXXydy = FPSFy.FXXydy;
    FYYydx = FPSFy.FYYydx;
    FYYydy = FPSFy.FYYydy;
    FZZydx = FPSFy.FZZydx;
    FZZydy = FPSFy.FZZydy;

    %joint x and y channel PSFs
    FXX(:, :, 1) = FXXx;
    FXX(:, :, 2) = FXXy;
    FYY(:, :, 1) = FYYx;
    FYY(:, :, 2) = FYYy;
    FZZ(:, :, 1) = FZZx;
    FZZ(:, :, 2) = FZZy;

    FXXdx(:, :, 1) = FXXxdx;
    FXXdx(:, :, 2) = FXXydx;
    FXXdy(:, :, 1) = FXXxdy;
    FXXdy(:, :, 2) = FXXydy;
    FYYdx(:, :, 1) = FYYxdx;
    FYYdx(:, :, 2) = FYYydx;
    FYYdy(:, :, 1) = FYYxdy;
    FYYdy(:, :, 2) = FYYydy;
    FZZdx(:, :, 1) = FZZxdx;
    FZZdx(:, :, 2) = FZZydx;
    FZZdy(:, :, 1) = FZZxdy;
    FZZdy(:, :, 2) = FZZydy;

    FXY(:, :, 1) = FXYx;
    FXY(:, :, 2) = FXYy;
    FXZ(:, :, 1) = FXZx;
    FXZ(:, :, 2) = FXZy;
    FYZ(:, :, 1) = FYZx;
    FYZ(:, :, 2) = FYZy;

    FXXRdx = FXXdx; FXXRdx(:,:,1) = zeros(size(FXXRdx(:,:,1))); FXXRdx(:,:,2) = zeros(size(FXXRdx(:,:,1)));
    FYYRdx = FYYdx; FYYRdx(:,:,1) = zeros(size(FXXRdx(:,:,1)));
    FZZRdx = FZZdx; FZZRdx(:,:,1) = zeros(size(FXXRdx(:,:,1)));

    FXXRdy = FXXdx; FXXRdy(:,:,1) = zeros(size(FXXRdx(:,:,1))); FXXRdy(:,:,2) = zeros(size(FXXRdx(:,:,1))); 
    FYYRdy = FYYdx; FYYRdy(:,:,1) = zeros(size(FXXRdx(:,:,1)));
    FZZRdy = FZZdx; FZZRdy(:,:,1) = zeros(size(FXXRdx(:,:,1)));
    
    
end