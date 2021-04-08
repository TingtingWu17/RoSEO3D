%% Necessary information

% place following files in a directory called source_dir
% 1- transmissionRatio_name: a mat file that contains transmit_ratio_L2R as
% a variable
% transmit_ratio_L2R is a scalar indicating the light transmission ratio
% between left and right channels
%
% 2- zerothOrder_name: a mat file containing zerothOrder as a variable
% zerothOrder_RL is a 1*2 array specifying the zeroth order factor
% in the pupil plane according to [x_channel y_channel] or [right channel left
% channel] formt
%
% 3- image stack and background stack in binary format
% You may want to write your image and background stack into binary format
% with names img_stack_name.bin and backg_stack_name.bin accordingly
% use writeSMLMbackg2bin function

%% simulate images

[backg,SMLM_img,I,SM_GT,basis_matrix_opt1,basis_matrix_opt2,angle_GT] = generate_estimation_image;
NLL_GT = sum(sum(I+backg-SMLM_img.*log(I+backg+10^-16)));


figure();
subplot(2,1,1); imagesc(I); axis image; title('GT image'); colorbar
subplot(2,1,2); imagesc(SMLM_img); axis image;  title('noise image'); colorbar
%% format data and computational resources

%----------------------------------------------
num_frames = 1400; % total number of frames in your stack
num_frames_backg = 1400; % number of background frames; it can be equal to ONE frame or the same as num_frames
num_frames_p = 100; % number of frames to be analyzed in parallel by all workers
number_of_workers = 4; % number of workers
imgSizeH = 570; 
imgSizeW = 2048;
%imgCenter = [262,369];
%imgCenterX = [302,369]; 
%imgCenterY = [302,369]; 
imgCenterX = [7.21*10^4,1.85*10^4]/58.5;  %in unit of pixel
imgCenterY = ([4.80*10^4,2.05*10^4]+[0,25])/58.5;  %in unit of pixel
imgCenterX(:) = imgCenterX(end:-1:1);
imgCenterY(:) = imgCenterY(end:-1:1);
imgSize = size(SMLM_img,1);
%imgSize = 357;
toPhoton = 1/3.4683;


%% analyze via RoSEO

%----------------------------------------------
% create nanoscope object analysis

maskName = fullfile('phasemask', 'pixOL_v12');
emitter_wavelength = 610; %nm
refractiveIndx=1.515;% objective refractive index
sampleRefractiveIndx=1.314;% sample refractive index

left_to_right_trans_ratio = 1;%1.1449; %transmit_ratio_L2R; for Di03-R488/561,523/610 bandpass; y channel.x channel
%1.5350;%transmit_ratio_L2R; for Di03-R488/561, 593/46; y channel.x channel
zerothOrder = [0,0];%zerothOrder_RL;
%construct phasemaskpara
phasemaskpara.zeroorder = zerothOrder;
phasemaskpara.maskname = maskName;

imgPara.pix_sizex = 58.5; %the griding unit for RoSEO3D, in unit of nm
imgPara.pix_sizez = 50;
imgPara.number_axial_pixels = 24;
imgPara.axial_grid_points = [-4:1:19]*imgPara.pix_sizez;

%%
n1 = Nanoscope('imageSize', imgSize*3+6,...
    'ADcount', 1/3.4683,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);
% create PSF matrix accounting for channel transmission ratio

%
[PSFx, PSFy] = n1.createPSFstruct3D(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'normal_focal_plane',-450,...
    'axial_grid_points',imgPara.axial_grid_points); %distance is is nm
%%
n1 = Nanoscope('imageSize', imgSize,...
    'ADcount', 1/3.4683,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);

[FPSFx, FPSFy] = n1.createPSFstruct(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'normal_focal_plane',-450,...
    'molecule_plane',[700,0]);




%%

Bx = cat(4,PSFx.XXx,PSFx.YYx,PSFx.ZZx,PSFx.XYx,PSFx.XZx,PSFx.YZx,...
     PSFx.XXxdx,PSFx.YYxdx,PSFx.ZZxdx,...
     PSFx.XXxdy,PSFx.YYxdy,PSFx.ZZxdy,...
     PSFx.XXxdz,PSFx.YYxdz,PSFx.ZZxdz);
By = cat(4,PSFy.XXy,PSFy.YYy,PSFy.ZZy,PSFy.XYy,PSFy.XZy,PSFy.YZy,...
     PSFy.XXydx,PSFy.YYydx,PSFy.ZZydx,...
     PSFy.XXydy,PSFy.YYydy,PSFy.ZZydy,...
     PSFy.XXydz,PSFy.YYydz,PSFy.ZZydz);

imgPara.img_size = imgSize;
imgPara.PSF_size = PSF_size_opt;
imgPara.Bx = Bx;
imgPara.By = By;

%%
%initialize a parallel pool
p = gcp('nocreate');
if isempty(p)
    parpool('local', 4)
end
%%
PSFx_isotropic = (1/3*PSFx.XXx(:,:,3)+1/3*PSFx.YYx(:,:,3)+1/3*PSFx.ZZx(:,:,3));
PSFy_isotropic = (1/3*PSFy.XXy(:,:,3)+1/3*PSFy.YYy(:,:,3)+1/3*PSFy.ZZy(:,:,3));
%count1 = 0;
for kk = 1:200
%SMLM_img = poissrnd(I);
%
SM_est = RoSEO3D(n1, SMLM_img, backg,n1,FPSFx, FPSFy,'PSF_size_opt',21);
if abs(SM_est(1))>50 || abs(SM_est(2))>50 || abs(SM_est(3)-350)>30
    aa=1; 
end
% if length(SM_est)==0
%     kk =kk-1;
%     count1 = count1+1;
%     continue
% end
SM_est_save{kk} = SM_est;
end
%%
for kk = 1:200 
    SM_est = SM_est_save{kk};
if length(SM_est)==0
    count1 = count1+1;
    kk = kk-1;
   continue 
end
if abs(SM_est(1,1))>100 || abs(SM_est(1,2))>100
    aaa=1;
    
end
%
saveAngle = [];
saveOren = [];
for ll = 1:size(SM_est,1)
secM(1:6)=SM_est(ll,5:10);
signal = SM_est(ll,4);

[mux,muy,muz,rotMobil] = secondM2SymmCone_for_experiment(double(secM),signal,2,double(basis_matrix_opt2));
if muz<=0
    mux = -mux;
    muy = -muy;
    muz = -muz;
end
saveOren(ll,:) = [mux,muy,muz,rotMobil];
[thetaD, phiD, alphaD] = symmCone2angle(mux,muy,muz,rotMobil);
saveAngle(ll,:) = [thetaD,phiD,alphaD,rotMobil,3*pi-sqrt(rotMobil*8*pi^2+pi^2)];
end

Angle_save{kk} = saveAngle;
%NLL_save{kk} = NLL;

end


%%
count = 0;
detactability = 0;
SM_est_SM1 = [];
SM_est_SM2 = [];
Angle_est_SM1 = [];
Angle_est_SM2 = [];
NLL_est = [];
for kk = 1:200
   SM_est_cur =  SM_est_save{kk};
   Angle_est_cur = Angle_save{kk};
   if size(SM_est_cur,1)==2
       detactability = detactability+1;
       if norm(SM_est_cur(1,1:3)-SM_GT(1,1:3))<norm(SM_est_cur(2,1:3)-SM_GT(1,1:3))
           SM_est_SM1(detactability,:) = SM_est_cur(1,:);
           SM_est_SM2(detactability,:) = SM_est_cur(2,:);
           Angle_est_SM1(detactability,:) = Angle_est_cur(1,:);
           Angle_est_SM2(detactability,:) = Angle_est_cur(2,:);
           
       else
           SM_est_SM1(detactability,:) = SM_est_cur(2,:);
           SM_est_SM2(detactability,:) = SM_est_cur(1,:);
           Angle_est_SM1(detactability,:) = Angle_est_cur(2,:);
           Angle_est_SM2(detactability,:) = Angle_est_cur(1,:);
           %NLL_est(detactability) = NLL_save{kk};
           
       end
   else
       count = count+1;
       SM_outlier{count} = SM_est_cur;
       Angle_outlier{count} = Angle_est_cur;
       %NLL_outlier{count} = NLL_save{kk};
   end
   
end
%%
for kk =1:200
SM_est_SM1(kk,:) = SM_outlier{kk};
SM_est_SM2(kk,:) = SM_outlier{kk};

Angle_est_SM1(kk,:) = Angle_outlier{kk};
Angle_est_SM2(kk,:) = Angle_outlier{kk};

%NLL_est(kk) = NLL_outlier{kk};
end


%%
figure();
hist(NLL_est(NLL_est<1000),100); title('negative loglikelihood'); hold on; plot([100,100],[0,100],'r--','LineWidth',2);
%%
options = {'FaceColor','b',...
            'EdgeColor','none'};
%indx = (SM_est_SM1(:,3)<80&abs(SM_est_SM1(:,1))<50&abs(SM_est_SM1(:,2))<50);
%Angle_est_SM1_phi = Angle_est_SM1(:,2);
%Angle_est_SM1_phi(Angle_est_SM1_phi<0)=Angle_est_SM1_phi(Angle_est_SM1_phi<0)+180;
%indx = abs(Angle_est_SM1_phi-45)<40;
indx = 1:length(SM_est_SM1(:,1));
%indx = SM_est_SM1(:,3)<200&abs(SM_est_SM1(:,1))<50&abs(SM_est_SM1(:,2))<50;
height = 10;
plotSize = [3,2];
figure();
subplot(plotSize(1),plotSize(2),1);
histogram(SM_est_SM1(indx,1),100,options{:});
hold on; plot([SM_GT(1,1),SM_GT(1,1)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(SM_est_SM1(indx,1))),' ','std=',num2str(std(SM_est_SM1(indx,1)))];
title({'x_{est} of SM1';test}); xlabel('x(nm)'); ylabel('count');

subplot(plotSize(1),plotSize(2),2);
histogram(SM_est_SM2(indx,1),100,options{:});
hold on; plot([SM_GT(2,1),SM_GT(2,1)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(SM_est_SM2(indx,1))),' ','std=',num2str(std(SM_est_SM2(indx,1)))];
title({'x_{est} of SM2';test}); xlabel('x(nm)'); ylabel('count');

subplot(plotSize(1),plotSize(2),3);
histogram(SM_est_SM1(indx,2),100,options{:});
hold on; plot([SM_GT(1,2),SM_GT(1,2)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(SM_est_SM1(indx,2))),' ','std=',num2str(std(SM_est_SM1(indx,2)))];
title({'y_{est} of SM1';test}); xlabel('y(nm)'); ylabel('count');

subplot(plotSize(1),plotSize(2),4);
histogram(SM_est_SM2(indx,2),100,options{:});
hold on; plot([SM_GT(2,2),SM_GT(2,2)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(SM_est_SM2(indx,2))),' ','std=',num2str(std(SM_est_SM2(indx,2)))];
title({'y_{est} of SM2';test}); xlabel('y(nm)'); ylabel('count');

subplot(plotSize(1),plotSize(2),5);
histogram(SM_est_SM1(indx,3),100,options{:});
hold on; plot([SM_GT(1,3),SM_GT(1,3)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(SM_est_SM1(indx,3))),' ','std=',num2str(std(SM_est_SM2(indx,3)))];
title({'z_{est} of SM1';test}); xlabel('z(nm)'); ylabel('count');

subplot(plotSize(1),plotSize(2),6);
histogram(SM_est_SM2(indx,3),100,options{:});
hold on; plot([SM_GT(2,3),SM_GT(2,3)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(SM_est_SM2(indx,3))),' ','std=',num2str(std(SM_est_SM2(indx,3)))];
title({'z_{est} of SM2';test}); xlabel('z(nm)'); ylabel('count');
%%

figure();
subplot(plotSize(1),plotSize(2),1);
histogram(Angle_est_SM1(indx,1),100,options{:});
hold on; plot([angle_GT(1,1),angle_GT(1,1)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(Angle_est_SM1(indx,1))),' ','std=',num2str(std(Angle_est_SM1(indx,1)))];
title({'\theta_{est} of SM1';test}); xlabel('\theta(\circ)'); ylabel('count');

subplot(plotSize(1),plotSize(2),2);
histogram(Angle_est_SM2(indx,1),100,options{:});
hold on; plot([angle_GT(2,1),angle_GT(2,1)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(Angle_est_SM2(indx,1))),' ','std=',num2str(std(Angle_est_SM2(indx,1)))];
title({'\theta_{est} of SM2';test}); xlabel('\theta(\circ)'); ylabel('count');

Angle_est_SM1_phi = Angle_est_SM1(indx,2);
Angle_est_SM1_phi(Angle_est_SM1_phi<0)=Angle_est_SM1_phi(Angle_est_SM1_phi<0)+180;
subplot(plotSize(1),plotSize(2),3);
histogram(Angle_est_SM1_phi,100,options{:});
hold on; plot([angle_GT(1,2),angle_GT(1,2)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(Angle_est_SM1_phi)),' ','std=',num2str(std(Angle_est_SM1_phi))];
title({'\phi_{est} of SM1';test}); xlabel('\phi(\circ)'); ylabel('count');

subplot(plotSize(1),plotSize(2),4);
histogram(Angle_est_SM2(indx,2),100,options{:});
hold on; plot([angle_GT(2,2),angle_GT(2,2)],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(Angle_est_SM2(indx,2))),' ','std=',num2str(std(Angle_est_SM2(indx,2)))];
title({'\phi_{est} of SM2';test}); xlabel('\phi(\circ)'); ylabel('count');

subplot(plotSize(1),plotSize(2),5);
histogram(Angle_est_SM1(indx,5),100,options{:});
hold on; plot([0,0],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(Angle_est_SM1(indx,5))),' ','std=',num2str(std(Angle_est_SM1(indx,5)))];
title({'\Omega_{est} of SM1';test}); xlabel('\gamma'); ylabel('count');

subplot(plotSize(1),plotSize(2),6);
histogram(Angle_est_SM2(indx,5),100,options{:});
hold on; plot([0,0],[0,height],'r--','LineWidth',2); 
test = ['mean=',num2str(mean(Angle_est_SM2(indx,5))),' ','std=',num2str(std(Angle_est_SM2(indx,5)))];
title({'\Omega_{est} of SM2';test}); xlabel('\gamma'); ylabel('count');


