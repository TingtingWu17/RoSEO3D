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

%% add folders to search path

%----------------------------------------------
folders = {'utils', 'phasemask'};

for ii = 1:length(folders)

    addpath(genpath(folders{ii}));

end

%% data information
%source_dir = 'F:\data\20210209lipid membrane\offset_subtracted_data\data30\';
source_dir = 'F:\data\20210302-3D lipid membrane\offset_subtracted\data27\'; % directory containing required data
img_stack_name = 'data27_1'; %name of the image stack to be analyzed (add .bin)
backg_stack_name = 'BGEst_lineFitWavelet\data27_1_BKG'; %name of the background stack (add .bin)
%tform_name = '20210209lipid membrane,Di03-R488_561,591_46 geometric transformation data.mat';
%load([source_dir tform_name])


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
imgSize = 41;
%imgSize = 357;
toPhoton = 1/3.4683;

Imgoptions = struct();

Imgoptions.pixSize = 58.5;
Imgoptions.registershifty = 0; %in unit of pixel
Imgoptions.registershiftx = 0;
Imgoptions.imgSize = imgSize;
Imgoptions.imgCenterX = imgCenterX;
Imgoptions.imgCenterY = imgCenterY;
Imgoptions.imgSizeH = imgSizeH;
Imgoptions.imgSizeW = imgSizeW;
Imgoptions.source_dir = source_dir;
Imgoptions.img_stack_name = img_stack_name;
Imgoptions.backg_stack_name = backg_stack_name;
Imgoptions.toPhoton = toPhoton;


%%
%initialize a parallel pool
p = gcp('nocreate');
if isempty(p)
    parpool('local', number_of_workers)
end
num_full_stack = floor(num_frames/num_frames_p);
%%
% map the address of the data to m_img & m_backg



if num_frames_backg == 1
    num_full_stack_backg = 1;
    num_frames_backg_p = 1;
else
    num_full_stack_backg = num_full_stack;
    num_frames_backg_p = num_frames_p;
end


%% analyze via RoSEO

%----------------------------------------------
% create nanoscope object analysis

maskName = fullfile('phasemask', 'pixOL_v12');
emitter_wavelength = 610; %nm
refractiveIndx=1.515;% objective refractive index
sampleRefractiveIndx=1.314;% sample refractive index

left_to_right_trans_ratio = 1.1449; %transmit_ratio_L2R; for Di03-R488/561,523/610 bandpass; y channel.x channel
%1.5350;%transmit_ratio_L2R; for Di03-R488/561, 593/46; y channel.x channel
zerothOrder = [0,0];%zerothOrder_RL;
%construct phasemaskpara
phasemaskpara.zeroorder = zerothOrder;
phasemaskpara.maskname = maskName;


n1 = Nanoscope('imageSize', imgSize,...
    'ADcount', 1/3.4683,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);
% create PSF matrix accounting for channel transmission ratio
%%
[FPSFx, FPSFy] = n1.createPSFstruct3D(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'number_axial_pixels',9,...
    'axial_pixel_size',100,...
    'normal_focal_plane',350);


%% check basis
B = n1.computeBasis_static(n1);
B = B./(sum(B(1:3,:),'all')/3);
%%

loc_data_tot = cell(1, num_full_stack*num_frames_p);


for ll = 1:num_full_stack

    num_char = 0; %initialize the progress bar
    cnt = 1;

    if exist('F', 'var')
        cancel(F);
        clear F
    end
    
    indx_stack = (ll-1)*num_frames_p+1:ll*num_frames_p;
    % load data on demand
    [SMLM_img,backg] = img_pre_processing(Imgoptions,indx_stack);



    %handle for recovery

    if num_frames_backg_p == 1

        RoSEO_h = @(img)RoSEO(n1, img, backg, FPSFx, FPSFy, 'regval', .22);
        parfor k = 1:num_frames_p
            %creat emitter structure

            F(k) = parfeval(RoSEO_h, 2, SMLM_img(:, :, k));

        end
    else
        RoSEO_h = @(img, backg)RoSEO(n1, img, backg, FPSFx, FPSFy, 'regval', 0.6);
        for k = 1:num_frames_p
            %creat emitter structure

            F(k) = parfeval(RoSEO_h, 2, SMLM_img(:, :, k), backg(:, :, k));

        end
    end


    %collect the results and monitor the progress
    loc_data = cell(1, num_frames_p);
    for k = 1:num_frames_p
        % fetchNext blocks until next results are available.
        [completedIndx, gammaf, recovS] = fetchNext(F);
        if ~isempty(gammaf)
            loc_data{completedIndx} = get_loc_data2(gammaf, recovS);
        end
        %free the memory
        clear gammaf;

        %display progress
        num_char = progress_bar(cnt/num_frames_p, num_char, 20);
        cnt = cnt + 1;
    end

    loc_data_tot((ll-1)*num_frames_p+1:ll*num_frames_p) = loc_data;
end

%%

i=71;
temp = loc_data_tot{i+900};
figure(); imagesc(SMLM_img(:,:,i)); axis image; hold on; colormap('hot'); axis off;
scatter(temp(:,2)/58.5-1.5+imgSize/2,temp(:,3)/58.5+imgSize/2-1.5,'w+');
scatter(temp(:,2)/58.5-1.5+imgSize/2+imgSize,temp(:,3)/58.5+imgSize/2-1.5,'w+');

%%
i=77;
figure(); imagesc(SMLM_img(:,1:imgSize,i)); axis image; hold on; colormap('hot');
figure(); imagesc(SMLM_img(:,(1:imgSize)+imgSize,i)); axis image; hold on; colormap('hot');
figure(); imagesc(SMLM_img(:,1:imgSize,i)./max(SMLM_img(:,1:imgSize,i),[],'all')+SMLM_img(:,(1:imgSize)+imgSize,i)./max(SMLM_img(:,(1:imgSize)+imgSize,i),[],'all')); axis image; hold on; colormap('hot');
%%
close all;
%%

% filter localization far from center; and too close
%----------------------------------------------
minimum_distance = 300; % nm
pixelSize = 58.5;
distance_from_center = (imgSize - 50) / 2; %pixels
br_thres = .8e3;
num_frames_p = numel(loc_data_tot);


x_est = cell(1, num_frames_p);
y_est = cell(1, num_frames_p);
secM = cell(1, num_frames_p);
br = cell(1, num_frames_p);
for i = 1:num_frames_p

    if ~isempty(loc_data_tot{i})

        x_est_t = loc_data_tot{i}(:, 2);
        y_est_t = loc_data_tot{i}(:, 3);
        secM_t = loc_data_tot{i}(:, 5:end);
        br_t = loc_data_tot{i}(:, 4);
        %filter based on distance from center
        indx_valid = (x_est_t > -distance_from_center * pixelSize) .* (x_est_t < distance_from_center * pixelSize) .* ...
            (y_est_t > -distance_from_center * pixelSize) .* (y_est_t < distance_from_center * pixelSize);

        x_est_t = x_est_t(indx_valid > 0);
        y_est_t = y_est_t(indx_valid > 0);
        secM_t = secM_t(indx_valid > 0, :);
        br_t = br_t(indx_valid > 0);

        indx = zeros(1, numel(x_est_t));

        %         filter based on intera-distance
        if ~isempty(x_est_t)
            for j = 1:numel(x_est_t)


                dist_t = sqrt(((x_est_t(j) - x_est_t).^2+(y_est_t(j) - y_est_t).^2));
                indx_t = dist_t > 0;
                if ~any(indx_t) || min(dist_t(indx_t > 0)) > minimum_distance

                    if br_t(j) > br_thres
                        indx(j) = 1;
                    end
                end
            end

            x_est{i} = x_est_t(indx > 0);
            y_est{i} = y_est_t(indx > 0);
            secM{i} = secM_t(indx > 0, :);
            br{i} = br_t(indx > 0);
        end
    end
end

%% save results
background_mean = 2.8;
saveOren = [];
saveAngle = [];
for i = 1:num_frames_p
    if ~isempty(loc_data_tot{i})
        loc_data = loc_data_tot{i};
        NSMs = size(loc_data,1);
        for ll=1:NSMs
            
            secM_cur(1:6)=loc_data(ll,5:10);
            signal = loc_data(ll,4);
            [mux,muy,muz,rotMobil] = secondM2SymmCone_for_experiment(secM_cur,loc_data(ll,4),background_mean,B.');
            if muz<=0
                mux = -mux;
                muy = -muy;
                muz = -muz;
            end
            saveOren= cat(1,saveOren,[loc_data(ll,1:end),mux,muy,muz,rotMobil]);
            [thetaD, phiD, alphaD] = symmCone2angle(mux,muy,muz,rotMobil);
            omega = 2*pi*(3/2-sqrt(2*rotMobil+1/4));
            saveAngle = cat(1,saveAngle,[i,loc_data(ll,2:end),thetaD,phiD,alphaD,omega]);

       end
%
 end
end

%%

save(['est','img_stack_name','_imgSize',num2str(imgSize),'_imgCenter',num2str(imgCenter(1)),'_',num2str(imgCenter(2)),'_v1','.mat'])
%%
saveAngle_temp = saveAngle;
saveAngle_temp(saveAngle_temp(:,4)<380,:) = [];
figure();
subplot(2,4,1);
hist(saveAngle_temp(:,end-3));

subplot(2,4,2);
hist(saveAngle_temp(:,end));

subplot(2,4,3)
scatter(saveAngle_temp(:,end-3),saveAngle_temp(:,end),5,'filled','MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);

subplot(2,4,4)
hist(saveAngle(:,4),20);


subplot(2,4,5)
scatter(saveAngle_temp(:,end-3),saveAngle_temp(:,4),5,'filled','MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

subplot(2,4,6)
scatter(saveAngle_temp(:,end),saveAngle_temp(:,4),5,'filled','MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

%%

figure(); hold on;
I_set = [0,380,1000,2000,100000];
color = {'b','r','g','k','y','m','c'};
for ii = 1:length(I_set)-1
    saveAngle_temp = saveAngle;
    saveAngle_temp(saveAngle_temp(:,4)<I_set(ii) ,:) = [];
    saveAngle_temp(saveAngle_temp(:,4)>I_set(ii+1),:) = [];
   

    scatter(rand(length(saveAngle_temp),1)*100+100*ii,saveAngle_temp(:,end-3),20,'filled','MarkerFaceColor',color{ii},'MarkerEdgeColor',color{ii},...
    'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);


end
