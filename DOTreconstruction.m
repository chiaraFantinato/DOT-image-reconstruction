clearvars
close all
clc

%%

addpath(genpath('iso2mesh-master'))
addpath(genpath('homer2'))
addpath(genpath('S111_Motor'))
load('S111_Motor.nirs','-mat');

nCh = size(SD.MeasList,1)/2;

%% 1 - Plot array configuration

figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
end
title('array configuration')
legend('sources', 'detectors', 'channels')
xlabel('x')
ylabel('y')
zlabel('z')

clear aux iCh src det

%% 2 - Compute SD distance and plot it as histogram

distCh = zeros(nCh,1);
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    distCh(iCh) = sqrt(sum((src-det).^2));
end
figure;
histogram(distCh, 16) 
title('source-detectors distances histogram')
xlabel('[mm]')

clear distCh iCh src det nbins

%% Visualize intensity data at the two wavelengths (not required)

figure;
plot(t,d(:,1:nCh))
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
xlim([t(1) t(end)])
title('intensity data at the first wavelength')

figure;
plot(t,d(:,nCh+1:end))
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
xlim([t(1) t(end)])
title('intensity data at the second wavelength')

%% 3 - Remove noisy channels

% Plot stimuli (not required)
figure;
stem(t,s)
xlabel('Time [s]')
ylabel('A.U.')
title('stimuli')
xlim([t(1) t(end)])
legend('condition 1','condition 2', 'Location', 'southeast')

% Remove from the data (only for these computations) the time points where
% no stimuli were presented
idx_start_first_cond = find(s(:,1),1);
idx_end_first_cond = find(s(:,1),1,'last');
t_end_first_cond = t(idx_end_first_cond) + 15;
idx_end_first_cond = find(t<t_end_first_cond,1,'last');

idx_start_second_cond = find(s(:,2),1);
idx_end_second_cond = find(s(:,2),1,'last');
t_end_second_cond = t(idx_end_second_cond) + 15;
idx_end_second_cond = find(t<t_end_second_cond,1,'last');

d_red = d([idx_start_first_cond:idx_end_first_cond,idx_start_second_cond:idx_end_second_cond],:); % d reduced
dRange = [0.03 3]; % average intensity range
SNRrange = 15;
K = 100;
remCh = removeNoisyChannels(d_red, dRange, SNRrange, K); % 1 = good channel
                                                         % 0 = bad channel
SD.MeasListAct = remCh; 

clear idx_start_first_cond idx_end_first_cond t_end_first_cond
clear idx_start_second_cond idx_end_second_cond t_end_second_cond
clear d_red dRange SNRrange K

% Plot array configuration
figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    if remCh(iCh) == 0 % bad channels 
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'m')
    else % good channels
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    end
end
xlabel('x')
ylabel('y')
zlabel('z')
title('array configuration after having selected bad channels (magenta)')

clear iCh src det

%% Visualize intensity data at the two wavelengths of non-removed channels (not required)

dGood = d(:,remCh == 1);

figure; 
plot(t,dGood(:,1:end/2))
title('intensity data at the first wavelength of the non-removed channels')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')

figure; 
plot(t,dGood(:,end/2+1:end))
title('intensity data at the second wavelength of the non-removed channels')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')

clear dGood

%% Visualize intensity data at the two wavelengths of removed channels (not required)

dBad = d(:,remCh==0);
figure; 
plot(t,dBad(:,1:end/2))
title('intensity data at the first wavelength of the removed channels')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')

figure; 
plot(t,dBad(:,end/2+1:end))
title('intensity data at the second wavelength of the removed channels')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')

clear dBad

%% 4 - Pre-process the fNIRS data

%% 4a - Convert to optical density changes
meanValue = mean(d);
dodConv = -log(abs(d)./meanValue);

clear meanValue

%% Visualize optical density changes at the two wavelengths of non-removed channels (not required)

dodConvGood = dodConv(:,remCh==1);

figure; 
plot(t,dodConvGood(:,1:end/2))
title('optical density chances at the first wavelength of non-removed channels')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

figure; 
plot(t,dodConvGood(:,end/2+1:end))
title('optical density chances at the second wavelength of non-removed channels')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

figure; 
plot(t,dodConvGood(:,1)) 
title('optical density changes, channel 1, wavelength 1')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')
% example of big baseline changes in the optical density signal, which is
% an optimal case for spline interpolation as motion correction technique

clear dodConvGood

%% 4b - Motion correction performed with the spline interpolation function 

% Detect motion artifacts in signal (first step of motion correction with spline intepolation)

% Motion artifacts can be identified using the following parameters: 
% tMotion = 1, tMask = 3, STDThresh = 5 and AmpThresh = 0.5. 

fs = 1/(t(2)-t(1));
tMotion = 1;
tMask = 3;
SDThresh = 5;
AmpThresh = 0.5;

tIncMan = ones(length(t),1); % set it to vectors of ones (this is a vector 
% used to remove manually parts of the data if needed)
 
[tInc,tIncCh] = hmrMotionArtifactByChannel(dodConv, fs, SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);
% tIncCh is a matrix number of samples x twice n of channels which contains 
% for each channel (column) the information about whether an artifact was present (0s) or not (1s).
% tInc is a vector which contains information on whether at that time sample in any of the channel
% was present an artifact (0s) or not (1s). tInc can therefore be obtained
% from tIncCh by setting to 0 every row that contains at least one 0.

%% Visualize detected motion artifacts in good channels at first wavelength (not required)

dodConvGood = dodConv(:,remCh==1);
figure;
plot(t,dodConvGood(:,1:end/2))
hold on;
for i = 1:length(tInc) % here we use tInc since we are displaying all channels and we want to see all time points that have been detected as motion in at least one channel
    if tInc(i) == 0 % For each sample identified as motion artifact
        lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5); % Plot vertical red line
        lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
    end
end
title('wavelength 1')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

clear i lh dodConvGood

%% Spline interpolation

p = 0.99;
dodSpline = hmrMotionCorrectSpline(dodConv,t,SD,tIncCh,p);

clear tMotion tMask SDThresh AmpThresh tIncMan p

%% Compare uncorrected vs. spline corrected data at each good channel with
% superimposed vertical lines for where artifacts were detected (not required)

for iCh = 1:nCh
    if remCh(iCh) == 1 % display only good channels
        figure;
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodSpline(:,iCh))
        title(['spline intepolation, channel ', num2str(iCh), ', wavelength 1'])
        xlim([t(1) t(end)])
        hold on;
        for i = 1:size(tIncCh,1) 
            if tIncCh(i,iCh) == 0 % here we use tIncCh since we are displaying channels one at a time and we are interested in evaluating whether spline was able to correct the artifacts detected specifically in each single channel
                lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5);
                lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
            end
        end
        xlabel('Time [s]')
        ylabel('Optical density [A.U.]')
        legend('uncorrected','corrected')
        pause
        close
    end
end

clear iCh lh i tInc tIncCh 

%% 4c - Band-pass filtering with cut-off frequency 0.01 and 2.5 Hz

lowerCutOff = 0.01;
higherCutOff = 2.5; % when applying the GLM to correct for serial correlation, 
% it is better to leave the entire serial correlation (and therefore physiology) 
% in the signal (e.g., heart beat)

dodFilt = hmrBandpassFilt(dodSpline,fs,lowerCutOff,higherCutOff);

clear lowerCutOff higherCutOff

%% Plot filtered optical density data at one selected good channel (not required)

dodSplineGood = dodSpline(:,remCh == 1);
dodFiltGood = dodFilt(:,remCh == 1);

figure;
plot(t,dodSplineGood(:,1))
hold on;
plot(t,dodFiltGood(:,1))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('band-pass filtering, channel 1, wavelength 1')
legend('unfiltered','filtered')

clear dodSplineGood dodFiltGood

%% 4d - Convert to concentration changes

DPF = [6.24 5.18]; % differential pathlength factor

% Convert to concentration changes
dc = hmrOD2Conc(dodFilt,SD,DPF); % #samples x #chromofore x #channels

clear DPF

%% Plot HbO concentration changes for the non-removed channels (not required)
figure;
plot(t,squeeze(dc(:,1,remCh(1:nCh)==1)))
xlabel('Time [s]')
ylabel('Concentration changes [M]')
xlim([t(1) t(end)])
title('HbO concentration changes of non-removed channels')

%% 5 - Compute the average hemodynamic response via GLM approach

% Compute the average HbO/HbR hemodynamic response across trials for each
% condition using the GLM approach with two different settings:

%% a. Use glmSolveMethod = 2, tRange = [-2 12], Gaussian basis functions
% with steps of 2 s and width of 2 s, driftOrder = 0 and regress SS signals
% choosing for each standard channel the SS with the highest correlation

glmSolveMethod = 2; % least square with correction of serially correlated error
tRange = [-2 12]; % time interval of the HRF (from two seconds before stimulus to 12 after it)
idxBasis = 1; % use as HRF model gaussian functions (the alternative is gamma functions)
paramsBasis = [2 2]; % mean and sd of gaussian function (parameters for gamma functions are dispersion and peak)
rhoSD_ssThresh = 15; % distance in mm lower than which a channel is defined an SS channel (we chose 15 mm looking at the source-detector distances histogram)
flagSSmethod = 1; % for each standard channel the SS channel is chosen as the one with the highest correlation
driftOrder = 0; % the drift order that could be removed is set to 0 because we glmSolveMethod = 2
flagMotionCorrect = 0; % we don't want to do motion correction within GLM (because it has already been performed)
[yavg_regressSS, ~, tHRF, ~, ~, ~, ~, ~, ~] = hmrDeconvHRF_DriftSS(dc,s,t,SD,[],[],tRange,glmSolveMethod,idxBasis,paramsBasis,rhoSD_ssThresh,flagSSmethod,driftOrder,flagMotionCorrect);


%% b. Use glmSolveMethod = 2, tRange = [-2 12], Gaussian basis functions
% with steps of 2 s and width of 2 s, driftOrder = 0 and do not regress SS
% channel signals

glmSolveMethod = 2; % least square with correction of serially correlated error
tRange = [-2 12]; % time interval of the HRF (from two seconds before stimulus to 12 after it)
idxBasis = 1; % use as HRF model gaussian functions (the alternative is gamma functions)
paramsBasis = [2 2]; % mean and sd of gaussian function (parameters for gamma functions are dispersion and peak)
rhoSD_ssThresh = 0; % without SS channels regression
flagSSmethod = 1; % since we are not doing SS regression we don't care about this parameter because it will be not kept into account
driftOrder = 0; % the drift order that could be removed is set to 0 because we glmSolveMethod = 2
flagMotionCorrect = 0; % we don't want to do motion correction within GLM (because it has already been performed)
[yavg, ~, ~, ~, ~, ~, ~, ~, ~] = hmrDeconvHRF_DriftSS(dc,s,t,SD,[],[],tRange,glmSolveMethod,idxBasis,paramsBasis,rhoSD_ssThresh,flagSSmethod,driftOrder,flagMotionCorrect);

clear d dodConv dodSpline dodFilt dc
clear idxBasis paramsBasis rhoSD_ssThresh flagSSmethod driftOrder glmSolveMethod flagMotionCorrect

%% Compare the average hemodynamic response for HbO and HbR computed
% via GLM regressing and not regressing the SS signals at the selected
% channels
     
%chSel = 1:nCh; % to decide which channels select for the comparison
chSel = 3; % selected channels
for iCh = 1:length(chSel)
    figure

    subplot(121)
    plot(tHRF,squeeze(yavg(:,1,chSel(iCh),1)),'r','LineWidth',2)
    hold on;
    plot(tHRF,squeeze(yavg(:,1,chSel(iCh),2)),'b','LineWidth',2)
    plot(tHRF,squeeze(yavg(:,2,chSel(iCh),1)),'m--','LineWidth',2)
    plot(tHRF,squeeze(yavg(:,2,chSel(iCh),2)),'c--','LineWidth',2)
    legend('HbO grasping','HbO squeezing','HbR grasping','HbR squeezing','Location', 'northeast')
    title(['hemodynamic response, GLM without SS regression, channel ', num2str(chSel(iCh))])
    xlabel('Time [s]')
    ylabel('\DeltaHb [M]')
    xlim([tHRF(1) tHRF(end)])
    %ylim([-3e-7 4e-7])
    ylim([-1e-7 2e-7])

    subplot(122)
    plot(tHRF,squeeze(yavg_regressSS(:,1,chSel(iCh),1)),'r','LineWidth',2)
    hold on;
    plot(tHRF,squeeze(yavg_regressSS(:,1,chSel(iCh),2)),'b','LineWidth',2)
    plot(tHRF,squeeze(yavg_regressSS(:,2,chSel(iCh),1)),'m--','LineWidth',2)
    plot(tHRF,squeeze(yavg_regressSS(:,2,chSel(iCh),2)),'c--','LineWidth',2)
    legend('HbO grasping','HbO squeezing','HbR grasping','HbR squeezing', 'Location', 'northeast')
    title(['hemodynamic response, GLM with SS regression, channel ', num2str(chSel(iCh))])
    xlabel('Time [s]')
    ylabel('\DeltaHb [M]')
    xlim([tHRF(1) tHRF(end)])
    %ylim([-3e-7 4e-7])
    ylim([-1e-7 2e-7])

    pause
    close
end

clear chSel iCh

% channel 3: without SS regression we can see activation for grasping while
% without not

%% Convert back the HbO/HbR hemodynamic responses to average optical density hemodynamic responses (end of point 5)

% we need to do this to solve image reconstruction problem later

% the conversion should be performed separately for each condition
yavg_cond1 = yavg(:,:,:,1);
yavg_cond2 = yavg(:,:,:,2); 
yavg_regressSS_cond1 = yavg_regressSS(:,:,:,1);
yavg_regressSS_cond2 = yavg_regressSS(:,:,:,2);

DPF = [6.24 5.18]; % differential pathlength factor

dod(:,:,1) = hmrConc2OD(yavg_cond1,SD,DPF); % #time points x #channels
dod(:,:,2) = hmrConc2OD(yavg_cond2,SD,DPF);
dod_regrSS(:,:,1) = hmrConc2OD(yavg_regressSS_cond1,SD,DPF);
dod_regrSS(:,:,2) = hmrConc2OD(yavg_regressSS_cond2,SD,DPF);

% dod & dod_regrSS: % #time points x #channels x #conditions

clear yavg yavg_cond1 yavg_cond2 yavg_regressSS yavg_regressSS_cond1 yavg_regressSS_cond2 DPF

%% Plot the average optical density hemodynamic responses

dodGood = dod(:,remCh==1); 

figure; 
plot(tHRF,dodGood(:,1:end/2))
title('average optical density hemodynamic responses at the first wavelength of non-removed channels')
xlim([tHRF(1) tHRF(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

figure; 
plot(tHRF,dodGood(:,end/2+1:end))
title('average optical density hemodynamic responses at the second wavelength of non-removed channels')
xlim([tHRF(1) tHRF(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

clear dodGood

%% 6 - Display whole array sensitivity for the first wavelength on GM mesh with all channels

% Load data
load(fullfile('HeadVolumeMesh.mat'))
load(fullfile('GMSurfaceMesh.mat'))
load('S111_Motor.jac','-mat')

% Substitute the tissue type of the nodes with their sensitivity
HeadVolumeMesh.node(:,4) = (sum(J{1}.vol)); 

figure;
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
view([0 90]) % Set the view angle
caxis([-3 0])
colorbar
title('whole array sensitivity on the GM volume mesh, wavelength 1, all channels')

%% 6 - Display whole array sensitivity for the first wavelength on GM mesh without bad channels

% Remove bad channels from Jacobian
JCropped=cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) % For each wavelength (we have two Js)
    tmp = J{i}.vol;
    JCropped{i} = tmp(SD.MeasListAct(SD.MeasList(:,4)==i)==1,:);
end
HeadVolumeMesh_good = HeadVolumeMesh;
HeadVolumeMesh_good.node(:,4) = (sum(JCropped{1}));

figure;
plotmesh(HeadVolumeMesh_good.node,HeadVolumeMesh_good.elem(HeadVolumeMesh_good.elem(:,5)==4,1:4))
view([0 90]) % Set the view angle
caxis([-3 0])
colorbar
title('whole array sensitivity on the GM volume mesh, wavelength 1, non-removed channels')

clear i tmp HeadVolumeMesh_good

%% 7 - Reconstruct HbO/HbR images for both conditions mapped to the GM surface for both GLM settings

% Compute inverse of Jacobian
lambda1 = 0.1;
invJ = cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) % for each Jacobian (so for each of the wavelength)
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    S=svd(JJT);
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda1*max(S)));
end

% Data to reconstruct are optical density changes compared to a baseline.
% In our case the baseline is 0, therefore we want to reconstruct 0-our data
datarecon = -dod;
datarecon_regr = -dod_regrSS;
% dod & dod_regrSS: % #time points x #channels x #conditions

% Inizialize matrices and load useful stuff
nNodeVol = size(HeadVolumeMesh.node,1);  %The node count of the volume mesh
nNodeGM = size(GMSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1); % Number of samples to reconstruct
load('vol2gm.mat')
wavelengths = SD.Lambda; % wavelengths of the system
nWavs = length(wavelengths); % n of wavelengths
nCond = size(s,2); % number of condition

% Initialize final results matrices
hbo.vol = zeros(nFrames,nNodeVol);
hbr.vol = zeros(nFrames,nNodeVol);
hbo.gm = zeros(nFrames,nNodeGM);
hbr.gm = zeros(nFrames,nNodeGM);
hbo_reg.vol = zeros(nFrames,nNodeVol);
hbr_reg.vol = zeros(nFrames,nNodeVol);
hbo_reg.gm = zeros(nFrames,nNodeGM);
hbr_reg.gm = zeros(nFrames,nNodeGM);

% Obtain specific absorption coefficients
Eall = zeros(nWavs,2);
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    Etmp = Etmp(1:2); %HbO and HbR only
    Eall(i,:) = Etmp./1e7; %This will be nWavs x 2;
end

% Reconstruct HbO and HbR images for both condition 1 and 2 mapped to the
% surface GM mesh for GLM setting without SS regression
for cond = 1:nCond
    
    % For each frame
    for frame = 1:nFrames
        
        % Reconstruct absorption changes
        muaImageAll = zeros(nWavs,nNodeVol);
        for wav = 1:nWavs
            dataTmp = squeeze(datarecon(frame,SD.MeasList(:,4)==wav & SD.MeasListAct==1,cond));
            invJtmp = invJ{wav};
            tmp = invJtmp * dataTmp';
            muaImageAll(wav,:) = tmp; %This will be nWavs * nNode
        end
        
        % Convert to concentration changes
        hbo_tmpVol = (Eall(2,2)*muaImageAll(1,:) - Eall(1,2)*muaImageAll(2,:))/(Eall(1,1)*Eall(2,2)-Eall(1,2)*Eall(2,1));
        hbr_tmpVol = (muaImageAll(2,:)-Eall(1,2)*hbo_tmpVol)/Eall(2,2);        
        
        % Map to GM surface mesh
        hbo_tmpGM = (vol2gm*hbo_tmpVol');
        hbr_tmpGM = (vol2gm*hbr_tmpVol');
        
        % Book-keeping and saving
        hbo.vol(frame,:,cond) = hbo_tmpVol;
        hbr.vol(frame,:,cond) = hbr_tmpVol;
        hbo.gm(frame,:,cond) = hbo_tmpGM;
        hbr.gm(frame,:,cond) = hbr_tmpGM;
        
    end
end

% Reconstruct HbO and HbR images for both condition 1 and 2 mapped to the
% surface GM mesh for GLM setting with SS regression
for cond = 1:nCond
    
    % For each frame
    for frame = 1:nFrames
        
        % Reconstruct absorption changes
        muaImageAll = zeros(nWavs,nNodeVol);
        for wav = 1:nWavs
            dataTmp = squeeze(datarecon_regr(frame,SD.MeasList(:,4)==wav & SD.MeasListAct==1,cond));
            invJtmp = invJ{wav};
            tmp = invJtmp * dataTmp';
            muaImageAll(wav,:) = tmp; %This will be nWavs * nNode
        end
        
        % Convert to concentration changes
        hbo_tmpVol = (Eall(2,2)*muaImageAll(1,:) - Eall(1,2)*muaImageAll(2,:))/(Eall(1,1)*Eall(2,2)-Eall(1,2)*Eall(2,1));
        hbr_tmpVol = (muaImageAll(2,:)-Eall(1,2)*hbo_tmpVol)/Eall(2,2);        
        
        % Map to GM surface mesh
        hbo_tmpGM = (vol2gm*hbo_tmpVol');
        hbr_tmpGM = (vol2gm*hbr_tmpVol');
        
        % Book-keeping and saving
        hbo_reg.vol(frame,:,cond) = hbo_tmpVol;
        hbr_reg.vol(frame,:,cond) = hbr_tmpVol;
        hbo_reg.gm(frame,:,cond) = hbo_tmpGM;
        hbr_reg.gm(frame,:,cond) = hbr_tmpGM;
        
    end
end

clear lambda1 invJ i Jtmp JJT S nNodeGM wavelengths Etmp datarecon datarecon_regr Eall nNodeVol
clear cond frame muaImageAll wav dataTmp invJtmp tmp hbo_tmpVol hbr_tmpVol hbo_tmpGM hbr_tmpGM vol2gm

%% Plot the reconstructed images at different time points for each condition (GLM setting without SS regression)
tRecon = [0 4 9];
baseline = abs(tRange(1)); % two seconds of baseline
sRecon = fix(tRecon*fs)+fix(baseline*fs); % Convert to samples
load('greyJet.mat')

%% condition 1
close all
iCond = 1;
for iT = 1:length(sRecon) % for each time point
      
    figure
    
    subplot(121)
    % Assign image to fourth column of node
    GMSurfaceMesh.node(:,4) = hbo.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    caxis([-0.08 0.08]) % Set the limit of the colorbar
    view([0 90]) % Set the view angle
    title(['HbO cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, without SS regression'])
    colormap(greyJet) % set the loaded colormap
    hb = colorbar;
    hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
    axis off % remove axis

    subplot(122)
    % Assign image to fourth column of node
    GMSurfaceMesh.node(:,4) = hbo_reg.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    caxis([-0.08 0.08]) % Set the limit of the colorbar
    %caxis([-0.2 0.2]) % Set the limit of the colorbar
    view([0 90]) % Set the view angle
    title(['HbO cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, with SS regression'])
    colormap(greyJet) % set the loaded colormap
    hb = colorbar;
    hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
    axis off % remove axis
    
    figure

    subplot(121)
    GMSurfaceMesh.node(:,4) = hbr.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    view([0 90])
    caxis([-0.08 0.08])
    title(['HbR cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, without SS regression'])
    colormap(greyJet)
    hb = colorbar;
    hb.Label.String = {'\DeltaHbR [\muM]'};
    axis off

    subplot(122)
    GMSurfaceMesh.node(:,4) = hbr_reg.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    view([0 90])
    caxis([-0.08 0.08])
    title(['HbR cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, with SS regression'])
    colormap(greyJet)
    hb = colorbar;
    hb.Label.String = {'\DeltaHbR [\muM]'};
    axis off

end

%% condition 2
close all
iCond = 2;
for iT = 1:length(sRecon) % for each time point
    
    figure
    
    subplot(121)
    % Assign image to fourth column of node
    GMSurfaceMesh.node(:,4) = hbo.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    caxis([-0.15 0.15]) % Set the limit of the colorbar
    view([0 90]) % Set the view angle
    title(['HbO cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, without SS regression'])
    colormap(greyJet) % set the loaded colormap
    hb = colorbar;
    hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
    axis off % remove axis

    subplot(122)
    % Assign image to fourth column of node
    GMSurfaceMesh.node(:,4) = hbo_reg.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    caxis([-0.15 0.15]) % Set the limit of the colorbar
    view([0 90]) % Set the view angle
    title(['HbO cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, with SS regression'])
    colormap(greyJet) % set the loaded colormap
    hb = colorbar;
    hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
    axis off % remove axis
    
    figure

    subplot(121)
    GMSurfaceMesh.node(:,4) = hbr.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    view([0 90])
    caxis([-0.15 0.15])
    title(['HbR cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, without SS regression'])
    colormap(greyJet)
    hb = colorbar;
    hb.Label.String = {'\DeltaHbR [\muM]'};
    axis off

    subplot(122)
    GMSurfaceMesh.node(:,4) = hbr_reg.gm(sRecon(iT),:,iCond);
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    view([0 90])
    caxis([-0.15 0.15])
    title(['HbR cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s, with SS regression'])
    colormap(greyJet)
    hb = colorbar;
    hb.Label.String = {'\DeltaHbR [\muM]'};
    axis off

end

%%
clear baseline iCond iT hb greyJet

