%% Figure 2
% Written 05-08-2024 ZBS (zbsinger@mit.edu)

%% Neural

% Figures:
% 2.A Example PSTH (split by contrast and preferred & non preferrerd ori) with raster
% inset for three example jointly recorded units  from example population
% 2.C Linear stimulus sign decoder vs linear choice decoder statistic
% (performance difference above chance vs quality of stimulus decoding)
% 2.D  Linear stimulus sign decoder vs linear confidence decoder statistic
% (correlation between two decoders vs quality of stimulus decoding)
% 2.F Linear vs non-linear decoding choice and confidence 

%% Notes:
% Expt 1: All stimuli & noise seeds with equal probability
% Expt 2: Three stimuli-noise seeds over represented per session

%% Start with clean slate
clearvars -except masterNeuralStruct
close all, clc

% Set Paths
this                = '/Users/zbsinger/MIT Dropbox/Zoe Boundy-Singer/Projects/Conf_phys_public';
thisPath            = user;
structurePath       = strcat([thisPath,'/Structures/Neural']);
functionPath        = strcat([thisPath, '/Functions/']);
addpath([thisPath,structurePath,functionPath])

if ~exist('masterNeuralStruct')  % only load if has been cleared (structure is large)
    load(strcat([structurePath, '/masterNeuralStruct.mat']))
end

% List of recordings:
expt1List       = {'05-28-2021', '06-01-2021', '06-03-2021' , '06-08-2021', '06-10-2021', ...
                   '08-17-2021' , '08-21-2021' ,'08-24-2021', '08-27-2021' , ...
                   '02-05-2023','02-08-2023', '02-11-2023', '02-13-2023', '02-15-2023'};
observerList1   = [repmat({'ziggy'},1,9), repmat({'friedrich'},1,5)];
expt2List       = { '06-15-2021','06-17-2021','06-19-2021','07-01-2021', ... 
                    '09-08-2021','09-15-2021', '09-22-2021','09-25-2021', ... 
                    '04-18-2023','04-20-2023','04-23-2023','04-25-2023','04-27-2023', '05-04-2023', '05-07-2023'};
observerList2   = [repmat({'ziggy'},1,8), repmat({'friedrich'},1,7)];
exptList        = [expt1List expt2List];
exptNum         = [repmat(1,1,numel(expt1List)) repmat(2,1,numel(expt2List))]; 
observerList    = [observerList1 observerList2];
observer{1}     = 'ziggy';
observer{2}     = 'friedrich';
observerInd{1}  = find(strcmp('ziggy',observerList)); % ZIGGY
observerInd{2}  = find(strcmp('friedrich',observerList)); % FRIEDRICH
obsInd          = strcmp('friedrich',observerList);

%% Plotting conventions:
colors(:,:,1)              = [1 0 0; 0 0 1];
colors(:,:,2)              = [1 .75 .75; .75 .75 1];
colors(:,:,3)              = [.1 .1 .1; .1 .1 .1];
colors(:,:,4)              = [.75 .75 .75; .75 .75 .75];
confColor                  = flipud([1 0 0; 0 1 0; 1 1 0; 0 0 1]);

flag.plotPSTH              = 1;

%% 2.A PSTH and raster
exampleDate                = 7;
dataStruct                 = masterNeuralStruct{exampleDate}.dataStructTask; % select data
trialMatrix_index          = masterNeuralStruct{exampleDate}.behavior.trialMatrix.trialMatrix_index; % select behavior data
trialMatrix                = masterNeuralStruct{exampleDate}.behavior.trialMatrix.trialMatrix; % select behavior data

%% Get indices from master-trialMatrix
cont_ind                    = [];
cont_ind(:,1)               = logical(trialMatrix(:,trialMatrix_index.CONTRAST_IND)==1);
cont_ind(:,2)               = logical(trialMatrix(:,trialMatrix_index.CONTRAST_IND)==2);
conf                        = logical(trialMatrix(:,trialMatrix_index.CONFIDENCE)); % 1 = HC ; 0 = LC
side                        = trialMatrix(:,trialMatrix_index.SIDE); %recode S.T 0 = CCW ; 1 = CW
side(side==-1)              = 0;
zeroInd                     = logical(trialMatrix(:,trialMatrix_index.ORIENTATION)==0);

for iR=1:2 % loop through reliability levels (contrast levels)
    uOri{iR} = unique(trialMatrix(:,trialMatrix_index.ORIENTATION));
end

if flag.plotPSTH
    %% PSTH
    bin_size_sec            = 0.020; % 20 ms
    stim_onset_t            = 0;
    edge_0                  = stim_onset_t -.04;
    edge_1                  = stim_onset_t + dataStruct.stimulus_on_duration + .08;
    resp_edges              = edge_0 : bin_size_sec/2 : edge_1;
    num_bins                = length(resp_edges)-1;
    
  
    set(figure(1), 'OuterPosition', [200, 200, 1800, 1200]);
    figSize             = [.3 .3];
    xax                 = linspace(.05,1,4);
    yax                 = xax;
    exampleUnit         = [62   67  38];
    
    for iN = 1:length(exampleUnit)
        
        % Plot example PSTHs split by choice
        psth                = [];
        st                  = [];
        st                  = [];
        spikeTimes          = dataStruct.unit_spike_times(isfinite(dataStruct.stimulus_offsets),exampleUnit(iN)); 
        
        for iR=1:2 % loop through two contrasts
            st{1}{iR}{1}                = spikeTimes(logical(cont_ind(:,iR)),:); % split by high and low contrast
        end
        
        % Extreme CCW ind
        ccwInd                          = trialMatrix(:,trialMatrix_index.ORIENTATION) == uOri{iR}(1) |  trialMatrix(:,trialMatrix_index.ORIENTATION) == uOri{iR}(2);
        st{2}{1}{1}                     = spikeTimes(logical(cont_ind(:,1)) & ccwInd ,:); % extreme CCW ori
        
        % Extreme CCW ind
        cwInd                           = trialMatrix(:,trialMatrix_index.ORIENTATION) == uOri{iR}(10) |  trialMatrix(:,trialMatrix_index.ORIENTATION) == uOri{iR}(11);
        st{2}{1}{2}                     = spikeTimes(logical(cont_ind(:,1))& cwInd ,:); % extreme CW ori
        
        %% Raster:
        figure; box off; axis square; hold on
        counter = 0;
        for iContext = 1:2 % contexts are 1) extreme CCW or 2) extreme CW ori
            for iT =  1:size( st{2}{1}{iContext} ,1)
                unitSpikeTimes = st{2}{1}{iContext}{iT};
                counter=counter+1;
                for iST=1:length(unitSpikeTimes)
                    line([unitSpikeTimes(iST) unitSpikeTimes(iST)],[counter-1 counter],'Color',colors(iContext,:,1),'linewidth',.5)
                end
            end
        end
        line([0   0  ],[0 210],'Color','m','linewidth',1)
        line([.5 .5  ],[0 210],'Color','m','linewidth',1)
        ylim([0 210])
        xlim([-.05 .6])
        xlabel('Time (seconds)', 'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 12)
        ylabel('Trial', 'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 12)
        title(['Unit #=' num2str(exampleUnit(iN))], 'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 12)
        
        %%
        figure(1)
        for iT=1:2 % two types of splits
            for iR=1:size(st{iT},2) % two contrasts (only high contrast used for figure)
                for iS = 1:size(st{iT}{1},2) % type of split
                    unit_st                     = st{iT}{iR}{iS}(:);
                    y                           = cellfun(@isempty, unit_st);
                    unit_st(y)                  = [];
                    valid_trials                = cellfun(@(x) ~isnan(x(1)), unit_st);
                    num_valid_trials            = sum(valid_trials);
                    psth{iT}{iR}{iS}(1,:)       = (1/bin_size_sec) * histcounts(cell2mat(unit_st), resp_edges ) ./ num_valid_trials;
                end
            end
        end
        
        for iT = 1:2 % two types of  splits
            for iR = 1:size(st{iT},2)
                for iS = 1:size(st{iT}{1},2) % type of split determins # of lines
                    axes('position', [xax(iN) yax(iT) figSize]);
                    phyplot([], [], 'k-', 'xticks', [edge_0 edge_1], 'yticks', [0, 80], 'linewidth', 1, 'width', 1, 'fontsize', 12);
                    hold on;
                    if iT==1
                        pl(iS) = plot(resp_edges(1:end-1), nanmean(psth{iT}{iR}{iS},1), '-','color',colors(1,:,iR),'linewidth',3);
                    elseif iT==2
                        pl(iS) = plot(resp_edges(1:end-1), nanmean(psth{iT}{iR}{iS},1), '-','color',colors(iT,:,iS),'linewidth',3);
                    end
                    
                    stim_on_vec     = [0 ];
                    stim_off_vec    = [stim_on_vec+ dataStruct.stimulus_on_duration];
                    
                    for iV=1:length(stim_on_vec)
                        plot([ stim_on_vec(iV) stim_on_vec(iV)],  [0  80],'k--','linewidth',1)
                        plot([ stim_off_vec(iV) stim_off_vec(iV)],  [0  80],'k--','linewidth',1)
                    end
                    
                    title(['Unit = ' num2str(exampleUnit(iN)) '; '],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
                    xlabel('time(s)','fontsize', 14);
                    ylabel('spike rate (ips)','fontsize', 14);
                end
            end
        end
        
        HC_task_theta               = masterNeuralStruct{exampleDate}.dataStructTask.orientations{1};
        LC_task_theta               = masterNeuralStruct{exampleDate}.dataStructTask.orientations{2};
        yaxLim                      = round( max(masterNeuralStruct{exampleDate}.dataStructTask.mu_rates{1}(:,exampleUnit(iN)))) +5;
        
        axes('position', [xax(iN) yax(3) figSize]); % Tuning curve (high vs low contrast)
        axisHandle=phyplot([], [], 'k-', 'xticks', [-20:10:20], 'yticks', [0, yaxLim], 'linewidth', 1, 'width', 1, 'fontsize', 12);
        hold on;
        
        plot( HC_task_theta,masterNeuralStruct{exampleDate}.dataStructTask.mu_rates{1}(:,exampleUnit(iN)),'-','color',colors(1,:,1),'linewidth',1)
        plot( LC_task_theta,masterNeuralStruct{exampleDate}.dataStructTask.mu_rates{2}(:,exampleUnit(iN)),'-','color',colors(1,:,2),'linewidth',1)
        
        plot( HC_task_theta,masterNeuralStruct{exampleDate}.dataStructTask.mu_rates{1}(:,exampleUnit(iN)),'ko','markerfacecolor',colors(1,:,1),'linewidth',1)
        plot( LC_task_theta,masterNeuralStruct{exampleDate}.dataStructTask.mu_rates{2}(:,exampleUnit(iN)),'ko','markerfacecolor',colors(1,:,2),'linewidth',1)
    end
end

%% Pre-saved structres:
NL = load(strcat([structurePath  '/decodeSummary_NL.mat'])); % includes non-linear predictions for score 4 trials (stimlus, choice, confidence) 
L  = load(strcat([structurePath  '/decodeSummary_L.mat'])); % includes linear predictions for score 4 trials (stimlus, choice, confidence) 

%% Fig 2.C/D/F LDA for stimulus 
% Loop through population and get percent correct neural and
% behavior (exclude zero, i.e focus on non-ambiguous trials only)

for iDate = 1:size(L.popStruct,2)
    trialMatrixIndex            = masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix_index;
    trialMatrix                 = [masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix(:,trialMatrixIndex.ORIENTATION) masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix(:,trialMatrixIndex.CONTRAST_IND) masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix(:,trialMatrixIndex.SIDE)  masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix(:,trialMatrixIndex.CONFIDENCE)  masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix(:,trialMatrixIndex.TRIAL) masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix(:,trialMatrixIndex.SCORE) ]; % GT ori; GT contrast;  observed choice; observed confidence
    
    % Behavior:
    obsCorrect                  = sign(masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix(:,masterNeuralStruct{iDate}.behavior.trialMatrix.trialMatrix_index.RESPONSE));
    obsCorrect                  = obsCorrect(L.popStruct{iDate}.trialNum);
    
    zeroInd                     = L.popStruct{iDate}.ori==0;
    
    % Stimulus decoder ability to predict stimulus sign:
    neuralPropCorrect(iDate)    =  mean(L.popStruct{iDate}.predStim(~zeroInd)==(L.popStruct{iDate}.ori(~zeroInd)>0));
    
    % Behavior proprtion correct choices:
    behaviorPropCorrect(iDate)  =  mean(obsCorrect(~zeroInd)>0);
    
    
    % Prop correct behavior based on linear/non-linear decoders: 
    %(includes ambiguous stimuli)
    L_choice(iDate)             = L.popStruct{iDate}.percentCorrectChoice;
    L_conf(iDate)               = L.popStruct{iDate}.percentCorrectConf;
    NL_choice(iDate)            = NL.popStruct{iDate}.percentCorrectChoice;
    NL_conf(iDate)              = NL.popStruct{iDate}.percentCorrectConf;
    
    
    % Calculate behavior-stimulus decoding congruency/ chance/ relationship
    % between stimuus decoding and confidence behavior:
    choiceCongreuncy(iDate)     = mean(L.popStruct{iDate}.obsChoice(~zeroInd)==L.popStruct{iDate}.predStim(~zeroInd));
    chance(iDate)               = behaviorPropCorrect(iDate)*neuralPropCorrect(iDate) + (1- behaviorPropCorrect(iDate))*(1-neuralPropCorrect(iDate) );
    corrValChoice(iDate)        = corr(L.popStruct{iDate}.obsChoice(~zeroInd),L.popStruct{iDate}.predStim(~zeroInd));
    corrVal(iDate)              = corr(L.popStruct{iDate}.obsConf(~zeroInd),L.popStruct{iDate}.predStim(~zeroInd));
end


%% 2.B Choice congreuncy (behavior and stimulus sign decoding)
figure;
hold on; axis square;
plot([.5 1],[0 0])
ylim([-.1 .1])
plot(neuralPropCorrect(obsInd), choiceCongreuncy(obsInd)-chance(obsInd),'o','markerfacecolor','y','markeredgecolor','r','markersize',10); % F
plot(neuralPropCorrect(~obsInd), choiceCongreuncy(~obsInd)-chance(~obsInd),'o','markerfacecolor','b','markeredgecolor','r','markersize',10); % Z
ylabel('Consistancy (stimulus sign and behavior)- chance')
xlabel('Prop correct stimulus decoder')

% 2B stat
[corr_val,p]                         = corr(neuralPropCorrect',(choiceCongreuncy-chance)');
stats.fig2B.stimVchoice_corr_val     = corr_val;
stats.fig2B.stimVchoice_p_val        = p;

[p,h]                                = signrank(choiceCongreuncy-chance);
stats.fig2B.diff_p                   = p;
stats.fig2B.medDiff                  = median(choiceCongreuncy-chance);

%% 2.B (alternate)
figure;
hold on; axis square;
plot([.5 1],[.5 1])
plot(neuralPropCorrect(obsInd), chance(obsInd),'o','markerfacecolor','y','markeredgecolor','r','markersize',10); % F
plot(neuralPropCorrect(~obsInd), chance(~obsInd),'d','markerfacecolor','y','markeredgecolor','r','markersize',10); % Z
ylabel('Consistancy (stimulus sign and behavior)')
xlabel('Prop correct stimulus decoder')

%% 2.C  Confidence behavior correaltion with stimulus sign decoding
figure;
hold on; axis square;
plot([.5 1],[0 0])
plot(neuralPropCorrect(obsInd), corrVal(obsInd),'o','markerfacecolor','y','markeredgecolor','r','markersize',10); % F
plot(neuralPropCorrect(~obsInd), corrVal(~obsInd),'o','markerfacecolor','b','markeredgecolor','r','markersize',10); % F
ylabel('Correlation (predicted sign and behavior (confidence)')
xlabel('Prop correct stimulus decoder')

% Stat 2C
[p,h]                               = signrank(corrVal);
stats.fig2C.medDiff                 = median(corrVal);
stats.fig2C.p_val_signrank          = p;
%% Fig 2.E/F  Choice & confidence (linear vs non-linear)

figure; axis square; hold on
plot(L_choice(obsInd),NL_choice(obsInd),'o','markerfacecolor','k','markeredgecolor','r','markersize',10); % F
plot(L_choice(~obsInd),NL_choice(~obsInd),'d','markerfacecolor','k','markeredgecolor','r','markersize',10); % Z
plot([.5 1],[.5 1])
xlim([.5 1])
ylim([.5 1])
xlabel('Percent correct (linear)', 'fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Percent correct (non-linear)','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
title('Choice decoding','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);

figure;axis square; hold on
plot(L_conf(obsInd),NL_conf(obsInd),'o','markerfacecolor','k','markeredgecolor','r','markersize',10); % F
plot(L_conf(~obsInd),NL_conf(~obsInd),'d','markerfacecolor','k','markeredgecolor','r','markersize',10); % Z
plot([.5 1],[.5 1])
xlim([.5 1])
ylim([.5 1])
xlabel('Percent correct (linear)', 'fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Percent correct (non-linear)','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
title('Confidence decoding','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);

figure; hold on; axis square;
xlim([-0.2 0.2])
ylim([0 13])
edges                   = linspace(min(xlim),max(xlim),15);
binCenter               = edges + (edges(2)-edges(1))/2;
choice_diff             = hist(L_choice-NL_choice,binCenter);
bar(binCenter,choice_diff ,'r')
plot([0 0],[0 13])
plot(repmat(median(L_choice-NL_choice),2,1),[0 13],'b')
title(['Choice difference'])

figure; hold on; axis square;
xlim([-0.2 0.2])
ylim([0 13])
conf_diff               = hist(L_conf-NL_conf,binCenter);
bar(binCenter,conf_diff ,'r')
plot([0 0],[0 13])
plot(repmat(median(L_conf-NL_conf),2,1),[0 13],'b')
title(['Confidence difference'])

% Stats 2E/F
[p_choice,h]                = signrank(L_choice-NL_choice);
[p_conf,h]                  = signrank(L_conf-NL_conf);
[p_conf_choice_L,h]         = signrank(L_choice-NL_conf);
[p_conf_choice_NL,h]        = signrank(NL_choice-NL_conf);

stats.fig2E.p_val_choice_diff   =  p_choice;
stats.fig2F.p_val_conf_diff     =  p_conf;
stats.fig2E.med_choice_diff     =  median(L_choice-NL_choice);
stats.fig2F.med_conf_diff       =  median(L_conf-NL_conf);

stats.fig2E.medLchoice          = median(L_choice);
stats.fig2E.medNLchoice         = median(NL_choice);
stats.fig2F.medLconf            = median(L_conf);
stats.fig2F.medNLconf           = median(NL_conf);

stats.fig2.conf_choice_NL_diff = median(NL_choice-NL_conf);
stats.fig2.conf_choice_L_diff  = median(L_choice-NL_conf);
stats.fig2.conf_choice_NL_p    = p_conf_choice_NL;
stats.fig2.conf_choice_L_p     = p_conf_choice_L;

%% Stats for section "Predicting perceptual decision confidence from V1 activity"/ Figure 2

for iDate=1:29
nUnit(iDate)                        = NL.popStruct{iDate}.nUnit;
end

stats.fig2.unitMin                  =  min(nUnit);
stats.fig2.unitMax                  =  max(nUnit);
stats.fig2.unitAvg                  =  median(nUnit);

% Stimulus decoder performance:
stats.fig2.minStimDecoderCorrect    = min(neuralPropCorrect);
stats.fig2.maxStimDecoderCorrect    = max(neuralPropCorrect);
stats.fig2.avgStimDecoderCorrect    = median(neuralPropCorrect);

% Relationship between stimulus decoding performance and choice decoding
% performance
[corr_val,p]                        = corr(neuralPropCorrect',L_choice');
stats.fig2.stimVchoice_corr_val     = corr_val;
stats.fig2.stimVchoice_p_val        = p;