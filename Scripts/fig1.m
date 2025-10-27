%% Figure 1
% Behavior (human, monkey)
% Written 05-03-2024 ZBS (zbsinger@mit.edu)

%% Behavior

% Figures:
% 1.B Example psychometric function split (high contrast)
% 1.C Slope ratio cross human & session
% 1.D Meta uncertainty cross human & session
% 1.F Example confidence function (high & low contrast)

%% Notes:
%% Details of human and monkey behavioral data structure:
% .crossBlock   contains cross subject summary statistics and full model fit parameters
% .masterBehaviorStruct contains sorted data per subject/date. For humans
%               the order follow the list human.observerList. .trials
%               contains trials in a sorted format with paramter labels in
%               .parmEst.index{n}
% .sessionStruct  contains the full trial matrix  of n compleated trials
%                 with columns giving details about single trials- see
%                 .trialmatrix_description and trialMatrix_index for column meanings and
%                 index. .taskNames describes the data-splits (full, block 1 etc.) and gives
%                 their index.


%% Start with clean slate
clearvars -except masterNeuralStruct 
clc;
close all;

% Set Paths
user                = '/Users/zbsinger/MIT Dropbox/Zoe Boundy-Singer/Projects/V1_Confidence_Physiology';
thisPath            = user;
structurePath       = strcat([thisPath,'/Structures/']);
functionPath        = strcat([thisPath, '/Functions/']);
addpath([thisPath,structurePath,functionPath])

%% Human data:
human = load( strcat([structurePath '/Behavior/Human/masterBehaviorStruct_Human.mat']),'masterBehaviorStruct','crossBlock', 'observerList', 'taskList' ); % fixed issue of fit

uncPsychList            = strcmp('UncertaintyPsych',human.taskList);

%% 
%% Monkey data:
monkey                  = load(strcat([structurePath 'Behavior/Monkey/masterBehaviorStruct.mat'])); % load all monkey behavioral sessions

recordingInd            = monkey.crossBlock.recordingInd;
obsInd                  = monkey.crossBlock.obsInd; % 1 = Friedrich; 0 = Ziggy;
recInd                  = find(recordingInd); % 1 = recording ; 0 = behavioral session
nDate                   = sum(recordingInd);

%% Human values:
iT      = 1; % iT 1 pulls out full model fit on all data/ all compleated trials
for iO = 1:length(human.observerList)
    
    % Stimulus strength analysis:
    HCtrials                    = human.masterBehaviorStruct{iO}.trials{1}{iT}{1}; % High contrast
    uOri                        = unique(HCtrials(:,2));
    
    if length(uOri)==11
        largeOriInd             = HCtrials(:,2) == uOri(1) | HCtrials(:,2) == uOri(2) | HCtrials(:,2) == uOri(10) | HCtrials(:,2) == uOri(11);
        smallOriInd             = HCtrials(:,2) == uOri(4) | HCtrials(:,2) == uOri(5) | HCtrials(:,2) == uOri(7) | HCtrials(:,2) == uOri(8);
        
        percentHC_largeOri(iO)  = mean(abs(HCtrials(largeOriInd,3))==2); % Column 3 is signed confdeince; 2 = HC
        percentHC_smallOri(iO)  = mean(abs(HCtrials(smallOriInd,3))==2);
        
        % Low-contrast
        LCtrials                 = human.masterBehaviorStruct{iO}.trials{1}{iT}{2};
        uOriLC                   = unique(LCtrials(:,2));
        nOri(iO)                 = numel(uOriLC);
        
        largeOriInd_LC           = LCtrials(:,2) == uOriLC(1) | LCtrials(:,2) == uOriLC(2) | LCtrials(:,2) == uOriLC(10) | LCtrials(:,2) == uOriLC(11);
        smallOriInd_LC           = LCtrials(:,2) == uOriLC(4) | LCtrials(:,2) == uOriLC(5) | LCtrials(:,2) == uOriLC(7) | LCtrials(:,2) == uOriLC(8);
        
        % stim strength analysis for LC
        percentHC_largeOriLC(iO)  = mean(abs(LCtrials(largeOriInd_LC,3))==2); % Column 3 is signed confidence; 2 = HC
        percentHC_smallOriLC(iO)  = mean(abs(LCtrials(smallOriInd_LC,3))==2);
        
    else
        oriQuart                = round(linspace(1,length(uOri),5));
        largeOriInd             = HCtrials(:,2)<uOri(oriQuart(2)) | HCtrials(:,2)>uOri(oriQuart(4)) ;
        smallOriInd             = (HCtrials(:,2)>uOri(oriQuart(2))& HCtrials(:,2)<uOri(oriQuart(4))) ;
        
        percentHC_largeOri(iO)  = mean(abs(HCtrials(largeOriInd,3))==2); % Column 3 is signed confdeince; 2 = HC
        percentHC_smallOri(iO)  = mean(abs(HCtrials(smallOriInd,3))==2);
        oriQuartLC                  = round(linspace(1,length(uOriLC),5));
        largeOriIndLC           = LCtrials(:,2)<uOriLC(oriQuartLC(2)) | LCtrials(:,2)>uOriLC(oriQuartLC(4)) ;
        smallOriIndLC           = (LCtrials(:,2)>uOriLC(oriQuartLC(2))& LCtrials(:,2)<uOriLC(oriQuartLC(4))) ;
        
        percentHC_largeOriLC(iO)  = mean(abs(LCtrials(largeOriInd_LC,3))==2); % Column 3 is signed confidence; 2 = HC
        percentHC_smallOriLC(iO)  = mean(abs(LCtrials(smallOriInd_LC,3))==2);
        
    end
    
    split_HC(iO,1) = human.masterBehaviorStruct{iO}.paramEst.split{1}{1}{1}(2); % High conf - high cont;
    split_LC(iO,1) = human.masterBehaviorStruct{iO}.paramEst.split{1}{1}{1}(3); % Low conf  - high cont;
    split_HC(iO,2) = human.masterBehaviorStruct{iO}.paramEst.split{1}{1}{2}(2); % High conf - low cont;
    split_LC(iO,2) = human.masterBehaviorStruct{iO}.paramEst.split{1}{1}{2}(3); % Low cont - low cont;
    
    human.uncRatio_highCont(iO)     = log10( 1./(split_LC(iO,1)./split_HC(iO,1))) ; % ratio
    human.uncRatio_lowCont(iO)      = log10(1./(split_LC(iO,2)./split_HC(iO,2))) ;

    human.HC_sens(iO) = human.masterBehaviorStruct{iO}.paramEst.values{1}{iT}(3);
    human.LC_sens(iO) = human.masterBehaviorStruct{iO}.paramEst.values{1}{iT}(4);
    
    human.uncRatio_contrast(iO) = human.LC_sens(iO)./human.HC_sens(iO);
 
    
    %% Meta-uncertainty and sensitivity per block
    if size(human.masterBehaviorStruct{iO}.paramEst.values{1},2)>3 % If subject didn't have 1st block thrown out
        human_meta_unc_block1(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{2}(2);% block 1
        human_meta_unc_block2(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{3}(2); % block 2
        human_meta_unc_block3(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{4}(2); % block 3
        
        human_highCont_sens_block1(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{2}(3);% block 1
        human_highCont_sens_block3(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{4}(3);% block 3
        
    else
        human_meta_unc_block1(iO)       = nan;
        human_meta_unc_block2(iO)       = nan;
        human_meta_unc_block3(iO)       = nan;
        
        human_highCont_sens_block1(iO)  = nan;
        human_highCont_sens_block2(iO)  = nan;
    end
end

%% 1.B/F  Example behavior
% 1B = Z210817 (6)
% 1F = Z210601 (2)

for exampleDate=[6 2]
    iT                  = 1; % iT 1 pulls out full model fit on all data/ all compleated trials
    sessionStruct       = monkey.masterBehaviorStruct{recInd(exampleDate)}.sessionStruct;
    paramEst            = monkey.masterBehaviorStruct{recInd(exampleDate)}.paramEst;
    trials              = monkey.masterBehaviorStruct{recInd(exampleDate)}.trials;
    nRel                = size(trials{1}{iT},2);
    asymFlag            = 1;
    % Save parameter names: 
    index_num                                  = 1;
    
    paramEst.index{iT}.LAPSE_RATE                 = index_num; index_num = index_num + 1;
    paramEst.index{iT}.META_UNC                   = index_num; index_num = index_num + 1;
    paramEst.index{iT}.STIM_SENS                  = index_num:index_num+nRel-1; index_num = index_num+nRel;
    paramEst.index{iT}.STIM_CRIT                  = index_num:index_num+nRel-1; index_num = index_num+nRel;
    paramEst.index{iT}.CONF_CRIT                  = index_num:index_num+asymFlag;
   
    % Intialize plot & set ploting conditions:
    
    set(figure(exampleDate+100), 'OuterPosition', [100 100 1500 1000])
    colors(:,:,1) = [0 0 255; 153 204 255;]/255;
    colors(:,:,2) = [255 0 0; 255 100 100;]/255;
    
    %% Step 1: unpack data
    nRel          = numel(trials{1}{1});
    relList       = 1:nRel;
    nConfCrit     = 1;
    asymFlag      = 1;
    if asymFlag == 1
        nParams       = 2 + nRel + nRel + nConfCrit+1; % [Guess rate, meta-uncertainty], [stimulus sensitivity], [stimulus criterion], [confidence criteria]
    else
        nParams       = 2 + nRel + nRel + nConfCrit;
    end
    iO              = 1; % one observer
    
    for iR=1:nRel
        % Set experiment parameters
        stimCat{iR}   = trials{iO}{iT}{iR}(:,1);
        stimValue{iR} = trials{iO}{iT}{iR}(:,2);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
        choice        = trials{iO}{iT}{iR}(:,3);
        
        if iR ==1 && exampleDate == 6 % 1B
            stats.fig1B.nTrials     = length(trials{iO}{iT}{iR});
        elseif iR ==1 && exampleDate == 2 % 1F
            stats.fig1F.nTrials     = length(trials{iO}{iT}{1})+length(trials{iO}{iT}{2});
        end
        
        % Recode responses
        respOptions                         = [-(nConfCrit+1):-1, 1:(nConfCrit+1)];
        choiceVec{iR}                       = zeros(2*(nConfCrit+1), numel(stimValue{iR}));
        
        for iC = 1:numel(respOptions)
            selTrials                       = (choice == respOptions(iC));
            choiceVec{iR}(iC, selTrials)    = 1;
        end
        
        oriPF{iR}{iT}           = unique(stimValue{iR});
        nPoints{iR}{iT}         = numel(oriPF{iR}{iT});
        
        confidence                                                      = zeros(length(choice),1);
        confidence(sum(choiceVec{iR}([1 4],:))==1)                      = 1; % 1 = HC
        
        %%
        for iP = 1:nPoints{iR}{iT}
            oriInd                          = stimValue{iR} == oriPF{iR}{iT}(iP);
            nTrials{iR}{iT}(iP)             = sum(oriInd);
            obsPF{iR}{iT}(iP)               = mean(choice(oriInd) > 0);
            obsCF{iR}{iT}(iP)               = mean(abs(choice(oriInd))>1);
            
            % Psychometric for low and high confidence trials
            indHighConf                     = abs(choice) >= 2;
            obsPF_HC{iR}{iT}(iP)            = mean(choice(oriInd & indHighConf) > 0);
            obsPF_LC{iR}{iT}(iP)            = mean(choice(oriInd & ~indHighConf) > 0);
            nTrials_HC{iR}{iT}(iP)          = sum(oriInd & indHighConf);
            nTrials_LC{iR}{iT}(iP)          = sum(oriInd & ~indHighConf);
            
        end
    end
    
    %% Step 3: Plot model
    % First select appropriate subset of parameters
    % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta noise, conf criteria]
    %%
    xax         = linspace(0,1,5);
    yax         = linspace(0,1,4);
    sizeFig     = [.3 .3];
    
    for iR=1:nRel
        
        stimPlot                         = -20:.2:20;% -25:.2:25;
        params                           = paramEst.values{iO}{iT}([1, 2+iR, 2+nRel+iR, 2, 3+(2*nRel):end]);
        fitLlh{iR}{iT}                   = getLlhChoice(stimPlot, params, 100,asymFlag);
        paramEst.predPFplot{iR}{iT}      = sum(fitLlh{iR}{iT}(size(fitLlh{iR}{iT}, 1)/2+1:end, :));
        paramEst.predCFplot{iR}{iT}      = getLlhChoice(stimPlot, params, 100,asymFlag)' * [1 0 0 1]';
        paramEst.predXaxis{iR}{iT}       = stimPlot;
        
        % Get split PF by confidence
        params       = paramEst.split{1}{1}{iR};
        noiseInt_HC  = params(2);
        noiseInt_LC  = params(3);
        lapseRate    = params(1);
        stimCrit     = params(4);
        if iR==1
            stat.exampleHC_slope=  noiseInt_HC ;
            stat.exampleLC_slope=  noiseInt_LC ;
        end
        
        propCW_H = lapseRate + (1 - 2*lapseRate)*normcdf((stimPlot - stimCrit)/(2*noiseInt_HC));
        propCW_L = lapseRate + (1 - 2*lapseRate)*normcdf((stimPlot - stimCrit)/(2*noiseInt_LC));

        %% Plot panel 1:
        xaxis           = [-20:5:20];%[-20:10:20];
        contrast_ind    = sessionStruct.relList(iR,1);
        dispersion_ind  = sessionStruct.relList(iR,2);
        
        figure(exampleDate+100);
        axes('position', [xax(1) yax(3) sizeFig])
        phyplot([], [], 'k-', 'xticks', [xaxis], 'yticks', [0:.5:1], 'linewidth', 1, 'width', 1, 'fontsize', 20);
        hold on
        
        plot(stimPlot, paramEst.predPFplot{iR}{iT}, '-', 'linewidth', 2, 'color',colors(contrast_ind,:,dispersion_ind))
        xlabel('Stimulus value','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 25)
        ylabel('Proportion CW','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 25)
        
        markerSize = round(nTrials{iR}{iT}./(max(nTrials{1}{iT}))*15)+1;
        for iP = 1:nPoints{iR}{iT}
            plot(oriPF{iR}{iT}(iP), obsPF{iR}{iT}(iP), 'ko', 'markerfacecolor',colors(contrast_ind,:,dispersion_ind), 'markersize', markerSize(iP));%max([round(nTrials{iR}(iP)/2),1]));
        end
        %% Plot panel 2:
        axes('position', [xax(1) yax(2) sizeFig])
        phyplot([], [], 'k-', 'xticks', [xaxis], 'yticks', [0:.5:1], 'linewidth', 1, 'width', 1, 'fontsize', 20);
        hold on
        
        plot(stimPlot, paramEst.predCFplot{iR}{iT}, '-', 'linewidth', 2, 'color', colors(contrast_ind,:,dispersion_ind))
        xlabel('Stimulus value','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 25)
        ylabel('Proportion HC','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 25)
        
        for iP = 1:nPoints{iR}{iT}
            plot(oriPF{iR}{iT}(iP), obsCF{iR}{iT}(iP), 'ko', 'markerfacecolor', colors(contrast_ind,:,dispersion_ind), 'markersize', markerSize(iP));%max([round(nTrials{iR}(iP)/2),1]));
        end
        
        %% Plot panel 3: text
        axes('position', [xax(1)+.1 yax(1) sizeFig]);axis off;
        
        titleText       = [monkey.masterBehaviorStruct{recInd(exampleDate)}.date];
        text1           = ['Meta-uncertainty: ' num2str(paramEst.values{iO}{iT}(paramEst.index{iT}.META_UNC))];
        text2           = ['Confidence criterion: ' num2str(paramEst.values{iO}{iT}(paramEst.index{iT}.CONF_CRIT(1)))];
        text(0, .85,   titleText, 'fontsize', 25);
        text(0, .85-.15, text1, 'fontsize', 20);
        text(0, .55-.15, text2, 'fontsize', 20);
        
        if  numel(paramEst.index{iT}.CONF_CRIT) ==2
            text3           = ['Confidence criterion (2): ' num2str(paramEst.values{iO}{iT}(end))];
            text(0, .4-.15, text3, 'fontsize', 20);
        end
        
        %% Plot panel 3: Plot PF split out for high and low confidence trials
        xaxis           = [-15:15/2:15];
        axes('position', [xax(iR+2) yax(3) sizeFig])
        phyplot([], [], 'k-', 'xticks', [xaxis], 'yticks', [0:.5:1], 'linewidth', 1, 'width', 1, 'fontsize', 20);
        hold on
        plot(stimPlot, propCW_H, '-', 'linewidth', 2, 'color', [0 1 0])
        plot(stimPlot, propCW_L, '-', 'linewidth', 2, 'color', [1 0 0])
        
        xlabel('Stimulus value','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 25)
        ylabel('Proportion CW','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 25)
        
        markerSize_HC = round(nTrials_HC{iR}{iT}./nTrials{iR}{iT} *10)+2;
        markerSize_LC = round(nTrials_LC{iR}{iT}./nTrials{iR}{iT} *10)+2;
        
        for iP = 1:nPoints{iR}{iT}
            if isnan(obsPF_HC{iR}{iT}(iP))
                obsPF_HC{iR}{iT}(iP)=0;
            end
            leg(2) = plot(oriPF{iR}{iT}(iP), obsPF_LC{iR}{iT}(iP), 'ko', 'markerfacecolor', [1 0 0], 'markersize', markerSize_LC(iP));%max([round(nTrials_LC{iR}(iP)),1]));
            
            leg(1) = plot(oriPF{iR}{iT}(iP), obsPF_HC{iR}{iT}(iP), 'ko', 'markerfacecolor', [0 1 0], 'markersize', markerSize_HC(iP)); %max([round(nTrials_HC{iR}(iP)),1]));
        end
        if iR==1
            legend(leg,{'HC trials', 'LC trials'})%, 'location', 'NorthEastOutside')%,'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 25)
        end
        title({'Split diff = ' num2str(params(3) -params(2))},'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
    end
end


%% 1.C Confidence split ratio
% Group low and high contrast confidence split difference
uncRatio_highCont        = log10(1./(monkey.crossBlock.split_LC(:,1)./monkey.crossBlock.split_HC(:,1))); % LC - HC for High contrast stim
uncRatio_lowCont         = log10(1./(monkey.crossBlock.split_LC(:,2)./monkey.crossBlock.split_HC(:,2))); % LC - HC for low contrast stim

uncRatio_F               = [uncRatio_highCont(obsInd&recordingInd) ]; % F recording sessions
uncRatio_Z               = [uncRatio_highCont(~obsInd&recordingInd)]; % Z recording sessions
uncRatio_H               = [ human.uncRatio_highCont]; % Human


% Calculate percentiles:
h_uncRatio_p = prctile(uncRatio_H,[25 50 75]);
z_uncRatio_p = prctile(uncRatio_Z,[25 50 75]);
f_uncRatio_p = prctile(uncRatio_F,[25 50 75]);

figure;hold on; axis square
xlim([0 4])
ylim([-4 1])
plot(ones(length(uncRatio_H),1)+(rand(length(uncRatio_H),1)-.5)/5,uncRatio_H,'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length(uncRatio_Z),1)+(rand(length(uncRatio_Z),1)-.5)/5,uncRatio_Z,'o','markerfacecolor','r','markeredgecolor','k','markersize',10) % Ziggy
plot(2+ones(length(uncRatio_F),1)+(rand(length(uncRatio_F),1)-.5)/5,uncRatio_F,'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Friedrich

plot([1 1],h_uncRatio_p([1 3]))
plot([1 1]+1,z_uncRatio_p([1 3]))
plot([1 1]+2,f_uncRatio_p([1 3]))

plot([1 ],h_uncRatio_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,z_uncRatio_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,f_uncRatio_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Observer','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Ratio High to Low contrast','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({'Human',  'Ziggy', 'Friedrich'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Log 10 slope ratio','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)


% Stats
[p_f,h_f]           = signrank(uncRatio_F); % recheck number in paper 
[p_z,h_z]           = signrank(uncRatio_Z);
[p_h,h_h]           = signrank(uncRatio_H); % because took log- should be diff from zero

stats.fig1C.diffOne_f  = p_f;
stats.fig1C.diffOne_z  = p_z;
stats.fig1C.diffOne_h  = p_h;


stats.fig1C.uncRatio_f = 10.^median(uncRatio_F);
stats.fig1C.uncRatio_z = 10.^median(uncRatio_Z);
stats.fig1C.uncRatio_h = 10.^median(uncRatio_H);

[p_fh,h_f]           = ranksum(10.^uncRatio_F,10.^uncRatio_H);
[p_zh,h_z]           = ranksum(10.^uncRatio_Z,10.^uncRatio_H);


stats.fig1C.fh       = p_fh;
stats.fig1C.zh       = p_zh;

stats.fig1B.splitRatio = 10.^uncRatio_highCont(exampleDate);
%% 1.D Meta-uncertainty 

human_MU_b1 = human_meta_unc_block1(isfinite(human_meta_unc_block1));
human_MU_b3 = human_meta_unc_block3(isfinite(human_meta_unc_block3));

MU_F = monkey.crossBlock.meta_unc(logical(obsInd)); % meta uncertainty Freddy
MU_Z = monkey.crossBlock.meta_unc(~logical(obsInd)); % meta uncertainty Ziggy

% Calculate percentiles:
h_m1_p = prctile(human_MU_b1,[25 50 75]);
h_m3_p = prctile(human_MU_b3,[25 50 75]);

z_m_p = prctile(MU_Z,[25 50 75]);
f_m_p = prctile(MU_F,[25 50 75]);

figure;hold on; axis square
xlim([0 5])
for iO=1:length(human_MU_b1)
human_rand_vec = [(rand(length(human_MU_b1),1)-.5)/5 (rand(length(human_MU_b1),1)-.5)/5];
plot(ones(length(human_MU_b1),1)+human_rand_vec(iO,1),human_MU_b1(iO),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length(human_MU_b3),1)+human_rand_vec(iO,2),human_MU_b3(iO),'o','markerfacecolor','m','markeredgecolor','k','markersize',10)
plot([1+human_rand_vec(iO,1) 2+human_rand_vec(iO,2)],[human_MU_b1(iO) human_MU_b3(iO)],'k-')
end

plot(2+ones(length(MU_Z),1)+(rand(length(MU_Z),1)-.5)/5,MU_Z,'o','markerfacecolor','r','markeredgecolor','k','markersize',10) % Ziggy
plot(3+ones(length(MU_F),1)+(rand(length(MU_F),1)-.5)/5,MU_F,'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Friedrich

plot([1 1],h_m1_p([1 3]))
plot([1 1]+1,h_m3_p([1 3]))

plot([1 1]+2,z_m_p([1 3]))
plot([1 1]+3,f_m_p([1 3]))

plot([1 ],h_m1_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,h_m3_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

plot([1 ]+2,z_m_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+3,f_m_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Observer','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Meta-uncertainty','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3 4])
xticklabels({'Human Block 1', 'Human Block 3', 'Ziggy', 'Friedrich'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Meta-uncertainty','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

% Stats
[p_m,h_m]               = ranksum(MU_F,MU_Z);
[p_hm_b1,h_hmb1]        = ranksum([MU_F ,MU_Z], human_MU_b1);
[p_hm_b3,h_hm_b3]       = ranksum([MU_F ,MU_Z], human_MU_b3);
[p_hz_b3]               = ranksum([MU_Z], human_MU_b3);
[p_hf_b3]               = ranksum([MU_F],human_MU_b3);
[p_hh,h_hh]             = ranksum(human_MU_b1, human_MU_b3);
stats.fig1D.MU_m        = p_m;     % monkey difference
stats.fig1D.MU_hm_b1    = p_hm_b1; % monkey human difference
stats.fig1D.MU_hm_b3    = p_hm_b3; % monkey human difference
stats.fig1D.MU_hz_b3    = p_hz_b3; % monkey human difference
stats.fig1D.MU_hf_b3    = p_hf_b3; % monkey human difference
stats.fig1D.MU_h_b1_b3  = p_hh;    % human-human difference

stats.fig1D.MU_median_mf    = median(MU_F);
stats.fig1D.MU_median_mz    = median(MU_Z);
stats.fig1D.MU_median_m     = median([MU_F MU_Z]);

stats.fig1D.MU_median_b1    = median(human_MU_b1);
stats.fig1D.MU_median_b3    = median(human_MU_b3);

%% Text stats:
% Percent correct extreme stimuli:

for iRec = 1:length(recInd)
    % Stimulus strength analysis:
    HCtrials                = monkey.masterBehaviorStruct{recInd(iRec)}.trials{1}{1}{1};
    uOri                    = unique(HCtrials(:,2));
    nExtremeOri = sum(HCtrials(:,2)==uOri(1) | HCtrials(:,2)==uOri(11));
    nCorrect =  sum(sign(HCtrials(HCtrials(:,2)==uOri(1),3))==-1) + sum( sign(HCtrials(HCtrials(:,2)==uOri(11),3))==1) ;
    pCorrectExtreme(iRec) = nCorrect./nExtremeOri;
    
    extremeOri(iRec) = uOri(1);
    
    if strcmp(monkey.masterBehaviorStruct{recInd(iRec)}.type,'Internal')
        exptInd(iRec) = 1;
    else
        exptInd(iRec) = 0;
    end
    
    recObsInd{iRec} =  monkey.masterBehaviorStruct{recInd(iRec)}.observer;
end
stats.fig1.pCorrectExtreme_F =median(pCorrectExtreme(strcmp('friedrich', recObsInd)));
stats.fig1.pCorrectExtreme_Z =median(pCorrectExtreme(strcmp('ziggy', recObsInd)));
stats.fig1.extremeOri_F =   median(extremeOri(strcmp('friedrich', recObsInd)));
stats.fig1.extremeOri_Z =   median(extremeOri(strcmp('ziggy', recObsInd)));

%% Stats related to behavior/methods:
% Get proportion of trials in score 4 or low scores:
for iDate = 1:length(monkey.masterBehaviorStruct)
    nScore4(iDate)              = sum(cellfun(@length, monkey.masterBehaviorStruct{iDate}.trials{1}{3}));
    nScoreLow(iDate)            = sum(cellfun(@length, monkey.masterBehaviorStruct{iDate}.trials{1}{2}));
    nTrial(iDate)               = sum(cellfun(@length, monkey.masterBehaviorStruct{iDate}.trials{1}{1}));
    percentScore4(iDate)        = nScore4(iDate)./nTrial(iDate);
    highConf_ind                = monkey.masterBehaviorStruct{iDate}.sessionStruct.trialMatrix(:,monkey.masterBehaviorStruct{iDate}.sessionStruct.trialMatrix_index.CONFIDENCE)==1;
    correct_ind                 = monkey.masterBehaviorStruct{iDate}.sessionStruct.trialMatrix(:,monkey.masterBehaviorStruct{iDate}.sessionStruct.trialMatrix_index.FEEDBACK)==1; % random assignment of zero at time of expt
    score4_ind                  = monkey.masterBehaviorStruct{iDate}.sessionStruct.trialMatrix(:,monkey.masterBehaviorStruct{iDate}.sessionStruct.trialMatrix_index.SCORE)==4;
    
    percentHighConf(iDate)      = mean( highConf_ind); 
    percentCorrect(iDate)       = mean(correct_ind); 
    
    percentHighConf_score4(iDate) = mean(highConf_ind(score4_ind));
    percentHighConf_score13(iDate) = mean(highConf_ind(~score4_ind));
end

stats.fig1.nRecF                 = sum(recordingInd & obsInd);
stats.fig1.nRecZ                 = sum(recordingInd & ~obsInd);
stats.fig1.avgPscore4            = median(percentScore4(recordingInd ));
stats.fig1.avgPscore4_F          = median(percentScore4(recordingInd & obsInd));
stats.fig1.avgPscore4_Z          = median(percentScore4(recordingInd & ~obsInd));
stats.fig1.avg_trialScore4_F     = median(nScore4(recordingInd & obsInd));
stats.fig1.avg_trialScore4_Z     = median(nScore4(recordingInd & ~obsInd));

stats.fig1.avg_p_HC_F            = median(percentHighConf(recordingInd&obsInd));
stats.fig1.avg_p_HC_Z            = median(percentHighConf(recordingInd&~obsInd));
