%% Supplemental Figures
% Written 05-29-2024 ZBS (zbsinger@mit.edu)

% Supplement 
% S1.A 
% S1.B 
% S1.C 
% S1.D 
% S1.E 

% S3.A 
% S3.B 
%% Start with clean slate
clearvars -except masterNeuralStruct
clc;
close all;

% Set Paths
user                    = '/Users/zbsinger/MIT Dropbox/Zoe Boundy-Singer/Projects/V1_Confidence_Physiology';
thisPath                = user;
structurePath           = strcat([thisPath '/Structures/']);
addpath(thisPath,structurePath)

%% Human data:
human                   = load( strcat([structurePath '/Behavior/Human/masterBehaviorStruct_Human.mat']),'masterBehaviorStruct','crossBlock', 'observerList', 'taskList' );
uncPsychList            = strcmp('UncertaintyPsych',human.taskList);
%% Monkey data:
monkey                  = load([structurePath '/Behavior/Monkey/masterBehaviorStruct.mat']); % load all monkey behavioral sessions
nDate                   = numel(monkey.masterBehaviorStruct);

recordingInd            = monkey.crossBlock.recordingInd;
obsInd                  = monkey.crossBlock.obsInd; % 1 = Friedrich; 0 = Ziggy;
recInd                  = find(recordingInd);
exDate                  = recInd(19); % 02-08-2023 (this is the 20th entry of recInd and 83rd overall entry)
nDate                   = sum(recordingInd);


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
observer{1}             = 'ziggy';
observer{2}             = 'friedrich';
observerInd{1}          = find(strcmp('ziggy',observerList)); % ZIGGY
observerInd{2}          = find(strcmp('friedrich',observerList)); % FRIEDRICH
obsIndRec               = strcmp('friedrich',observerList);

%% Human values:
iT=1; % iT 1 pulls out full model fit on all data
for iO = 1:length(human.observerList)
    
    % Stimulus strength analysis:
    HCtrials                = human.masterBehaviorStruct{iO}.trials{1}{iT}{1};
    uOri                    = unique(HCtrials(:,2));
    
    if length(uOri)==11
        largeOriInd             = HCtrials(:,2) == uOri(1) | HCtrials(:,2) == uOri(2) | HCtrials(:,2) == uOri(10) | HCtrials(:,2) == uOri(11);
        smallOriInd             = HCtrials(:,2) == uOri(4) | HCtrials(:,2) == uOri(5) | HCtrials(:,2) == uOri(7) | HCtrials(:,2) == uOri(8);
        
        percentHC_largeOri(iO)  = mean(abs(HCtrials(largeOriInd,3))==2); % Column 3 is signed confdeince; 2 = HC
        percentHC_smallOri(iO)  = mean(abs(HCtrials(smallOriInd,3))==2);
        
        % Low-contrast
        LCtrials                 =  human.masterBehaviorStruct{iO}.trials{1}{iT}{2};
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
    
    human.uncRatio_highCont(iO) =log10( 1./(split_LC(iO,1)./split_HC(iO,1))) ; % ratio
    human.uncRatio_lowCont(iO) = log10(1./(split_LC(iO,2)./split_HC(iO,2))) ;


    human.HC_sens(iO) = human.masterBehaviorStruct{iO}.paramEst.values{1}{iT}(3);
    human.LC_sens(iO) =  human.masterBehaviorStruct{iO}.paramEst.values{1}{iT}(4);
    
    human.uncRatio_contrast(iO) = human.LC_sens(iO)./human.HC_sens(iO);
    %% GOF:
    %[GOFval]                                    = behaviorGOF( human.masterBehaviorStruct{iO}.trials{1}{1}, human.masterBehaviorStruct{iO}.paramEst.values{1}{1});
    %GOF(iO,:)
    
%     %% Meta-uncertainty
%     if size(human.masterBehaviorStruct{iO}.paramEst.values{1},2)>3 % If subject didn't have 1st block thrown out
%         human_meta_unc_block1(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{4}(2);% block 1
%         human_meta_unc_block2(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{5}(2); % block 2
%         human_meta_unc_block3(iO)       = human.masterBehaviorStruct{iO}.paramEst.values{1}{6}(2); % block 3
%      else
%         human_meta_unc_block1(iO)       = nan;
%         human_meta_unc_block2(iO)       = nan;
%         human_meta_unc_block3(iO)       = nan;
%     end
end

%% S1.A  Example behavior for supplements (F230208)

for exampleDate = [19]
    iT=1;
    sessionStruct           = monkey.masterBehaviorStruct{recInd(exampleDate)}.sessionStruct;
    paramEst                = monkey.masterBehaviorStruct{recInd(exampleDate)}.paramEst;
    trials                  = monkey.masterBehaviorStruct{recInd(exampleDate)}.trials;
    nRel                    = size(trials{1}{1},2);
    asymFlag                = 1;
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
        nParams       = 2 + nRel + nRel + nConfCrit+1; % [Guess rate, meta-noise], [stimulus sensitivity], [stimulus criterion], [confidence criteria]
    else
        nParams       = 2 + nRel + nRel + nConfCrit;
    end
    iO              = 1; % one observer
    
    stats.Sfig1A.nTrials   = length(trials{iO}{iT}{1})+length(trials{iO}{iT}{2}); % number of total trials for example session
    for iR=1:nRel
        % Set experiment parameters
        stimCat{iR}   = trials{iO}{iT}{iR}(:,1);
        stimValue{iR} = trials{iO}{iT}{iR}(:,2);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
        choice        = trials{iO}{iT}{iR}(:,3);
        
        % Recode responses
        respOptions   = [-(nConfCrit+1):-1, 1:(nConfCrit+1)];
        choiceVec{iR} = zeros(2*(nConfCrit+1), numel(stimValue{iR}));
        
        for iC = 1:numel(respOptions)
            selTrials                       = (choice == respOptions(iC));
            choiceVec{iR}(iC, selTrials)    = 1;
        end
        
        oriPF{iR}{iT}       = unique(stimValue{iR});
        nPoints{iR}{iT}         = numel(oriPF{iR}{iT});
        
        confidence                                                      = zeros(length(choice),1);
        confidence(sum(choiceVec{iR}([1 4],:))==1)                      = 1; % 1 = HC
        
        %%
        for iP = 1:nPoints{iR}{iT}
            oriInd                         = stimValue{iR} == oriPF{iR}{iT}(iP);
            nTrials{iR}{iT}(iP)             = sum(oriInd);
            obsPF{iR}{iT}(iP)               = mean(choice(oriInd) > 0);
            obsCF{iR}{iT}(iP)               = mean(abs(choice(oriInd))>1);
            
            % Psychometric for low and high confidence trials
            indHighConf                 = abs(choice) >= 2;
            obsPF_HC{iR}{iT}(iP)            = mean(choice(oriInd & indHighConf) > 0);
            obsPF_LC{iR}{iT}(iP)            = mean(choice(oriInd & ~indHighConf) > 0);
            nTrials_HC{iR}{iT}(iP)          = sum(oriInd & indHighConf);
            nTrials_LC{iR}{iT}(iP)          = sum(oriInd & ~indHighConf);
            
        end
    end
    
    %% Step 3: Plot model
    % First select appropriate subset of parameters
    % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta uncertianty, conf criteria]
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
        params = paramEst.split{1}{1}{iR};
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
        
        titleText          = [monkey.masterBehaviorStruct{recInd(exampleDate)}.date];
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
stat.Sfig1A.exRatio = monkey.crossBlock.LC_sens(recInd(exampleDate))./monkey.crossBlock.HC_sens(recInd(exampleDate));

%% Monkey effect stimulus strength on confidence
monkey.percentHC_largeOriHC    = nan(length(recordingInd),1);
monkey.percentHC_smallOriHC    = nan(length(recordingInd),1);
monkey.percentHC_largeOriLC    = nan(length(recordingInd),1);
monkey.percentHC_smallOriLC    = nan(length(recordingInd),1);

for iO = 1:length(recordingInd)
    % High contrast; stim strength analysis
    HCtrials                = monkey.masterBehaviorStruct{iO}.trials{1}{1}{1};
    uOri                    = unique(HCtrials(:,2));
    
    largeOriInd             = HCtrials(:,2) == uOri(1) | HCtrials(:,2) == uOri(2) | HCtrials(:,2) == uOri(10) | HCtrials(:,2) == uOri(11);
    smallOriInd             = HCtrials(:,2) == uOri(4) | HCtrials(:,2) == uOri(5) | HCtrials(:,2) == uOri(7) | HCtrials(:,2) == uOri(8);
    monkey.percentHC_largeOriHC(iO)  = mean(abs(HCtrials(largeOriInd,3))==2); % Column 3 is signed confdeince; 2 = HC
    monkey.percentHC_smallOriHC(iO)  = mean(abs(HCtrials(smallOriInd,3))==2);
    
    % Low-contrast
    LCtrials                 = monkey.masterBehaviorStruct{iO}.trials{1}{1}{2};
    uOriLC                   = unique(LCtrials(:,2));
    nOri(iO)                 = numel(uOriLC);
    
    if nOri(iO)==11
        smallOriInd_HC      = HCtrials(:,2) == uOri(4) | HCtrials(:,2) == uOri(5) | HCtrials(:,2) == uOri(7) | HCtrials(:,2) == uOri(8);
        smallOriInd_LC      = LCtrials(:,2) == uOriLC(4) | LCtrials(:,2) == uOriLC(5) | LCtrials(:,2) == uOriLC(7) | LCtrials(:,2) == uOriLC(8);
        largeOriInd_LC      = LCtrials(:,2) == uOriLC(1) | LCtrials(:,2) == uOriLC(2) | LCtrials(:,2) == uOriLC(10) | LCtrials(:,2) == uOriLC(11);
        
        percentHC_HC(iO)    = mean(abs(HCtrials(smallOriInd_HC,3))==2);
        avgHCori(iO)        = mean(abs(HCtrials( smallOriInd_HC ,2)));
        
        percentHC_LC(iO)    = mean(abs(LCtrials(smallOriInd_LC,3))==2);
        avgLCori(iO)        =  mean(abs(LCtrials( smallOriInd_LC ,2)));
        
        % stim strength analysis for LC
        monkey.percentHC_largeOriLC(iO)  = mean(abs(LCtrials(largeOriInd_LC,3))==2); % Column 3 is signed confidence; 2 = HC
        monkey.percentHC_smallOriLC(iO)  = mean(abs(LCtrials(smallOriInd_LC,3))==2);
    elseif nOri(iO)==9
        
        smallOriInd_HC      = HCtrials(:,2) == uOri(4) | HCtrials(:,2) == uOri(5) | HCtrials(:,2) == uOri(7) | HCtrials(:,2) == uOri(8);
        smallOriInd_LC      = LCtrials(:,2) == uOriLC(4) | LCtrials(:,2) == uOriLC(3) | LCtrials(:,2) == uOriLC(6) | LCtrials(:,2) == uOriLC(7);
        largeOriInd_LC      = LCtrials(:,2) == uOriLC(1) | LCtrials(:,2) == uOriLC(2) | LCtrials(:,2) == uOriLC(8) | LCtrials(:,2) == uOriLC(9);
        
        percentHC_HC(iO)    = mean(abs(HCtrials(smallOriInd_HC,3))==2);
        avgHCori(iO)        = mean(abs(HCtrials( smallOriInd_HC ,2)));
        
        percentHC_LC(iO)    = mean(abs(LCtrials(smallOriInd_LC,3))==2);
        avgLCori(iO)        =  mean(abs(LCtrials( smallOriInd_LC ,2)));
        
        % stim strength analysis for LC
        monkey.percentHC_largeOriLC(iO)  = mean(abs(LCtrials(largeOriInd_LC,3))==2); % Column 3 is signed confdeince; 2 = HC
        monkey.percentHC_smallOriLC(iO)  = mean(abs(LCtrials(smallOriInd_LC,3))==2);
    end
    
    %% Analysis of proporition CW choices given CW stimuli vs CCW stimuli
    CW_ori_ind      = monkey.masterBehaviorStruct{iO}.sessionStruct.trialMatrix(:,monkey.masterBehaviorStruct{iO}.sessionStruct.trialMatrix_index.ORIENTATION)>0;
    CCW_ori_ind     = monkey.masterBehaviorStruct{iO}.sessionStruct.trialMatrix(:,monkey.masterBehaviorStruct{iO}.sessionStruct.trialMatrix_index.ORIENTATION)<0;
    choice_colum    = monkey.masterBehaviorStruct{iO}.sessionStruct.trialMatrix_index.SIDE;
    pCW_CW          = mean(monkey.masterBehaviorStruct{iO}.sessionStruct.trialMatrix(CW_ori_ind,choice_colum)==1);
    pCCW_CW         = mean(monkey.masterBehaviorStruct{iO}.sessionStruct.trialMatrix(CCW_ori_ind,choice_colum)==1);
    pCW_diff(iO)    = pCW_CW -pCCW_CW;
end
m_Z     = mean(pCW_diff(recordingInd& ~obsInd));
m_F     = mean(pCW_diff(recordingInd& obsInd));
m_rec   = mean(pCW_diff(recordingInd));

%% Human stimulus strength effect on confidence
human.percentHC_largeOriHC    = nan(sum(uncPsychList),1);
human.percentHC_smallOriHC    = nan(sum(uncPsychList),1);
human.percentHC_largeOriLC    = nan(sum(uncPsychList),1);
human.percentHC_smallOriLC    = nan(sum(uncPsychList),1);

for iO = 1:11 % only the uncPsych (task) subjects
    % High contrast; stim strength analysis
    HCtrials                = human.masterBehaviorStruct{iO}.trials{1}{1}{1};
    uOri                    = unique(HCtrials(:,2));
    
    largeOriInd             = HCtrials(:,2) == uOri(1) | HCtrials(:,2) == uOri(2) | HCtrials(:,2) == uOri(10) | HCtrials(:,2) == uOri(11);
    smallOriInd             = HCtrials(:,2) == uOri(4) | HCtrials(:,2) == uOri(5) | HCtrials(:,2) == uOri(7) | HCtrials(:,2) == uOri(8);
    human.percentHC_largeOriHC(iO)  = mean(abs(HCtrials(largeOriInd,3))==2); % Column 3 is signed confdeince; 2 = HC
    human.percentHC_smallOriHC(iO)  = mean(abs(HCtrials(smallOriInd,3))==2);
    
    % Low-contrast
    LCtrials                 = human.masterBehaviorStruct{iO}.trials{1}{1}{2};
    uOriLC                   = unique(LCtrials(:,2));
    
    largeOriIndLC             = LCtrials(:,2) == uOriLC(1) | LCtrials(:,2) == uOriLC(2) | LCtrials(:,2) == uOriLC(10) | LCtrials(:,2) == uOriLC(11);
    smallOriIndLC             = LCtrials(:,2) == uOriLC(4) | LCtrials(:,2) == uOriLC(5) | LCtrials(:,2) == uOriLC(7) | LCtrials(:,2) == uOriLC(8);
    human.percentHC_largeOriLC(iO)  = mean(abs(LCtrials(largeOriIndLC,3))==2); % Column 3 is signed confdeince; 2 = HC
    human.percentHC_smallOriLC(iO)  = mean(abs(LCtrials(smallOriIndLC,3))==2);
    
end

%% S1.B High contrast perceptual uncertainty (monkey F vs Z vs Human) (supplement)
% Sensitivity comes from CASANDRE fit
HC_sens_F       = 1./monkey.crossBlock.HC_sens(obsInd&recordingInd);
HC_sens_Z       = 1./monkey.crossBlock.HC_sens(~obsInd&recordingInd);
HC_sens_H       = 1./human.HC_sens(uncPsychList);
HC_sens         = [HC_sens_F,HC_sens_Z,HC_sens_H]';

% Calculate percentiles:
h_HC_sens_p = prctile(HC_sens_H,[25 50 75]);
z_HC_sens_p = prctile(HC_sens_Z,[25 50 75]);
f_HC_sens_p = prctile(HC_sens_F,[25 50 75]);

figure;hold on; axis square
xlim([0 4])
ylim([0 5])
plot(ones(length(HC_sens_H),1)+(rand(length(HC_sens_H),1)-.5)/5,HC_sens_H,'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length(HC_sens_Z),1)+(rand(length(HC_sens_Z),1)-.5)/5,HC_sens_Z,'o','markerfacecolor','r','markeredgecolor','k','markersize',10) % Ziggy
plot(2+ones(length(HC_sens_F),1)+(rand(length(HC_sens_F),1)-.5)/5,HC_sens_F,'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Friedrich

plot([1 1],h_HC_sens_p([1 3]))
plot([1 1]+1,z_HC_sens_p([1 3]))
plot([1 1]+2,f_HC_sens_p([1 3]))

plot([1 ],h_HC_sens_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,z_HC_sens_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,f_HC_sens_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Observer','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Slope (high contrast)','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({'Human',  'Ziggy', 'Friedrich'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Uncertainty param (high contrast)','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)


% Stats
[p_m,h_m]           = ranksum(HC_sens_F,HC_sens_Z);
[p_hm,h_hm]         = ranksum([HC_sens_F,HC_sens_Z], HC_sens_H);
[p_hz]              = ranksum([HC_sens_Z], HC_sens_H);
[p_hf]              = ranksum([HC_sens_F], HC_sens_H);


stats.Sfig1B.C_sens_m        = p_m; % monkey difference
stats.Sfig1B.HC_sens_hm      = p_hm; % monkey (grouped) human difference
stats.Sfig1B.HC_sens_hz      = p_hz; % monkey Z human difference
stats.Sfig1B.HC_sens_hf      = p_hf; % monkey F human difference


%% S1.C low:high contrast ratio
uncRatio_cont_F               = monkey.crossBlock.LC_sens(obsInd&recordingInd)./monkey.crossBlock.HC_sens(obsInd&recordingInd);
uncRatio_cont_Z               = monkey.crossBlock.LC_sens(~obsInd&recordingInd)./monkey.crossBlock.HC_sens(~obsInd&recordingInd);
uncRatio_cont_H               = human.uncRatio_contrast;

% Calculate percentiles:
h_uncRatio_cont_p = prctile(uncRatio_cont_H,[25 50 75]);
z_uncRatio_cont_p = prctile(uncRatio_cont_Z,[25 50 75]);
f_uncRatio_cont_p = prctile(uncRatio_cont_F,[25 50 75]);

figure;hold on; axis square
xlim([0 4])
ylim([0 1.5])
plot(ones(length(uncRatio_cont_H),1)+(rand(length(uncRatio_cont_H),1)-.5)/5,uncRatio_cont_H,'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length(uncRatio_cont_Z),1)+(rand(length(uncRatio_cont_Z),1)-.5)/5,uncRatio_cont_Z,'o','markerfacecolor','r','markeredgecolor','k','markersize',10) % Ziggy
plot(2+ones(length(uncRatio_cont_F),1)+(rand(length(uncRatio_cont_F),1)-.5)/5,uncRatio_cont_F,'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Friedrich

plot([1 1],h_uncRatio_cont_p([1 3]))
plot([1 1]+1,z_uncRatio_cont_p([1 3]))
plot([1 1]+2,f_uncRatio_cont_p([1 3]))

plot([1 ],h_uncRatio_cont_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,z_uncRatio_cont_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,f_uncRatio_cont_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Observer','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Ratio Low:high contrast','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({'Human',  'Ziggy', 'Friedrich'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Slope ratio','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

[p_hz]              = ranksum([uncRatio_cont_Z], uncRatio_cont_H);
[p_hf]              = ranksum([uncRatio_cont_F], uncRatio_cont_H);
stats.Sfig1C.ratio_sens_hz      = p_hz; % monkey Z human difference
stats.Sfig1C.ratio_sens_hf      = p_hf; % monkey F human difference

%% S1.F Effect of stimulus strength
stimStrength_F               = [(monkey.percentHC_largeOriHC(obsInd&recordingInd) - monkey.percentHC_smallOriHC(obsInd&recordingInd)); (monkey.percentHC_largeOriLC(obsInd&recordingInd) -monkey.percentHC_smallOriLC(obsInd&recordingInd))];
stimStrength_Z               = [(monkey.percentHC_largeOriHC(~obsInd&recordingInd) - monkey.percentHC_smallOriHC(~obsInd&recordingInd)); (monkey.percentHC_largeOriLC(~obsInd&recordingInd) -monkey.percentHC_smallOriLC(~obsInd&recordingInd))];
stimStrength_H               = [(human.percentHC_largeOriHC-human.percentHC_smallOriHC); (human.percentHC_largeOriLC -human.percentHC_smallOriLC)];

% Calculate percentiles:
h_stimStrength_p = prctile(stimStrength_H,[25 50 75]);
z_stimStrength_p = prctile(stimStrength_Z,[25 50 75]);
f_stimStrength_p = prctile(stimStrength_F,[25 50 75]);

figure;hold on; axis square
xlim([0 4])
ylim([-.5 1])
plot(ones(length(stimStrength_H),1)+(rand(length(stimStrength_H),1)-.5)/5,stimStrength_H,'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length(stimStrength_Z),1)+(rand(length(stimStrength_Z),1)-.5)/5,stimStrength_Z,'o','markerfacecolor','r','markeredgecolor','k','markersize',10) % Ziggy
plot(2+ones(length(stimStrength_F),1)+(rand(length(stimStrength_F),1)-.5)/5,stimStrength_F,'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Friedrich

plot([1 1],h_stimStrength_p([1 3]))
plot([1 1]+1,z_stimStrength_p([1 3]))
plot([1 1]+2,f_stimStrength_p([1 3]))

plot([1 ],h_stimStrength_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,z_stimStrength_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,f_stimStrength_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Observer','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({'Human',  'Ziggy', 'Friedrich'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Diff large- small stim strength % HC','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

%% version 2.0
obsColor        = {'y','b','r'};
stimStrengthLarge_F               = [monkey.percentHC_largeOriHC(obsInd&recordingInd); monkey.percentHC_largeOriLC(obsInd&recordingInd)];
stimStrengthSmall_F               = [monkey.percentHC_smallOriHC(obsInd&recordingInd); monkey.percentHC_smallOriLC(obsInd&recordingInd)];

stimStrengthLarge_Z               = [monkey.percentHC_largeOriHC(~obsInd&recordingInd); monkey.percentHC_largeOriLC(~obsInd&recordingInd)];
stimStrengthSmall_Z               = [monkey.percentHC_smallOriHC(~obsInd&recordingInd); monkey.percentHC_smallOriLC(~obsInd&recordingInd)];

stimStrengthLarge_H               = [human.percentHC_largeOriHC; human.percentHC_largeOriLC];
stimStrengthSmall_H               = [human.percentHC_smallOriHC; human.percentHC_smallOriLC];

% Calculate percentiles:
h_stimStrengthLarge_p = prctile(stimStrengthLarge_H,[25 50 75]);
h_stimStrengthSmall_p = prctile(stimStrengthSmall_H,[25 50 75]);

f_stimStrengthLarge_p = prctile(stimStrengthLarge_F,[25 50 75]);
f_stimStrengthSmall_p = prctile(stimStrengthSmall_F,[25 50 75]);

z_stimStrengthLarge_p = prctile(stimStrengthLarge_Z,[25 50 75]);
z_stimStrengthSmall_p = prctile(stimStrengthSmall_Z,[25 50 75]);


figure;hold on; axis square
xlim([0 7])
ylim([-0.1 1])
plot([0 6],[0 0])

for iO=1:length(stimStrengthSmall_H  )
    % human
randVal             = (rand(2,1)-.5)/5;
plot(1+randVal(1),  [stimStrengthLarge_H(iO)],'o','markerfacecolor',obsColor{3},'markeredgecolor','k','markersize',10)
plot(2+randVal(2),  [stimStrengthSmall_H(iO)],'o','markerfacecolor',obsColor{3},'markeredgecolor','k','markersize',10)
plot([1+randVal(1) 2+randVal(2)],  [stimStrengthLarge_H(iO) stimStrengthSmall_H(iO)],'-','color',obsColor{3})
end
for iO=1:length(stimStrengthSmall_F  )
% F
randVal             = (rand(2,1)-.5)/5;
plot(3+randVal(1),  [stimStrengthLarge_F(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot(4+randVal(2),  [stimStrengthSmall_F(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot([3+randVal(1) 4+randVal(2)],  [stimStrengthLarge_F(iO) stimStrengthSmall_F(iO)],'-','color',obsColor{obsInd(iO)+1})
end
for iO=1:length(stimStrengthSmall_Z  )
%Z
randVal             = (rand(2,1)-.5)/5;
plot(5+randVal(1),  [stimStrengthLarge_Z(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+2},'markeredgecolor','k','markersize',10)
plot(6+randVal(2),  [stimStrengthSmall_Z(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+2},'markeredgecolor','k','markersize',10)
plot([5+randVal(1) 6+randVal(2)],  [stimStrengthLarge_Z(iO) stimStrengthSmall_Z(iO)],'-','color',obsColor{obsInd(iO)+2})
end

plot([1 1],h_stimStrengthLarge_p([1 3]))
plot([1 1]+1,h_stimStrengthSmall_p([1 3]))

plot([1 1]+2,f_stimStrengthLarge_p([1 3]))
plot([1 1]+3,f_stimStrengthSmall_p([1 3]))

plot([1 1]+4,z_stimStrengthLarge_p([1 3]))
plot([1 1]+5,z_stimStrengthSmall_p([1 3]))

ylabel('Stim strength','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3 4 5 6 ])
xticklabels({'Human large','H small','F large' ,'F small','Z large','Z small'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

%stats 1E
[p_f]              = signrank(stimStrengthLarge_F ,stimStrengthSmall_F);
[p_z]              = signrank(stimStrengthLarge_Z ,stimStrengthSmall_Z);
[p_h]              = signrank(stimStrengthLarge_H ,stimStrengthSmall_H);

stats.Sfig1E.ss_f     = p_f;
stats.Sfig1E.ss_z     = p_z; 
stats.Sfig1E.ss_h     = p_h;

%% Effect of stimulus reliability on confidence
stimRel_F               = [monkey.crossBlock.predCF(recordingInd & obsInd,1) - monkey.crossBlock.predCF(recordingInd & obsInd,2)];
stimRel_Z               = [monkey.crossBlock.predCF(recordingInd & ~obsInd,1) - monkey.crossBlock.predCF(recordingInd & ~obsInd,2)];
stimRel_H               = [human.crossBlock.predCF(uncPsychList,1) - human.crossBlock.predCF(uncPsychList,2)];

% Calculate percentiles:
h_stimRel_p = prctile(stimRel_H,[25 50 75]);
z_stimRel_p = prctile(stimRel_Z,[25 50 75]);
f_stimRel_p = prctile(stimRel_F,[25 50 75]);

figure;hold on; axis square
xlim([0 4])
ylim([-.5 1])
plot(ones(length(stimRel_H),1)+(rand(length(stimRel_H),1)-.5)/5,stimRel_H,'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length(stimRel_Z),1)+(rand(length(stimRel_Z),1)-.5)/5,stimRel_Z,'o','markerfacecolor','r','markeredgecolor','k','markersize',10) % Ziggy
plot(2+ones(length(stimRel_F),1)+(rand(length(stimRel_F),1)-.5)/5,stimRel_F,'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Friedrich

plot([1 1],h_stimRel_p([1 3]))
plot([1 1]+1,z_stimRel_p([1 3]))
plot([1 1]+2,f_stimRel_p([1 3]))

plot([1 ],h_stimRel_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,z_stimRel_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,f_stimRel_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Observer','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({'Human',  'Ziggy', 'Friedrich'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Diff pred HC high - low contrast(%) ','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)


%% Version 2.0
stimRel_HC_F               = [monkey.crossBlock.predCF(recordingInd & obsInd,1) ];
stimRel_LC_F               = [monkey.crossBlock.predCF(recordingInd & obsInd,2)];

stimRel_HC_Z               = [monkey.crossBlock.predCF(recordingInd & ~obsInd,1) ];
stimRel_LC_Z               = [monkey.crossBlock.predCF(recordingInd & ~obsInd,2)];

stimRel_HC_H               = [human.crossBlock.predCF(uncPsychList,1)];
stimRel_LC_H               = [human.crossBlock.predCF(uncPsychList,2)];

% Calculate percentiles:
h_stimRel_HC_p = prctile(stimRel_HC_H,[25 50 75]);
h_stimRel_LC_p = prctile(stimRel_LC_H,[25 50 75]);

f_stimRel_HC_p = prctile(stimRel_HC_F,[25 50 75]);
f_stimRel_LC_p = prctile(stimRel_LC_F,[25 50 75]);

z_stimRel_HC_p = prctile(stimRel_HC_Z,[25 50 75]);
z_stimRel_LC_p = prctile(stimRel_LC_Z,[25 50 75]);

figure;hold on; axis square
xlim([0 7])
ylim([-0.1 1])
plot([0 6],[0 0])

for iO=1:length(stimRel_HC_H  )
    % human
randVal             = (rand(2,1)-.5)/5;
plot(1+randVal(1),  [stimRel_HC_H(iO)],'o','markerfacecolor',obsColor{3},'markeredgecolor','k','markersize',10)
plot(2+randVal(2),  [stimRel_LC_H(iO)],'o','markerfacecolor',obsColor{3},'markeredgecolor','k','markersize',10)
plot([1+randVal(1) 2+randVal(2)],  [stimRel_HC_H(iO) stimRel_LC_H(iO)],'-','color',obsColor{3})
end
for iO=1:length(stimRel_LC_F  )
% F
randVal             = (rand(2,1)-.5)/5;
plot(3+randVal(1),  [stimRel_HC_F(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot(4+randVal(2),  [stimRel_LC_F(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot([3+randVal(1) 4+randVal(2)],  [stimRel_HC_F(iO) stimRel_LC_F(iO)],'-','color',obsColor{obsInd(iO)+1})
end
for iO=1:length(stimRel_LC_Z  )
%Z
randVal             = (rand(2,1)-.5)/5;
plot(5+randVal(1),  [stimRel_HC_Z(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+2},'markeredgecolor','k','markersize',10)
plot(6+randVal(2),  [stimRel_LC_Z(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+2},'markeredgecolor','k','markersize',10)
plot([5+randVal(1) 6+randVal(2)],  [stimRel_HC_Z(iO) stimRel_LC_Z(iO)],'-','color',obsColor{obsInd(iO)+2})
end

plot([1 1],h_stimRel_HC_p([1 3]))
plot([1 1]+1,h_stimRel_LC_p([1 3]))

plot([1 1]+2,f_stimRel_HC_p([1 3]))
plot([1 1]+3,f_stimRel_LC_p([1 3]))

plot([1 1]+4,z_stimRel_HC_p([1 3]))
plot([1 1]+5,z_stimRel_LC_p([1 3]))

ylabel('Rel test','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3 4 5 6 ])
xticklabels({'Human large','H small','F large' ,'F small','Z large','Z small'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

%stats 1E
[p_f]              = signrank(stimRel_HC_F ,stimRel_LC_F);
[p_z]              = signrank(stimRel_HC_Z ,stimRel_LC_Z);
[p_h]              = signrank(stimRel_HC_H ,stimRel_LC_H);

stats.Sfig1E.cont_f     = p_f;
stats.Sfig1E.cont_z      = p_z; 
stats.Sfig1E.cont_h     = p_h;
%% S3.A Vector angle 
NL  = load(strcat([structurePath  '/Neural//decodeSummary_NL.mat']));
L   = load(strcat([structurePath  '/Neural/decodeSummary_L.mat']));

choiceStim_p    = prctile(L.crossDate.avgChoiceStimAngle,[25 50 75]);
choiceConf_p    = prctile(L.crossDate.avgChoiceConfAngle,[25 50 75]);
confStim_p      = prctile(L.crossDate.avgConfStimAngle,[25 50 75]);
AmbchoiceStim_p = prctile(L.crossDate.avgAmbChoiceStimAngle,[25 50 75]);


figure;hold on; axis square
xlim([0 5])
ylim([0 180])
plot([0 5],[90 90])

plot(ones(length( L.crossDate.avgChoiceStimAngle(obsIndRec)),1)+(rand(length( L.crossDate.avgChoiceStimAngle(obsIndRec)),1)-.5)/5, L.crossDate.avgChoiceStimAngle(obsIndRec),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(ones(length( L.crossDate.avgChoiceStimAngle(~obsIndRec)),1)+(rand(length( L.crossDate.avgChoiceStimAngle(~obsIndRec)),1)-.5)/5, L.crossDate.avgChoiceStimAngle(~obsIndRec),'o','markerfacecolor','b','markeredgecolor','k','markersize',10)

plot(1+ones(length(L.crossDate.avgChoiceConfAngle(obsIndRec)),1)+(rand(length(L.crossDate.avgChoiceConfAngle(obsIndRec)),1)-.5)/5,L.crossDate.avgChoiceConfAngle(obsIndRec),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)  % F
plot(1+ones(length(L.crossDate.avgChoiceConfAngle(~obsIndRec)),1)+(rand(length(L.crossDate.avgChoiceConfAngle(~obsIndRec)),1)-.5)/5,L.crossDate.avgChoiceConfAngle(~obsIndRec),'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Z

plot(2+ones(length(L.crossDate.avgConfStimAngle(obsIndRec)),1)+(rand(length(L.crossDate.avgConfStimAngle(obsIndRec)),1)-.5)/5,L.crossDate.avgConfStimAngle(obsIndRec),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)  % F
plot(2+ones(length(L.crossDate.avgConfStimAngle(~obsIndRec)),1)+(rand(length(L.crossDate.avgConfStimAngle(~obsIndRec)),1)-.5)/5,L.crossDate.avgConfStimAngle(~obsIndRec),'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Z

plot(3+ones(length(L.crossDate.avgAmbChoiceStimAngle(obsIndRec)),1)+(rand(length(L.crossDate.avgAmbChoiceStimAngle(obsIndRec)),1)-.5)/5,L.crossDate.avgAmbChoiceStimAngle(obsIndRec),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)  % F
plot(3+ones(length(L.crossDate.avgAmbChoiceStimAngle(~obsIndRec)),1)+(rand(length(L.crossDate.avgAmbChoiceStimAngle(~obsIndRec)),1)-.5)/5,L.crossDate.avgAmbChoiceStimAngle(~obsIndRec),'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Z


plot([1 1],choiceStim_p([1 3]))
plot([1 1]+1,choiceConf_p([1 3]))
plot([1 1]+2,confStim_p([1 3]))
plot([1 1]+3,AmbchoiceStim_p([1 3]))

plot([1 ],choiceStim_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,choiceConf_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,confStim_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+3,AmbchoiceStim_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Comparison','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3 4])
xticklabels({'Choice-stim',  'Choice-conf','Conf-stim','Amb-stim'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Avg vec angle','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

% Stats S2A
[p_choice_stim,h]                   = signrank(90-L.crossDate.avgChoiceStimAngle);
[p_choice_conf,h]                   = signrank(90-L.crossDate.avgChoiceConfAngle);
[p_conf_stim,h]                     = signrank(90-L.crossDate.avgConfStimAngle);
[p_amb_stim,h]                      = signrank(90-L.crossDate.avgAmbChoiceStimAngle);

stats.Sfig2A.p_choice_stim          = p_choice_stim;
stats.Sfig2A.p_choice_conf          = p_choice_conf;
stats.Sfig2A.p_conf_stim            = p_conf_stim;
stats.Sfig2A.p_amb_stim             = p_amb_stim;
%% 3.B N unit vs decoding performance
for iDate=1:29
  nUnit(iDate)= NL.popStruct{iDate}.nUnit ;
end

figure;hold on; axis square
%xlim([0.5 1])
ylim([0.5 1])
pl(1) = plot(nUnit(obsIndRec),L.crossDate.stimPerf(obsIndRec),'o','markerfacecolor','y','markeredgecolor','k','markersize',10);
pl(2) = plot(nUnit(~obsIndRec),L.crossDate.stimPerf(~obsIndRec),'o','markerfacecolor','b','markeredgecolor','k','markersize',10);
legend(pl,{'F','Z'},'fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xlabel('N unit','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Stim decoding performance','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)


[r_unit,p_unit]                 = corr(nUnit',L.crossDate.stimPerf');
stats.sfig2B.corr_unit_stimD    = r_unit;
stats.sfig2B.p_unit_stimD       = p_unit;

%% 3.B Choice probability vs stimulus decoding (using non-linear decoders)
% Confidence
figure;hold on; axis square
xlim([.4 .7])
ylim([0.5 1])
plot(NL.CP_conf(obsIndRec),L.crossDate.stimPerf(obsIndRec),'o','markerfacecolor','y','markeredgecolor','k','markersize',10);
plot(NL.CP_conf(~obsIndRec),L.crossDate.stimPerf(~obsIndRec),'o','markerfacecolor','b','markeredgecolor','k','markersize',10);
xlabel('Stimulus decoding performance','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Confidence CP','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
title(['Corr = ' num2str(round(corr(NL.CP_conf',L.crossDate.stimPerf'),2))],'fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

% Choice
figure;hold on; axis square
xlim([.4 .7])
ylim([0.5 1])
plot(NL.CP_choice(obsIndRec),L.crossDate.stimPerf(obsIndRec),'o','markerfacecolor','y','markeredgecolor','k','markersize',10);
plot(NL.CP_choice(~obsIndRec),L.crossDate.stimPerf(~obsIndRec),'o','markerfacecolor','b','markeredgecolor','k','markersize',10);
xlabel('Stimulus decoding performance','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Choice CP','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
title(['Corr = ' num2str(round(corr(NL.CP_choice',L.crossDate.stimPerf'),2))],'fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

% Stats
[r_conf,p_conf]                     = corr(NL.CP_conf',L.crossDate.stimPerf');
[r_choice,p_choice]                 = corr(NL.CP_choice',L.crossDate.stimPerf');

stats.Sfig2B.corr_conf_stimD        = r_conf;
stats.Sfig2B.p_conf_stimD           = p_conf;
stats.Sfig2B.corr_choice_stimD      = r_choice;
stats.Sfig2B.p_choice_stimD         = p_choice;