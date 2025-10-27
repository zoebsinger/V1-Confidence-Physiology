%% Figure 3
% Written 05-08-2024 ZBS (zbsinger@mit.edu)

%% Neural

% Figures:
% 3.A Example "neurometric" function - using conf DV to condition observed
% choices for example population
% 3.B Stat plot for 3.A
% 3.C Example population confidence decoder output (simulated activity)
% 3.D Stat plots for 3.C

%% Start with clean slate
clearvars -except masterNeuralStruct
clc;
close all;

% Set Paths
user                = '/Users/zbsinger/MIT Dropbox/Zoe Boundy-Singer/Projects/V1_Confidence_Physiology';
thisPath            = user;
structurePath       = strcat([thisPath, '/Structures/Neural']);
CSVpath             = strcat([thisPath, '/Structures/CSV']);
functionPath        = strcat([thisPath, '/Functions/']);
addpath([thisPath,structurePath,CSVpath, functionPath])

if ~exist('masterNeuralStruct')  % only load if has been deleted
    load(strcat([structurePath, '/masterNeuralStruct.mat']))
end

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

exampleDate     = 7; % 08-21-2021


%% 3.A/B Median split analysis
% If split confidence DV by median (high vs low) - is there is difference
% in the slope of the observed psychometric function (behavior)

NL                              = load(strcat([structurePath  '/decodeSummary_NL.mat']));
popStruct                       = NL.popStruct;
exampleDateMedSplit             = 8; % to plot split
exampleContrast                 = 1;
ratioMedianSplit                = nan(2,29);

for iDate=1:29
    for iC=1:2
        ind         = popStruct{iDate}.contrast==iC;% Is there a psychometric function split
        
        dvConf      = popStruct{iDate}.dvConf(ind);
        choice      = popStruct{iDate}.obsChoice(ind);
        
        ori         = popStruct{iDate}.ori(ind);
        cont        = popStruct{iDate}.contrast(ind);
        splitInd    = median(dvConf)<=dvConf; % 0 if less than average
        
        uOri        = unique(ori);
       
        if length(uOri)>=7
            for iO = 1:length(uOri) % For each unique orientation
                oriInd          = ori == uOri(iO);
                nH(iO)          = sum(splitInd & oriInd);
                nL(iO)          = sum(~splitInd & oriInd);
                
                % Number of CW observed choices given dvConf conditioning:
                nCW_H(iO)       = sum(choice(splitInd==1 & oriInd)==1);
                nCW_L(iO)       = sum(choice(splitInd==0 & oriInd)==1);
                
                nCCW_H(iO)      = sum(choice(splitInd==1  & oriInd)==0);
                nCCW_L(iO)      = sum(choice(splitInd==0  & oriInd)==0);
            end
            
            nHigh(iC,iDate)     = sum((nCW_H+nCCW_H));
            nLow(iC,iDate)      = sum((nCW_L+nCCW_L));
            propCW_H            = nCW_H./(nCW_H+nCCW_H);
            propCW_L            = nCW_L./(nCW_L+nCCW_L);

            paramsSplit                       = fitPF(uOri,  nCCW_H, nCW_H,nCCW_L, nCW_L);
            ratioMedianSplit(iC,iDate)        = (paramsSplit(2))./(paramsSplit(3)); % high conf: low conf
            
            HC(iC,iDate)                      = paramsSplit(2);
            LC(iC,iDate)                      = paramsSplit(3);
            
            if iDate==exampleDateMedSplit && iC==exampleContrast
                propCW_H_plot                 = propCW_H ;
                propCW_L_plot                 = propCW_L ;
                nPointH                       = nCW_H+nCCW_H;
                nPointL                       = nCW_L+nCCW_L;
                paramPlot                     = paramsSplit;
                oriPlot                       = uOri;
                stats.fig3A.nTrials           = sum(nCW_H+nCCW_H +nCW_L+nCCW_L); 
            end
            
            % Clear vars:
            nCW_L       = [];
            nCCW_L      = [];
            nCW_H       = [];
            nCCW_H      = [];
            uOri        = [];
            propCW      = [];
            nH          = [];
            nL          = [];
        end
    end 
end
%% 3.A Example psychometric function- partioned by confidence median split

figure;hold on; axis square;

oriVals         = linspace(min(oriPlot),max(oriPlot),200);
predH_PF        = paramPlot(1)+ (1 - 2*paramPlot(1))*normcdf((oriVals - paramPlot(4))/(2*paramPlot(2)));
predL_PF        = paramPlot(1)+ (1 - 2*paramPlot(1))*normcdf((oriVals - paramPlot(4))/(2*paramPlot(3)));

plot(oriVals,predH_PF,'g-','linewidth',3)
plot(oriVals,predL_PF,'r-','linewidth',3)
for iO=1:length(oriPlot)
    plot(  oriPlot(iO),    propCW_H_plot(iO) ,'ko','markerfacecolor','g','markersize',round(nPointH(iO)/2))
    plot(  oriPlot(iO),    propCW_L_plot(iO) ,'ko','markerfacecolor','r','markersize',round(nPointL(iO)/2))
end
title(['Example date: ' num2str(exampleDate) '; contrast = ' num2str(exampleContrast)],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
xlabel('Orientation','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
ylabel('Proportion CW','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);

% Stat 3A
stats.fig3A.exSlopeRatio = ratioMedianSplit(exampleContrast,exampleDate);

%% 3.B Summary slope ratio
% Group low and high contrast confidence split difference
uncRatio_highCont        = (ratioMedianSplit(1,:)); % High:Low conf. for High contrast stim
uncRatio_lowCont         = (ratioMedianSplit(2,:)); % High:Low conf. for low contrast stim

uncRatio_F               = [uncRatio_highCont(obsInd) uncRatio_lowCont(obsInd)] ;
uncRatio_Z               = [uncRatio_highCont(~obsInd)  uncRatio_lowCont(~obsInd)];
comb_FZ                  = [ uncRatio_F  uncRatio_Z]; % combine two observers

% Calculate percentiles:
z_uncRatio_p             = prctile(log10(uncRatio_Z),[25 50 75]);
f_uncRatio_p             = prctile(log10(uncRatio_F),[25 50 75]);
comb_uncRatio_p          = prctile(log10(comb_FZ) ,[25 50 75]);

figure;hold on; axis square
xlim([0 4])
ylim([-1 1])
plot([0 4],[0 0])
plot(0+ones(length(uncRatio_Z),1)+(rand(length(uncRatio_Z),1)-.5)/5,log10(uncRatio_Z),'o','markerfacecolor','r','markeredgecolor','k','markersize',10) % Ziggy
plot(1+ones(length(uncRatio_F),1)+(rand(length(uncRatio_F),1)-.5)/5, log10(uncRatio_F),'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Friedrich
plot(2+ones(length( comb_FZ),1)+(rand(length(     comb_FZ),1)-.5)/5, log10(comb_FZ),'o','markerfacecolor','g','markeredgecolor','k','markersize',10) % Friedrich

plot([1 1]+0,z_uncRatio_p([1 3]))
plot([1 1]+1,f_uncRatio_p([1 3]))
plot([1 1]+2,comb_uncRatio_p([1 3]))

plot([1 ]+0,z_uncRatio_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,f_uncRatio_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,comb_uncRatio_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Observer','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Ratio High to Low contrast','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({ 'Ziggy', 'Friedrich','Combined'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Log 10 slope ratio','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

% Stats 3B
[p_f,h_f]                   = signrank(1 - uncRatio_F);
[p_z,h_z]                   = signrank(1 - uncRatio_Z);
[p_c,h_z]                   = signrank(1 - comb_FZ);


stats.fig3B.diffOne_f       = p_f;
stats.fig3B.diffOne_z       = p_z;
stats.fig3B.diffOne_comb    = p_c;

stats.fig3B.med_ratio_f     = nanmedian(uncRatio_F);
stats.fig3B.med_ratio_z     = nanmedian( uncRatio_Z);
stats.fig3B.med_ratio_comb  = nanmedian(comb_FZ);

%% 3.C Effect of stimulus strength and gain

colorsTemp      = linspace(0,1,11);
colors          = [fliplr(colorsTemp)' zeros(11,1) colorsTemp'];
gainColors      = flipud(repmat(linspace(0,.8,7),3,1)');

for iDate = 1:29
    exDate                  = iDate;
    dateText                = masterNeuralStruct{exDate}.expt.date;
    
    choice                  = [];
    confidence              = [];
    choiceCorrected         = [];
    confidenceCorrected     = [];
    values                  = [];
    
    choice                  = readtable([CSVpath '/Prediction/Choice' dateText '_sim_gain_score4_wRepeat_muScore4.csv']);
    confidence              = readtable([CSVpath '/Prediction/Confidence' dateText '_sim_gain_score4_wRepeat_muScore4.csv']);
    values                  = readtable([CSVpath '/' dateText '_sim_gain_score4.csv']);
    
    choice                  = table2array(choice);
    confidence              = table2array(confidence);
    
    values                  = table2array(values);
    nTrial                  = length(values);
    
    incInd                  = 1:length(values);
    gainVec                 = values(incInd,4);
    oriVec                  = values(incInd,7);
    
    oriVal                  = unique(oriVec);
    gain                    = unique(gainVec);
    
    predChoice              = choice(incInd,4);
    predConf                = confidence(incInd,4);
    
    % DV in probability:
    dvChoice                = choice(incInd,5);
    dvConf                  = confidence(incInd,5);
    
    % Norm inverse of probability:
    dvNIChoice                      = norminv(choice(incInd,5)); % take norminv immediatly
    dvNIConf                        = norminv(confidence(incInd,5));
    uOri                            = unique(oriVal);
    sori(iDate)                     = length(uOri);
    largeOriInd                     = oriVec == uOri(1) |oriVec == uOri(2) | oriVec == uOri(10) | oriVec == uOri(11);
    smallOriInd                     = oriVec == uOri(4) | oriVec == uOri(5) | oriVec == uOri(7) | oriVec == uOri(8);
    
    gain_one                        = gainVec==1;
    gain_two                        = gainVec==2.0;
    gain_half                       = gainVec==0.5;
    
    stimStrengthDiff_HC(iDate)      = mean(dvConf(gain_one & largeOriInd))-mean(dvConf(gain_one & smallOriInd));
    contrastDiff_HC(iDate)          = mean(dvConf(gain_two ))-mean(dvConf(gain_half ));
    
    stimStrengthLarge(iDate)        = mean(dvConf(gain_one & largeOriInd)) ; % large stim strength
    stimStrengthSmall(iDate)        = mean(dvConf(gain_one & smallOriInd));  % small stim strength
    
    gainLarge(iDate)                = mean(dvConf(gain_two )) ; % large gain
    gainSmall(iDate)                = mean(dvConf(gain_half )); % small gain;
    
    if iDate==exampleDate
        figure; box off; axis square; hold on
        phyplot([], [], 'k-', 'xticks', [min(oriVec) 0 max(oriVec)], 'yticks', [0:.5:1], 'linewidth', 1, 'width', 1, 'fontsize', 20);
        hold on;
        for iG=1:length(gain)
            oris                        = oriVec(gainVec==gain(iG));
            confPlot                    = dvConf(gainVec==gain(iG));
            crossDateConf(:,iG,iDate)   = dvNIConf(gainVec==gain(iG));
            p(iG)                       = plot(oris, confPlot ,'-','color',gainColors(iG,:),'linewidth',2) ;
            leg{iG}                     = ['Gain  ' num2str(gain(iG))];
        end
        legend(p,leg,'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16 )
        xlabel('Orientation','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
        ylabel('Avg probability HC ','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
    end
end


%% Calculate percentiles:
stimStrength_p          = prctile( stimStrengthDiff_HC,[25 50 75]);
contrastEffect_p        = prctile(contrastDiff_HC,[25 50 75]);

figure;hold on; axis square
xlim([0 3])
ylim([-.5 1])
plot(ones(length( stimStrengthDiff_HC(obsInd)),1)+(rand(length( stimStrengthDiff_HC(obsInd)),1)-.5)/5, stimStrengthDiff_HC(obsInd),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(ones(length( stimStrengthDiff_HC(~obsInd)),1)+(rand(length( stimStrengthDiff_HC(~obsInd)),1)-.5)/5, stimStrengthDiff_HC(~obsInd),'o','markerfacecolor','b','markeredgecolor','k','markersize',10)
plot(1+ones(length(contrastDiff_HC(obsInd)),1)+(rand(length(contrastDiff_HC(obsInd)),1)-.5)/5,contrastDiff_HC(obsInd),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)  % F
plot(1+ones(length(contrastDiff_HC(~obsInd)),1)+(rand(length(contrastDiff_HC(~obsInd)),1)-.5)/5,contrastDiff_HC(~obsInd),'o','markerfacecolor','b','markeredgecolor','k','markersize',10) % Z

plot([1 1],stimStrength_p([1 3]))
plot([1 1]+1,contrastEffect_p([1 3]))

plot([1 ],stimStrength_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,contrastEffect_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

xlabel('Stat','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({'Stim strength',  'Contrast (gain) effect'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Diff % HC','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)


% Stats
[p_ss,h_ss]                         = signrank(stimStrengthDiff_HC);
[p_c,h_c]                           = signrank(contrastDiff_HC);

stats.fig3D.stat_stimStregnth       = p_ss;
stats.fig3D.statGain                = p_c;
stats.fig3D.avgStimStregnth_diff    = median(stimStrengthDiff_HC);
stats.fig3D.avgGain_diff            = median(contrastDiff_HC);


%% Stim strength and gain effect version 2 (connected lines)
stimStrengthLarge_p          = prctile( stimStrengthLarge,[25 50 75]);
stimStrengthSmall_p          = prctile( stimStrengthSmall,[25 50 75]);

gainLarge_p                  = prctile( gainLarge,[25 50 75]);
gainSmall_p                  = prctile( gainSmall,[25 50 75]);

obsColor        = {'y','b'};

figure;hold on; axis square
xlim([0 6])
ylim([-.1 1])

for iO=1:length(stimStrengthLarge  )
    % stim strength:
    randVal             = (rand(2,1)-.5)/5;
    plot(1+randVal(1),  [stimStrengthLarge(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
    plot(2+randVal(2),  [stimStrengthSmall(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
    plot([1+randVal(1) 2+randVal(2)],  [stimStrengthLarge(iO) stimStrengthSmall(iO)],'-','color',obsColor{obsInd(iO)+1})
    
    % gain:
    randVal             = (rand(2,1)-.5)/5;
    plot(4+randVal(1),  [gainLarge(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
    plot(5+randVal(2),  [gainSmall(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
    plot([4+randVal(1) 5+randVal(2)],  [gainLarge(iO) gainSmall(iO)],'-','color',obsColor{obsInd(iO)+1})
    
end

plot([1 1],stimStrengthLarge_p ([1 3]))
plot([1 1]+1,stimStrengthSmall_p ([1 3]))

plot([1 1]+3,gainLarge_p ([1 3]))
plot([1 1]+4,gainSmall_p ([1 3]))


xlabel('Stat','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3 4 5])
xticklabels({'Stim strength large',  'small', '','Gain high' ,'gain low'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
ylabel('Stat val','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

% Stads fig 3D
[p_val_stimStrength,h]              = ranksum(stimStrengthLarge,stimStrengthSmall);
[p_val_gain,h]                      = ranksum(gainLarge,gainSmall);

stats.fig3D.p_val_stimStrength      = p_val_stimStrength;
stats.fig3D.avg_SS_diff             = median(stimStrengthLarge-stimStrengthSmall);
stats.fig3D.p_val_gain              = p_val_gain;
stats.fig3D.avg_gain_diff           = median(gainLarge-gainSmall);



%% Percent correctly predicted choices and confidence
NL = load(strcat([structurePath  '/decodeSummary_NL.mat']));

for iDate=1:29
    pChoiceCorrect(iDate)   = mean(NL.popStruct{iDate}.predChoice== NL.popStruct{iDate}.obsChoice);
    pConfCorrect(iDate)     = mean(NL.popStruct{iDate}.predConf== NL.popStruct{iDate}.obsConf);
end

stat.fig3.avgChoiceCorrect      = median(pChoiceCorrect); 
stat.fig3.avgConfCorrect        = median(pConfCorrect);

%% Functions
function [P] = fitPF(oriVals, nCCW_H, nCW_H,nCCW_L, nCW_L) % Fit psychometric function
options = optimset('Display', 'off', 'Maxiter', 1000, 'MaxFuneval', 2000, 'Algorithm', 'interior-point');

%% FitPF
lBounds(1) = 0.0000;             uBounds(1)    = .02;                         % [lapse rate]
lBounds(2) = 0.0001;             uBounds(2)    = 10;                          % [internal noise High conf]
lBounds(3) = 0.0001;             uBounds(3)    = 10;                          % [internal noise Low conf]
lBounds(4) = -10;                uBounds(4)    = 10;                          % [criterion]

startValues = [0 , .1, .1, 0];

obFun       = @(params) getPFNLL(params, oriVals, nCCW_H, nCW_H,nCCW_L, nCW_L);
P           = fmincon(obFun, startValues, [], [], [], [], lBounds, uBounds, [], options);

end

function [NLL] = getPFNLL(params, oriVals, nCCW_H, nCW_H,nCCW_L, nCW_L) % Get negative log-likelihood of PF fit

noiseInt_HC  = params(2);
noiseInt_LC  = params(3);
lapseRate    = params(1);
stimCrit     = params(4);

propCW_H     = lapseRate + (1 - 2*lapseRate)*normcdf((oriVals - stimCrit)/(2*noiseInt_HC));
propCW_L     = lapseRate + (1 - 2*lapseRate)*normcdf((oriVals - stimCrit)/(2*noiseInt_LC));

NLL_HC       = -sum(log(binopdf(nCW_H', (nCW_H + nCCW_H)', propCW_H)));
NLL_LC       = -sum(log(binopdf(nCW_L', (nCW_L + nCCW_L)', propCW_L)));

NLL          = NLL_HC + NLL_LC;
end
