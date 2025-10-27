%% Figure 4
% Written 05-08-2024 ZBS (zbsinger@mit.edu)

%%
% 4.A Example pop avg choice DV split by choice at high and low cotrast
% 4.B Population choice probability 
% 4.C Example pop avg confidence DV split by choice at high and low cotrast
% 4.D Slope ratio and offset difference cross population
% 4.E Population confidence choice probability 
% 4.F Offset difference (contrast and confidence) cross population

%% Start with clean slate
clearvars -except masterNeuralStruct
clc;
close all;

% Set Paths
user                = '/Users/zbsinger/MIT Dropbox/Zoe Boundy-Singer/Projects/V1_Confidence_Physiology';
thisPath            = user;
structurePath       = strcat([thisPath, '/Structures/Neural']);
functionPath        = strcat([thisPath, '/Functions/']);
addpath(thisPath,functionPath,structurePath )

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
exampleDate     = 7;

%% Plotting conventions:
pCol                    = [0 0 0; 153 0 6]./255; % HC
pCol(:,:,2)             = [170 170 170; 255 153 204]./255; % LC
cCol                    = [0 0 0; .9 .9 .9]; % contrast color
cCol2                   = [ 0 1 0;1 0 0;];
cCol2(:,:,2)            = [ 177, 224, 187;230, 154, 164]./255;
contColor               = [0 0 0; .5 .5 .5]; % contrast color


%% Presaved decoding sumary:
NL = load(strcat([structurePath  '/decodeSummary_NL.mat']));
%% 4.S Choice DV for example pop split by choice (high and low contrast)

OPTIONS.Display             = 'off';

%% Descriptive fit key: 
%% 1.1 Choice: effect of contrast on slpoe
%% 1.2 Choice: effect of perceptual decision on offset
%% 2.1 Confidence: effect of contrast on offset
%% 2.2 Confidence: effect of confidence on offset

for iDate = 1:29
    %% Section 1: Choice
    
    for iC=1:2
        y11{iC}                                 = NL.dvChoicePerOri{iC,iDate};
        nObs11{iC}                              = NL.NdvChoicePerOri{iC,iDate};
        oriVal11{iC}                            = NL.oriVal{iC,iDate};
        
        if iC==1 % high contrast
            nObs12(:,1)                           = NL.NdvChoicePerOri_CW{1,iDate};
            nObs12(:,2)                           = NL.NdvChoicePerOri_CCW{1,iDate};
            y12{1}                                = NL.dvChoicePerOri_CW{1,iDate};
            y12{2}                                = NL.dvChoicePerOri_CCW{1,iDate};
        end
        
        nObs12_2=[];
        if iC==2 % low contrast
            nObs12_2(:,1)                           = NL.NdvChoicePerOri_CW{2,iDate};
            nObs12_2(:,2)                           = NL.NdvChoicePerOri_CCW{2,iDate};
            y12_2{1}                                = NL.dvChoicePerOri_CW{2,iDate};
            y12_2{2}                                = NL.dvChoicePerOri_CCW{2,iDate};
        end
    end

    %% 1.1 Effect of contrast
    obFunLinear11                           = @(paramsLinear11) getLinearDiffSlope(paramsLinear11, oriVal11,   y11, nObs11 );
    [paramLinear11]                         = fmincon(obFunLinear11, [1 1 0  ], [], [], [], [], [0 0 -1  ],[100 100 100 ],[],OPTIONS);
    crossDateParam11(:,iDate)               = paramLinear11;% effect of contrast on slope
    
    %% 1.2 Linear: effect choice on offset (free or fixed) - high contrast
    % linear regression fit- joint with diff offset
    obFunLinearJoint                        = @(paramsLinearJoint) getLinearDiffOffset(paramsLinearJoint, NL.oriVal{1,iDate},   y12,nObs12 );
    [paramLinearJoint,MSE]                  = fmincon(obFunLinearJoint, [1 0 0 ], [], [], [], [], [0 -1 -1],[100 100 100],[],OPTIONS);
    crossDateParam12(:,iDate,1)             = paramLinearJoint; % effect of choice on slope
    
    %% 1.2 low contrast
    obFunLinearJoint                        = @(paramsLinearJoint) getLinearDiffOffset(paramsLinearJoint, NL.oriVal{2,iDate},   y12_2,nObs12_2 );
    [paramLinearJoint,MSE]                  = fmincon(obFunLinearJoint, [1 0 0 ], [], [], [], [], [0 -1 -1],[100 100 100],[],OPTIONS);
    crossDateParam12(:,iDate,2)             = paramLinearJoint; % effect of choice on slope
    
    %% Section 2: Confidence
    for iC=1:2
        oriVal21{iC}                        = NL.oriVal{iC,iDate};
        dvConfQuad21{iC}                    = NL.dvConfPerOri{iC,iDate};
        Nobs21{iC}                          = NL.NdvConfPerOri{iC,iDate};
    end
    
    %% 2.2 Quad - joint (diff offset) effect of confidence
    Nobs22                                  = [];
    Nobs22{1}                               = NL.NdvConfPerOri_HC{1,iDate};
    Nobs22{2}                               = NL.NdvConfPerOri_LC{1,iDate};
    
    % High contrast
    oriQuad22{1}                            = NL.oriVal{1,iDate};
    oriQuad22{2}                            = NL.oriVal{1,iDate};
    yConf22{1}                              = NL.dvConfPerOri_HC{1,iDate};
    yConf22{2}                              = NL.dvConfPerOri_LC{1,iDate};
    
    %Low contrast:
    Nobs22_2                                = [];
    Nobs22_2{1}                             = NL.NdvConfPerOri_HC{2,iDate};
    Nobs22_2{2}                             = NL.NdvConfPerOri_LC{2,iDate};
    
    oriQuad22_2{1}                          = NL.oriVal{2,iDate};
    oriQuad22_2{2}                          = NL.oriVal{2,iDate};
    yConf22_2{1}                            = NL.dvConfPerOri_HC{2,iDate};
    yConf22_2{2}                            = NL.dvConfPerOri_LC{2,iDate};
    
    
    %% 2.1 Contrast effect
    obFunQuad                               = @(paramsQuad) getQuadJointDiffOff(paramsQuad, oriVal21,   dvConfQuad21, Nobs21 );
    [paramQuad]                             = fmincon(obFunQuad, [1 0 0 0], [], [], [], [], [0 -1 -1 -1 ]*5,[100 100 100 100],[],OPTIONS);
    crossDateParam21(:,iDate)               = paramQuad; % compare param 1 vs param 2
    
    %% 2.1.5 Contrast effect - independent quad fit- high contrast
    obFunQuad                               = @(paramsQuad) getQuadJointDiffOff(paramsQuad, oriVal21(:,1),   dvConfQuad21(1), Nobs21(1) );
    [paramQuad]                             = fmincon(obFunQuad, [1 0 0 ], [], [], [], [], [-1 -1 -1  ]*5,[100 100 100 ],[],OPTIONS);
    crossDateParam215(:,iDate,1)            = paramQuad;
    
    % low contrast
    obFunQuad                               = @(paramsQuad) getQuadJointDiffOff(paramsQuad, oriVal21(:,2),   dvConfQuad21(2), Nobs21(2) );
    [paramQuad]                             = fmincon(obFunQuad, [1 0 0 ], [], [], [], [], [-1 -1 -1  ]*5,[100 100 100 ],[],OPTIONS);
    crossDateParam215(:,iDate,2)            = paramQuad;
    
    %% 2.2 Confidence effect- high contrast
    obFunQuadJoint                          = @(paramsQuadJoint) getQuadJointDiffOff(paramsQuadJoint, oriQuad22,   yConf22 ,Nobs22 );
    [paramsQuadJoint]                       = fmincon(obFunQuadJoint, [1 0 0 0], [], [], [], [], [0 -1 -1 -1 ]*5,[100 100 100 100],[],OPTIONS);
    crossDateParam22(:,iDate,1)             = paramsQuadJoint;
    
    %% 2.2 Low contrast
    obFunQuadJoint                          = @(paramsQuadJoint) getQuadJointDiffOff(paramsQuadJoint, oriQuad22_2,   yConf22_2 ,Nobs22_2 );
    [paramsQuadJoint]                       = fmincon(obFunQuadJoint, [1 0 0 0], [], [], [], [], [0 -1 -1 -1 ]*5,[100 100 100 100],[],OPTIONS);
    crossDateParam22(:,iDate,2)             = paramsQuadJoint;
end

%% Slope statistic:
slope_cont_ratio_11 = crossDateParam11(2,:)./ crossDateParam11(1,:); % low: high contrast slope

%% Offset statistic
offset_choice_diff_12 = crossDateParam12(2,:,:) - crossDateParam12(3,:,:); % CCW minus CW offset

%% Confidence: contrast Offset:
offset_cont_diff_21 = crossDateParam21(3,:) - crossDateParam21(4,:); % high minus low contrast slope

%% Confidence Offset confidene statistic
offset_conf_diff_22 = crossDateParam22(3,:,:) - crossDateParam22(4,:,:); % HC minus LC offset

%% Plot example population:
yxaisObs{1,1}       = [-3:1.5:3]; % this axis range is tailored to example population
yxaisObs{1,2}       = [ -1:1:2 ];
xax                 = linspace(.01,1,4);
xaxis               = [-20:10:20]; % stimulus axis

%% 4.X Plot choice DV split by contrast (not used)
figure; axis square; box off; hold on;
phyplot([], [], 'k-', 'xticks', xaxis, 'yticks', yxaisObs{1,1}, 'linewidth', 1, 'width', 1, 'fontsize', 20);
hold on;
ylim([-5 5])

for iC = 1:2
    markerS                 = NL.NdvChoicePerOri{iC,exampleDate}./max( NL.NdvChoicePerOri{iC,exampleDate})*10+1; % marker size tailored to # observations
    for iS=1:length(NL.oriVal{iC,exampleDate})
        plot(NL.oriVal{iC,exampleDate}(iS),NL.dvChoicePerOri{iC,exampleDate}(iS) ,'ko','color','w','markerfacecolor',pCol(1,:,iC),'linewidth',1,'markersize',round(markerS(iS)));
    end
    oriX                    = [NL.oriVal{iC,exampleDate}(1):.1:NL.oriVal{iC,exampleDate}(end)];
    pl1(iC)                 = plot(oriX,[oriX.*crossDateParam11(iC,exampleDate)+crossDateParam11(3,exampleDate)], '-','color',pCol(1,:,iC),'linewidth',3);
    
end
legend(pl1,{'High contrast','Low contrast'},'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20,'location','southeast');
title(['Choice DV;' masterNeuralStruct{exampleDate}.expt.date],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
xlabel('Orientation ','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
ylabel(['Choice DV'],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);

%% 4.A Plot high contrast choice DV split by choice
contText= {'High','Low'};
for iC=1:2
    figure; axis square; box off; hold on;
    phyplot([], [], 'k-', 'xticks', xaxis, 'yticks', yxaisObs{1,1}, 'linewidth', 1, 'width', 1, 'fontsize', 20);
    xlim([xaxis(1) xaxis(end)])
    ylim([-5 5])
    hold on;
    
    markerS_CW      = NL.NdvChoicePerOri_CW{iC,exampleDate}./max( [NL.NdvChoicePerOri_CW{iC,exampleDate}  NL.NdvChoicePerOri_CCW{iC,exampleDate}])*10+1;
    markerS_CCW     = NL.NdvChoicePerOri_CCW{iC,exampleDate}./max( [NL.NdvChoicePerOri_CW{iC,exampleDate}  NL.NdvChoicePerOri_CCW{iC,exampleDate}])*10+1;
    nObs(:,1)       = NL.NdvChoicePerOri_CW{iC,exampleDate};
    nObs(:,2)       = NL.NdvChoicePerOri_CCW{iC,exampleDate};
    
    for iS = 1:length(NL.oriVal{iC,exampleDate})
        plot(NL.oriVal{iC,exampleDate}(iS),NL.dvChoicePerOri_CW{iC,exampleDate}(iS) ,'ko','color','w','markerfacecolor',pCol(1,:,iC),'linewidth',1,'markersize',round(markerS_CW(iS)));
        plot(NL.oriVal{iC,exampleDate}(iS),NL.dvChoicePerOri_CCW{iC,exampleDate}(iS) ,'ko','color','w','markerfacecolor',pCol(2,:,iC),'linewidth',1,'markersize',round(markerS_CCW(iS)));
    end
    y{1}            = NL.dvChoicePerOri_CW{iC,exampleDate};
    y{2}            = NL.dvChoicePerOri_CCW{iC,exampleDate};
    oriX            = [NL.oriVal{2,exampleDate}(1):.1:NL.oriVal{2,exampleDate}(end)];
    
    if iC==1
        pl2(1)          = plot(oriX,[oriX.*crossDateParam12(1,exampleDate,1)+crossDateParam12(2,exampleDate,1)], '-','color',pCol(1,:,iC),'linewidth',3);
        pl2(2)          = plot(oriX,[oriX.*crossDateParam12(1,exampleDate,1)+crossDateParam12(3,exampleDate,1)], '-','color',pCol(2,:,iC),'linewidth',3);
    elseif iC==2
        pl2(1)          = plot(oriX,[oriX.*crossDateParam12(1,exampleDate,2)+crossDateParam12(2,exampleDate,2)], '-','color',pCol(1,:,iC),'linewidth',3);
        pl2(2)          = plot(oriX,[oriX.*crossDateParam12(1,exampleDate,2)+crossDateParam12(3,exampleDate,2)], '-','color',pCol(2,:,iC),'linewidth',3);
    end
    legend(pl2,{'CW','CCW choice'},'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20,'location','southeast');
    title(['Contrast = ', contText{iC} ' ;Choice DV;' masterNeuralStruct{exampleDate}.expt.date],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
    xlabel('Orientation ','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
    ylabel(['Choice DV high contrast split by choice'],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
end

% Stats 4A
stats.fig4A.highContSlope_all_ex    = crossDateParam11(1,exampleDate); 
stats.fig4A.lowContSlope_all_ex     = crossDateParam11(2,exampleDate);
stats.fig4A.highCont_slope_ex       = crossDateParam12(1,exampleDate,1); % as plotted in fig4A (i.e slope when fit to choice split)
stats.fig4A.lowCont_slope_ex        = crossDateParam12(1,exampleDate,2);
stats.fig4A.highCont_offset_diff    = offset_choice_diff_12(1,exampleDate,1);
stats.fig4A.lowCont_offset_ex       = offset_choice_diff_12(1,exampleDate,2);

%% 4.X: Confidence varaible split by contrast (cont used)
figure; axis square; box off; hold on;
phyplot([], [], 'k-', 'xticks', xaxis, 'yticks',  yxaisObs{1,2}, 'linewidth', 1, 'width', 1, 'fontsize', 20);
hold on;
ylim([-2 3])
for iC = 1:2
    markerS = NL.NdvConfPerOri{iC,exampleDate}./max( [NL.NdvConfPerOri{iC,exampleDate} ])*10+1;
    for iS = 1:length(NL.oriVal{iC,exampleDate})
        plot(NL.oriVal{iC,exampleDate}(iS), NL.dvConfPerOri{iC,exampleDate}(iS)  ,'o','color','w','markerfacecolor',pCol(1,:,iC),'linewidth',1,'markersize',round(markerS(iS)));
    end
    predQuad        = crossDateParam21(1,exampleDate)*(oriX-crossDateParam21(2,exampleDate)).^2 +crossDateParam21(2+iC,exampleDate);
    pl3(iC)         = plot(oriX, predQuad ,'-','color',pCol(1,:,iC),'linewidth',3);
end

title(['Confidence DV'],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
xlabel('Orientation ','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
ylabel('Confidence DV','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);

%% 4C: Confidence DV split by confidence
contText= {'High','Low'};
for iC=1:2
    figure; axis square; box off; hold on;
    phyplot([], [], 'k-', 'xticks', xaxis, 'yticks', [-2:1:3], 'linewidth', 1, 'width', 1, 'fontsize', 20);
    xlim([xaxis(1) xaxis(end)])
    ylim([-2 3])
    hold on;

    markerS_HC      = NL.NdvConfPerOri_HC{iC,exampleDate}./max( [NL.NdvConfPerOri_LC{iC,exampleDate}  NL.NdvConfPerOri_HC{iC,exampleDate}])*10+1;
    markerS_LC      = NL.NdvConfPerOri_LC{iC,exampleDate}./max( [NL.NdvConfPerOri_LC{iC,exampleDate}  NL.NdvConfPerOri_HC{iC,exampleDate}])*10+1;
    nObs(:,1)       = NL.NdvConfPerOri_HC{iC,exampleDate};
    nObs(:,2)       = NL.NdvConfPerOri_LC{iC,exampleDate};
    
    for iS = 1:length(NL.oriVal{iC,exampleDate})
        plot(NL.oriVal{iC,exampleDate}(iS),NL.dvConfPerOri_HC{iC,exampleDate}(iS) ,'ko','color','w','markerfacecolor',pCol(1,:,iC),'linewidth',1,'markersize',round(markerS_HC(iS)));
        plot(NL.oriVal{iC,exampleDate}(iS),NL.dvConfPerOri_LC{iC,exampleDate}(iS) ,'ko','color','w','markerfacecolor',pCol(2,:,iC),'linewidth',1,'markersize',round(markerS_LC(iS)));
    end
    y{1}            = NL.dvConfPerOri_HC{iC,exampleDate};
    y{2}            = NL.dvConfPerOri_LC{iC,exampleDate};
    oriX            = [NL.oriVal{2,exampleDate}(1):.1:NL.oriVal{2,exampleDate}(end)];
    
    if iC==1
        pl2(1)          = plot(oriX,[ crossDateParam22(1,exampleDate,1)*(oriX-crossDateParam22(2,exampleDate,1)).^2 + crossDateParam22(3,exampleDate,1)], '-','color',pCol(1,:,iC),'linewidth',3);
        pl2(2)          = plot(oriX,[ crossDateParam22(1,exampleDate,1)*(oriX-crossDateParam22(2,exampleDate,1)).^2 + crossDateParam22(4,exampleDate,1)], '-','color',pCol(2,:,iC),'linewidth',3);
    elseif iC==2
        pl2(1)          = plot(oriX,[ crossDateParam22(1,exampleDate,2)*(oriX-crossDateParam22(2,exampleDate,2)).^2 + crossDateParam22(3,exampleDate,2)], '-','color',pCol(1,:,iC),'linewidth',3);
        pl2(2)          = plot(oriX,[ crossDateParam22(1,exampleDate,2)*(oriX-crossDateParam22(2,exampleDate,2)).^2 + crossDateParam22(4,exampleDate,2)], '-','color',pCol(2,:,iC),'linewidth',3);
    end
    legend(pl2,{'HC','CCW'},'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20,'location','southeast');
    title(['Contrast = ', contText{iC} ' ;Conf DV;' masterNeuralStruct{exampleDate}.expt.date],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
    xlabel('Orientation ','fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
    ylabel(['Confidence DV high contrast split by conf'],'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 20);
end

% Stats 4C
stats.fig4C.highContQuad_all_ex    = crossDateParam215(:,exampleDate,1)    ;
stats.fig4C.lowContQuad_all_ex     = crossDateParam215(:,exampleDate,2) ;   
stats.fig4C.highCont_Quad_ex       = crossDateParam22(1,exampleDate,1); % as plotted in fig4A (i.e slope when fit to choice split)
stats.fig4C.lowCont_Quad_ex        = crossDateParam22(1,exampleDate,2);
stats.fig4C.cont_offset_diff       = offset_cont_diff_21(1,exampleDate,1);
stats.fig4C.highConf_offset_diff   = offset_conf_diff_22(1,exampleDate,1);
stats.fig4C.lowConf_offset_ex      = offset_conf_diff_22(1,exampleDate,2);

%% 4B Cross pop sumary
% Slope greater than 0?

highContSlope_11        = crossDateParam11(1,:);
lowContSlope_11         = crossDateParam11(2,:);
highContSlope_11_p      = prctile(highContSlope_11,[25 50 75]);
lowContSlope_11_p       = prctile(lowContSlope_11 ,[25 50 75]);

obsColor                = {'y','b'};

figure;hold on; axis square
xlim([0 3])
ylim([-.3 .3])
plot([0 3],[0 0])

for iO=1:length(highContSlope_11)
randVal             = (rand(2,1)-.5)/5;
plot(1+randVal(1),  [highContSlope_11(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot(2+randVal(2),  [lowContSlope_11(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot([1+randVal(1) 2+randVal(2)],  [highContSlope_11(iO) lowContSlope_11(iO)],'-','color',obsColor{obsInd(iO)+1})
end

plot([1 1],highContSlope_11_p([1 3]))
plot([1 ],highContSlope_11_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

plot([1 1]+1,lowContSlope_11_p([1 3]))
plot([1 ]+1,lowContSlope_11_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

ylabel('Slope','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 ])
xticklabels({'High cont.','Low cont'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)


% Quad term greater than 0?
highContQuad_215        = crossDateParam215(1,:,1); % high cont quad term
lowContQuad_215         = crossDateParam215(1,:,2); % low cont quad term
highContQuad_215_p      = prctile(highContQuad_215,[25 50 75]);
lowContQuad_215_p       = prctile(lowContQuad_215 ,[25 50 75]);


figure;hold on; axis square
xlim([0 3])
ylim([-.02 .02])
plot([0 3],[0 0])

for iO=1:length(highContQuad_215  )
randVal             = (rand(2,1)-.5)/5;
plot(1+randVal(1),  [highContQuad_215(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot(2+randVal(2),  [lowContQuad_215(iO)],'o','markerfacecolor',obsColor{obsInd(iO)+1},'markeredgecolor','k','markersize',10)
plot([1+randVal(1) 2+randVal(2)],  [highContQuad_215(iO) lowContQuad_215(iO)],'-','color',obsColor{obsInd(iO)+1})
end

plot([1 1],highContQuad_215_p([1 3]))
plot([1 ],highContQuad_215_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

plot([1 1]+1,lowContQuad_215_p([1 3]))
plot([1 ]+1,lowContQuad_215_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

ylabel('Quad term','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 ])
xticklabels({'High cont.','Low cont'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)


% Offset diff
offsetDiff_12_p         = prctile( offset_choice_diff_12(:),[25 50 75]);
offsetDiff_21_p         = prctile( offset_cont_diff_21,[25 50 75]);
offsetDiff_22_p         = prctile( offset_conf_diff_22(:),[25 50 75]);

offset_choice_diff_12_f = offset_choice_diff_12(:,obsInd,:);
offset_choice_diff_12_z = offset_choice_diff_12(:,~obsInd,:);

offset_cont_diff_21_f   = offset_cont_diff_21(obsInd);
offset_cont_diff_21_z   = offset_cont_diff_21(~obsInd);

offset_conf_diff_22_f   = offset_conf_diff_22(:,obsInd,:);
offset_conf_diff_22_z   = offset_conf_diff_22(:,~obsInd,:);

figure;hold on; axis square
xlim([0 4])
ylim([-.5 1])
plot([0 4],[0 0])
plot(ones(length( offset_choice_diff_12_f(:)),1)+(rand(length(  offset_choice_diff_12_f(:)),1)-.5)/5,  offset_choice_diff_12_f(:),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(ones(length( offset_choice_diff_12_z(:)),1)+(rand(length(  offset_choice_diff_12_z(:)),1)-.5)/5,  offset_choice_diff_12_z(:),'o','markerfacecolor','b','markeredgecolor','k','markersize',10)

plot(1+ones(length( offset_cont_diff_21_f(:)),1)+(rand(length(  offset_cont_diff_21_f(:)),1)-.5)/5,  offset_cont_diff_21_f(:),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length( offset_cont_diff_21_z(:)),1)+(rand(length(  offset_cont_diff_21_z(:)),1)-.5)/5,  offset_cont_diff_21_z(:),'o','markerfacecolor','b','markeredgecolor','k','markersize',10)

plot(2+ones(length( offset_conf_diff_22_f(:)),1)+(rand(length(  offset_conf_diff_22_f(:)),1)-.5)/5,  offset_conf_diff_22_f(:),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(2+ones(length( offset_conf_diff_22_z(:)),1)+(rand(length(  offset_conf_diff_22_z(:)),1)-.5)/5,  offset_conf_diff_22_z(:),'o','markerfacecolor','b','markeredgecolor','k','markersize',10)

plot([1 1],offsetDiff_12_p([1 3]))
plot([1 1]+1,offsetDiff_21_p([1 3]))
plot([1 1]+2,offsetDiff_22_p([1 3]))
plot([1 ],offsetDiff_12_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,offsetDiff_21_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+2,offsetDiff_22_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)

ylabel('Offset Diff','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 3])
xticklabels({'Offset diff CCW-CW','Offset diff high-low cont','Offset diff HC-LC'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 10)

% 4.B/F Stats

[p_11,h]                    = signrank(1-(10.^slope_cont_ratio_11));
stats.fig4B.stat_11         = p_11;
stats.fig4B.avg_11          = median(10.^slope_cont_ratio_11);

[p_highCont_11]             = signrank(highContSlope_11(:)); % slope diff from zero
[p_lowCont_11]              = signrank(lowContSlope_11(:));
[p_highCont_215]            = signrank(highContQuad_215(:)); % quad diff from zero
[p_lowCont_215]             = signrank(lowContQuad_215(:));
[p_12,h]                    = signrank(offset_choice_diff_12(:));
[p_21,h]                    = signrank(offset_cont_diff_21(:));
[p_22,h]                    = signrank(offset_conf_diff_22(:));
[p_contSlopeDiff,h]         = signrank(highContSlope_11(:),lowContSlope_11(:));
[p_contQuadDiff,h]          = signrank(highContQuad_215(:),lowContQuad_215(:));


stats.fig4B.p_highCont_11   = p_highCont_11;% slope diff from zero
stats.fig4B.p_lowCont_11    = p_lowCont_11;% slope diff from zero
stats.fig4B.stat_12         = p_12; % offset diff from zero
stats.fig4B.avg_12          = median(offset_choice_diff_12(:)); % avg offset diff

stats.fig4D.p_highCont_215  = p_highCont_215;
stats.fig4D.p_lowCont_215   = p_lowCont_215;
stats.fig4D.stat_21         = p_21; % diff offset for contrast
stats.fig4D.avg_21          = median(offset_cont_diff_21(:)); % avg offset diff
stats.fig4D.stat_22         = p_22; % diff offset for conf.
stats.fig4D.avg_22          = median(offset_conf_diff_22(:)); % avg  offset for conf.

stats.fig4B.medHighContslope = median(highContSlope_11(:));
stats.fig4B.medLowContslope  = median(lowContSlope_11(:));

stats.fig4B.p_contSlopeDiff = p_contSlopeDiff;
stats.fig4B.medDiffSlope    = median(highContSlope_11(:)- lowContSlope_11(:));


stats.fig4D.medHighContQuad = median(highContQuad_215(:));
stats.fig4D.medLowContQuad  = median(lowContQuad_215(:));

stats.fig4D.p_contQuadDiff  = p_contQuadDiff;
stats.fig4D.medDiffQuad     = median(highContQuad_215(:)- lowContQuad_215(:));

%% 4.B Choice probabiliy (cross population)

% Offset diff
cp_p         = prctile( NL.CP_choice,[25 50 75]);
ccp_p         = prctile( NL.CP_conf,[25 50 75]);

figure;hold on; axis square
xlim([0 4])
ylim([.4 .7])
plot([0 3],[.5 .5])
plot(ones(length(NL.CP_choice(obsInd)),1)+(rand(length(NL.CP_choice(obsInd)),1)-.5)/5,  NL.CP_choice(obsInd),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(ones(length(NL.CP_choice(~obsInd)),1)+(rand(length(NL.CP_choice(~obsInd)),1)-.5)/5,  NL.CP_choice(~obsInd),'o','markerfacecolor','b','markeredgecolor','k','markersize',10)

plot(1+ones(length(NL.CP_choice(obsInd)),1)+(rand(length(NL.CP_choice(obsInd)),1)-.5)/5,  NL.CP_conf(obsInd),'o','markerfacecolor','y','markeredgecolor','k','markersize',10)
plot(1+ones(length(NL.CP_choice(~obsInd)),1)+(rand(length(NL.CP_choice(~obsInd)),1)-.5)/5,  NL.CP_conf(~obsInd),'o','markerfacecolor','b','markeredgecolor','k','markersize',10)

plot([1 1],cp_p([1 3]))
plot([1 1]+1,ccp_p([1 3]))

plot([1 ],cp_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)
plot([1 ]+1,ccp_p([2]),'d','markerfacecolor','r','markeredgecolor','k','markersize',15)


ylabel('CP','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)
xticks([1 2 ])
xticklabels({'Choice','Confidence'}); set(gca,'TickLabelInterpreter', 'tex','fontName','Helvetica', 'fontAngle', 'oblique', 'fontSize', 16)

%% 4B/D Stats:
stats.fig4B.mean_obs1_CP        = nanmean(NL.CP_choice(observerInd{1}));
stats.fig4B.mean_obs2_CP        = nanmean(NL.CP_choice(observerInd{2}));
stats.fig4D.mean_obs1_CCP       = nanmean(NL.CP_conf(observerInd{1}));
stats.fig4D.mean_obs2_CCP       = nanmean(NL.CP_conf(observerInd{2})); % do we want mean or median...

[p,h]                           = signrank(NL.CP_choice(observerInd{1})-.5 );
stats.fig4B.pval_obs1_CP        = p;

[p,h]                           = signrank(NL.CP_choice(observerInd{2})-.5 );
stats.fig4B.pval_obs2_CP        = p;

[p,h]                           = signrank(NL.CP_conf(observerInd{1})-.5 );
stats.fig4D.pval_obs1_CCP       = p;

[p,h]                           = signrank(NL.CP_conf(observerInd{2})-.5 );
stats.fig4D.pval_obs2_CCP       = p;


% Combined CP:
[p,h]                           = signrank(NL.CP_choice-.5 );
stats.fig4B.pval_CP             = p;
stats.fig4B.avg_CP              = median(NL.CP_choice);
[p,h]                           = signrank(NL.CP_conf-.5 );
stats.fig4D.pval_CCP            = p;
stats.fig4D.avg_CCP             = median(NL.CP_conf);

%% 4.E Choice DV vs Conf DV for example date
choiceLabel                     = {'LC-CCW','HC-CCW','LC-CW','HC-CW'};  % in this order due to the scatterhistogramfunction

% Behavioral data
HC_CCW                          = NL.popStruct{exampleDate}.obsChoice==0 & NL.popStruct{exampleDate}.obsConf==1;
HC_CW                           = NL.popStruct{exampleDate}.obsChoice==1 & NL.popStruct{exampleDate}.obsConf==1;
LC_CCW                          = NL.popStruct{exampleDate}.obsChoice==0 & NL.popStruct{exampleDate}.obsConf==0 ;
LC_CW                           = NL.popStruct{exampleDate}.obsChoice==1 & NL.popStruct{exampleDate}.obsConf==0;
obsChoice                       = NL.popStruct{exampleDate}.obsChoice;
obsConf                         = NL.popStruct{exampleDate}.obsConf;

[idx,CWgroup,HCgroup]           = findgroups(obsChoice,obsConf);
[si,idx2]                       = sort(idx);
ConfChoice                      = strcat(choiceLabel(si));

set(figure, 'OuterPosition', [1800, 1000, 500, 500]); axis square
s                               = scatterhistogram( NL.popStruct{exampleDate}.dvChoice(idx2),  NL.popStruct{exampleDate}.dvConf(idx2),'GroupData',ConfChoice,'LegendVisible','on','HistogramDisplayStyle','smooth', 'LineStyle','-','linewidth',3);
s.Color                         = {'Green','Red','Yellow','Blue',};

%% 4.E Median Split choice DV effect on conf DV for example date
avgConfDV_greater   = nan(11,2,29);
avgConfDV_less      = nan(11,2,29);

stdConfDV_greater   = nan(11,2,29);
stdConfDV_less      = nan(11,2,29);

for iDate = 1:29
    for iC = 1:2
        for iS = 1:length(NL.oriVal{iC,iDate})
            condInd                            = [];
            condInd                            = NL.popStruct{iDate}.contrast== iC & NL.popStruct{iDate}.ori == NL.oriVal{iC,iDate}(iS); % index for a unique contrast-ori val
            
            median_choiceDV                     = median(abs(NL.popStruct{iDate}.dvChoice(condInd)));
            
            greaterInd                          = abs(NL.popStruct{iDate}.dvChoice)>median_choiceDV & condInd;
            lessInd                             = abs(NL.popStruct{iDate}.dvChoice)<=median_choiceDV & condInd;
            
            avgConfDV_greater(iS,iC,iDate)      = median(NL.popStruct{iDate}.dvConf(greaterInd));
            avgConfDV_less(iS,iC,iDate)         = median(NL.popStruct{iDate}.dvConf(lessInd));
            
            stdConfDV_greater(iS,iC,iDate)      = std(NL.popStruct{iDate}.dvConf(greaterInd));
            stdConfDV_less(iS,iC,iDate)         = std(NL.popStruct{iDate}.dvConf(lessInd));
        end
    end
end

%% 4.E - example population median choice Split effect conf confidence DV
figure; axis square; box off; hold on
xlim([-20 20])
ylim([-3 4])
plot(NL.oriVal{1,exampleDate}, avgConfDV_greater(:,1,exampleDate),'ko','color','w','markerfacecolor','g','linewidth',1,'markersize',10)
plot(NL.oriVal{1,exampleDate}, avgConfDV_less(:,1,exampleDate),'ko','color','w','markerfacecolor','r','linewidth',1,'markersize',10)
plot(NL.oriVal{1,exampleDate}, avgConfDV_greater(:,1,exampleDate),'g-')
plot(NL.oriVal{1,exampleDate}, avgConfDV_less(:,1,exampleDate),'r-')
xlabel('Orientation')
ylabel('Avg Conf DV (median split by choice DV)')
title('Example date; median split by choice DV')

%%  4.F per ori stat plot per observer
for iO=1:length(observer)

    if iO==1
        obsInd_temp = ~obsInd; % Ziggy
    elseif iO==2
        obsInd_temp = obsInd; % Friedrich
    end
    
    figure; axis square; box off; hold on
    xlim([-2 3])
    ylim([-2 3])
    plot(squeeze( avgConfDV_less(:,1,obsInd_temp)),squeeze( avgConfDV_greater(:,1,obsInd_temp)),'ko','color','w','markerfacecolor','k','linewidth',1,'markersize',10);
    plot(squeeze( avgConfDV_less(:,2,obsInd_temp)),squeeze( avgConfDV_greater(:,2,obsInd_temp)),'ko','color','w','markerfacecolor',[.5 .5 .5],'linewidth',1,'markersize',10);
    plot([-2 3],[-2 3])
    ylabel('Avg Conf DV (greater than median choice DV)')
    xlabel('Avg Conf DV (less than median choice DV)')
    title(['Observer ' observer{iO} '; Difference (Greater than median - less than median choice DV'])
    
    
    comb_greater            = [squeeze( avgConfDV_greater(:,1,obsInd_temp)) squeeze( avgConfDV_greater(:,2,obsInd_temp))];
    comb_less               = [squeeze( avgConfDV_less(:,1,obsInd_temp)) squeeze( avgConfDV_less(:,2,obsInd_temp))];
    
    figure; hold on; axis square; % histogram
    xlim([-2 2])
    ymax=150;
    ylim([0 ymax])
    edges                   = linspace(min(xlim),max(xlim),15);
    binCenter               = edges + (edges(2)-edges(1))/2;
    dv_diff                 = hist(comb_less(:)-comb_greater(:),binCenter);
    bar(binCenter,dv_diff ,'r')
    plot([0 0],[0 ymax])
    plot(repmat(nanmedian(comb_less(:)-comb_greater(:)),2,1),[0 ymax],'b')
    title(['Choice difference ' observer{iO}])
    
    % Stats 4F
    [p,h]                           = signrank(comb_less(:)-comb_greater(:));
    stats.fig4F.medianSplit_p(iO)   = p;
    stats.fig4F.avgSplit(iO)        = nanmedian(comb_less(:)-comb_greater(:));
end



%% Functions:
function [MSE]= getLinearDiffSlope(params, x,y,nObs)
for iC=1:size(y,2)
    xval        = x{iC}(isfinite(y{iC}));
    yval        = y{iC}(isfinite(y{iC}));
    nObsTemp    = nObs{iC}(isfinite(yval));
    predResp    = params(iC)*(xval) + params(3);
    
    MSE(iC)     = mean((predResp- yval' ).^2 .* nObsTemp');
end
MSE=sum(MSE);
end


function [MSE]= getLinearSameSlope(params, x,y,nObs)

for iC=1:size(y,2)
    xval        = x{iC}(isfinite(y{iC}));
    yval        = y{iC}(isfinite(y{iC}));
    nObsTemp    = nObs{iC}(isfinite(yval));
    predResp    = params(1)*(xval) + params(2);
    
    MSE(iC)     = mean((predResp- yval' ).^2 .* nObsTemp');
end
MSE=sum(MSE);
end


function [MSE]= getLinearDiffOffset(params, x,y,nObs)
for iC=1:2 % loop through choice options (CW/CCW)
    xval        = x;
    yval        = y{iC};
    
    xval        = xval(isfinite(yval));
    yval        = yval(isfinite(yval));
    nObsTemp    = nObs(isfinite(yval),iC);
    predResp    = params(1)*xval+ params(1+iC);
    MSE(iC)     = mean((predResp- yval').^2.* nObsTemp);
end
MSE             = sum(MSE);
end

function [MSE]= getLinearSameOffset(params, x,y,nObs)
for iC=1:2 % loop through choice options (CW/CCW)
    xval        = x;
    yval        = y{iC};
    
    xval        = xval(isfinite(yval));
    yval        = yval(isfinite(yval));
    nObsTemp    = nObs(isfinite(yval),iC);
    predResp    = params(1)*xval+ params(2);
    MSE(iC)     = mean((predResp- yval').^2.* nObsTemp);
end
MSE             = sum(MSE);
end

function [MSE]= getQuadJointDiffOff(params, x,y,nObs) % offset is contrast or confidence dependant
for iC=1:size(y,2)
    xval        = x{iC};
    yval        = y{iC};
    
    xval        = xval(isfinite(yval));
    yval        = yval(isfinite(yval));
    nObsTemp    = nObs{iC}(isfinite(yval));
    predResp    = params(1)*(xval-params(2)).^2 + params(2+iC);
    
    MSE(iC)     = mean((predResp- yval' ).^2 .* nObsTemp');
end
MSE             = sum(MSE);
end

function [MSE]= getQuadJointSameOff(params, x,y,nObs) % offset is contrast or confidence dependant
for iC=1:size(y,2)
    xval        = x{iC};
    yval        = y{iC};
    
    xval        = xval(isfinite(yval));
    yval        = yval(isfinite(yval));
    nObsTemp    = nObs{iC}(isfinite(yval));
    predResp    = params(1)*(xval-params(2)).^2 + params(3);
    
    MSE(iC)     = mean((predResp- yval' ).^2 .* nObsTemp');
end
MSE             = sum(MSE);
end
