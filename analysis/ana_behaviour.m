%% ana_behaviour.m
% Behavioural analysis script
%% TO DO:
% Major overhaul!!!
clear all
% close all
clc

%% Load participant data

a{1}='01_01';
% a{2}='02_01';
a{2}='02_02';
a{3}='03_01';
a{4}='04_01';
a{5}='05_01';
a{6}='06_01';
a{7}='07_01';
a{8}='08_01';
a{9}='09_01';
a{10}='10_01';
a{11}='11_01';
a{12}='12_01';
a{13}='13_01';
a{14}='14_01';
a{15}='15_01';
a{16}='16_01';
a{17}='17_01';
a{18}='18_01';
a{19}='19_01';
a{20}='20_01';
a{21}='21_01';
a{22}='22_01';
a{23}='23_01';

%% Settings

    % Specify participants of analysis
    Ptt = [7,9,19,20];
    Ptt = [1:6,8,10:18,20:22]; % bad ptt removed (7(bad and blinky),9(bad and med-blinky),19(bad and blinky), 20(maybe - ok but blinky), 23(maybe - bad but non-blinky))
    
    % Error Trial removal (1 = remove, 0 = leave in)
    errTrialRmvl = 1;

    % First Trials removal
    firstTrialRmvl = 1;

    % Valid location expectation only
    invalLoc = 1;

    % Filter RT outliers
    % 0 = none; 1 = MAD(median absolute deviation); 2 = 2.5% percentile; 3 = exclusion of absolute RT (i.e. longer than 1 second)
    outlierScheme = 1;

    % Individual plots
    indivPlots = 0;
    
    % Measures of interest
    
    plotRTMed = 1;
    plotRTMedDiff = 1;

    % Give T-test output:
    ttestOutput = 1;     

%% Individual data set analyses

PttCount=1;

for PttIterat = 1:length(Ptt);
    fid = fopen(a{Ptt(PttCount)});
    file = strcat('/Users/Yannik/Documents/MATLAB/Distracter_Exptectation/EEG_results/','s', a{Ptt(PttCount)},'/','s', a{Ptt(PttCount)},'Behav.mat');
    Bhv = []; % create Bhv analysis matrix
    Bhv = [Bhv;load(file)]; % loads Bhv output data file for Bhv matrix
    outputDir = ['/Users/Yannik/Documents/MATLAB/Distracter_Exptectation/EEG_results/'];
    file = fullfile(outputDir, 'BehavAnalysis.m');
    
    %% Identify error trials
    
    % clear previous error trials & outliers (redundant if clear all)
    clear Bhv.outputDat.errTrials
    
    if errTrialRmvl == 1;
        errorTrialCounter=1; % absolute errorTrial counter
        for errorTrialIteration=1:length(Bhv.outputDat.data(:,1)) % look at outputDat in Bhv
            if  Bhv.outputDat.data(errorTrialIteration,14) == 0 % column 14 specifies response accuracy (0=incorrect)
                Bhv.outputDat.errTrials(errorTrialCounter,:) = Bhv.outputDat.data(errorTrialIteration,:); % puts in the whole data matrix of the error trial
                errorTrialCounter=errorTrialCounter+1;
            end
        end
    end
    
    
    %% Identify invalid trials
    % affects only 75% validity scheme
    
    invalLocLoop=1;
    if invalLoc == 1;
        Bhv.outputDat.data(:,28) = zeros; % prepare column 28 to contain invalid location trials (makes for-loop faster)
        for invalLocIteration=1:length(Bhv.outputDat.data(:,1))
            if Bhv.outputDat.data(invalLocIteration,3) == 2; % if it is a 75% T or D trial
                if Bhv.outputDat.data(invalLocIteration,5) == 1; % If in T manipulation
                    if Bhv.outputDat.data(invalLocIteration,6) ~= Bhv.outputDat.data(invalLocIteration,7); % if expected location is NOT T loc
                        Bhv.outputDat.data(invalLocIteration,28) = 1; % column for invalid locations
                        Bhv.outputDat.invalLocs(invalLocLoop,:) = Bhv.outputDat.data(invalLocIteration,:); % make matrix containing all info on invalid trials
                        invalLocLoop = invalLocLoop+1;
                    end
                elseif Bhv.outputDat.data(invalLocIteration,5) == 2; % If in D manipulation
                    if Bhv.outputDat.data(invalLocIteration,6) ~= Bhv.outputDat.data(invalLocIteration,8) % if expected location is NOT D loc
                        Bhv.outputDat.data(invalLocIteration,28) = 1; % column for invalid locations
                        Bhv.outputDat.invalLocs(invalLocLoop,:) = Bhv.outputDat.data(invalLocIteration,:); % make matrix containing all info on invalid trials
                        invalLocLoop = invalLocLoop+1;
                    end
                end
%                 invalLocLoop = invalLocLoop+1;
            end
        end
    end
    
    
    %% Identify first trials of block
    % Time needed by Pt to learn statistical regularities
    
    if firstTrialRmvl == 1;
        firstTrialRmvlCount = 1;
        for firstTrialRmvlIteration = 1:length(Bhv.outputDat.data(:,1))
            if Bhv.outputDat.data(firstTrialRmvlIteration,2) < 6 % if it is within the first 5 trials of each block
                Bhv.outputDat.data(firstTrialRmvlIteration,29) = 1; % column for first trials
                Bhv.outputDat.firstTrialRmvl(firstTrialRmvlCount,:) = Bhv.outputDat.data(firstTrialRmvlIteration,:); % make matrix containing all info on first trials
                firstTrialRmvlCount = firstTrialRmvlCount + 1;
            end
        end
    end
    
    %% Identify outliers
    % NB: Outliers defined within subject across all conditions
    outlierCount = 1; % set outlier counter to 0
    
    % No outlier scheme
    if outlierScheme == 0
        Bhv.outputDat.data(:,30) = zeros; % column for outliers
        
        % Median Absolute Deviation (MAD) (from the median)
    elseif outlierScheme == 1
        medianRT = median(Bhv.outputDat.data(:,22)); % column 22 = RTs
        madRT = mad(Bhv.outputDat.data(:,22),1);
        for i=1:length(Bhv.outputDat.data(:,1))
            if Bhv.outputDat.data(i,22) > medianRT+4*madRT || Bhv.outputDat.data(i,22) < medianRT-4*madRT
                Bhv.outputDat.data(i,30) = 1; % column for outliers
                Bhv.outputDat.outlier(outlierCount,:) = Bhv.outputDat.data(i,:); % make matrix containing all info on outlier trials
                outlierCount = outlierCount+1;
            end
        end
        
        % 2.5th %-ile
    elseif outlierScheme == 2
        highestprctle = prctile(Bhv.outputDat.data(:,22), 97.5);
        lowestprctle = prctile(Bhv.outputDat.data(:,22), 2.5);
        for i=1:length(Bhv.outputDat.data(:,1))
            if Bhv.outputDat.data(i,22)> highestprctle || Bhv.outputDat.data(i,22)< lowestprctle
                Bhv.outputDat.data(i,30) = 1; % column for outliers
                Bhv.outputDat.outlier(outlierCount,:) = Bhv.outputDat.data(i,:); % make matrix containing all info on outlier trials
                %
                outlierCount=outlierCount+1;
            end
        end
        
        % Absolute RT lower and upper limits
    elseif outlierScheme == 3
        for i=1:length(Bhv.outputDat.data(:,1))
            if Bhv.outputDat.data(i,22) < 0.25 || Bhv.outputDat.data(i,22) > 1 %Only excluding trials shorter than .25 and those longer than 1 second.
                Bhv.outputDat.data(i,30) = 1; % column for outliers
                Bhv.outputDat.outlier(outlierCount,:) = Bhv.outputDat.data(i,:); % make matrix containing all info on outlier trials
                outlierCount=outlierCount+1;
            end
        end
    end
    
    %% Identify different conditions and valid trials
    
    % Separate counters per condition
    % (required to make trial matrices per condition; reset to 1 for every new ptt)
    T25Count=1; T75Count=1; T100Count=1; D25Count=1; D75Count=1; D100Count=1;
    % Error counters
    % (required because number of error trials not in valid trial per condition matrix)
    T25ErrCount=1; T75ErrCount=1; T100ErrCount=1; D25ErrCount=1; D75ErrCount=1; D100ErrCount=1;
    
    T75InvalLocCount=1; D75InvalLocCount=1;
    
    clear T25 T75 T100 D25 D75 D100 % v. important!
    
    for i=1:length(Bhv.outputDat.data(:,1))
        if Bhv.outputDat.data(i,5)==1 % T
            if Bhv.outputDat.data(i,3)==1 % T25
                if errTrialRmvl == 1; % If we want to remove errTrials
                    % Include trial only if: correct (i,14), valid location (i,28), not at the start of the block (i,29), no outlier (i,30)
                    if Bhv.outputDat.data(i,14)==1 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T25(T25Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        T25Tot(T25Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        T25Count = T25Count+1;
                        % elseif incorrect trial + other criteria
                    elseif Bhv.outputDat.data(i,14)==0 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T25ErrCount=T25ErrCount+1;
                    end
                elseif errTrialRmvl==0; % if we want to include errTrials (NB: this will not affect accuracy measure but RT measure as the latter is taken from corr and incorr trials)
                    if Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T25(T25Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        T25Tot(T25Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        T25Count = T25Count+1;
                        if Bhv.outputDat.data(i,14)==0; % count all error trials
                            T25ErrCount=T25ErrCount+1;
                        end
                    end
                end
            elseif Bhv.outputDat.data(i,3)==2 % T75
                if errTrialRmvl == 1; % If we want to remove errTrials
                    % Include trial only if: correct (i,14), valid location (i,28), not at the start of the block (i,29), no outlier (i,30)
                    if Bhv.outputDat.data(i,14)==1 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T75(T75Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        T75Tot(T75Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        T75Count = T75Count+1;
                        % elseif incorrect trial + other criteria
                    elseif Bhv.outputDat.data(i,14)==0 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T75ErrCount=T75ErrCount+1;
                    elseif Bhv.outputDat.data(i,28)==1 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0 % takes in invalid T75 trials that are not outliers, not within first block trials, correct or incorrect
                        T75InvalLoc(T75InvalLocCount,:) = Bhv.outputDat.data(i,:);
                        T75InvalLocCount=T75InvalLocCount+1;
                    end
                elseif errTrialRmvl==0; % if we want to include errTrials
                    if Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T75(T75Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        T75Tot(T75Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        T75Count = T75Count+1;
                        if Bhv.outputDat.data(i,14)==0; % count all error trials
                            T75ErrCount=T75ErrCount+1;
                        end
                    end
                end
            elseif Bhv.outputDat.data(i,3)==3 % T100
                if errTrialRmvl == 1; % If we want to remove errTrials
                    % Include trial only if: correct (i,14), valid location (i,28), not at the start of the block (i,29), no outlier (i,30)
                    if Bhv.outputDat.data(i,14)==1 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T100(T100Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        T100Tot(T100Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        T100Count = T100Count+1;
                        % elseif incorrect trial + other criteria
                    elseif Bhv.outputDat.data(i,14)==0 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T100ErrCount=T100ErrCount+1;
                    end
                elseif errTrialRmvl==0; % if we want to include errTrials
                    if Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        T100(T100Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        T100Tot(T100Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        T100Count = T100Count+1;
                        if Bhv.outputDat.data(i,14)==0; % count all error trials
                            T100ErrCount=T100ErrCount+1;
                        end
                    end
                end
            end
        elseif Bhv.outputDat.data(i,5)==2 % D
            if Bhv.outputDat.data(i,3)==1 % D25
                if errTrialRmvl == 1; % If we want to remove errTrials
                    % Include trial only if: correct (i,14), valid location (i,28), not at the start of the block (i,29), no outlier (i,30)
                    if Bhv.outputDat.data(i,14)==1 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D25(D25Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        D25Tot(D25Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        D25Count = D25Count+1;
                        % elseif incorrect trial + other criteria
                    elseif Bhv.outputDat.data(i,14)==0 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D25ErrCount=D25ErrCount+1;
                    end
                elseif errTrialRmvl==0; % if we want to include errTrials
                    if Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D25(D25Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        D25Tot(D25Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        D25Count = D25Count+1;
                        if Bhv.outputDat.data(i,14)==0; % count all error trials
                            D25ErrCount=D25ErrCount+1;
                        end
                    end
                end
            elseif Bhv.outputDat.data(i,3)==2 % D75
                if errTrialRmvl == 1; % If we want to remove errTrials
                    % Include trial only if: correct (i,14), valid location (i,28), not at the start of the block (i,29), no outlier (i,30)
                    if Bhv.outputDat.data(i,14)==1 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D75(D75Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        D75Tot(D75Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        D75Count = D75Count+1;
                        % elseif incorrect trial + other criteria
                    elseif Bhv.outputDat.data(i,14)==0 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D75ErrCount=D75ErrCount+1;
                    elseif Bhv.outputDat.data(i,28)==1 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0 % takes in invalid D75 trials that are not outliers, not within first block trials, correct or incorrect
                        D75InvalLoc(D75InvalLocCount,:) = Bhv.outputDat.data(i,:);
                        D75InvalLocCount=D75InvalLocCount+1;
                    end
                elseif errTrialRmvl==0; % if we want to include errTrials
                    if Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D75(D75Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        D75Tot(D75Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        D75Count = D75Count+1;
                        if Bhv.outputDat.data(i,14)==0; % count all error trials
                            D75ErrCount=D75ErrCount+1;
                        end
                    end
                end
            elseif Bhv.outputDat.data(i,3)==3 % D100
                if errTrialRmvl == 1; % If we want to remove errTrials
                    % Include trial only if: correct (i,14), valid location (i,28), not at the start of the block (i,29), no outlier (i,30)
                    if Bhv.outputDat.data(i,14)==1 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D100(D100Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        D100Tot(D100Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        D100Count = D100Count+1;
                        % elseif incorrect trial + other criteria
                    elseif Bhv.outputDat.data(i,14)==0 && Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D100ErrCount=D100ErrCount+1;
                    end
                elseif errTrialRmvl==0; % if we want to include errTrials
                    if Bhv.outputDat.data(i,28)==0 && Bhv.outputDat.data(i,29)==0 && Bhv.outputDat.data(i,30)==0
                        D100(D100Count,:) = Bhv.outputDat.data(i,:); % trials per condition per ptt ~~~ still nec if we have Total counter
                        D100Tot(D100Count,:, PttCount) = Bhv.outputDat.data(i,:); % Total trials across ptts
                        D100Count = D100Count+1;
                        if Bhv.outputDat.data(i,14)==0; % count all error trials
                            D100ErrCount=D100ErrCount+1;
                        end
                    end
                end
            end
        end
    end % end trial-per-condition loop
    
    % Subtract 1l from each condition counter to counter-act the very last addition at the end of this for-loop
    T25Count=T25Count-1; T75Count=T75Count-1; T100Count=T100Count-1; D25Count=D25Count-1; D75Count=D75Count-1; D100Count=D100Count-1;
    T25ErrCount=T25ErrCount-1; T75ErrCount=T75ErrCount-1; T100ErrCount=T100ErrCount-1; D25ErrCount=D25ErrCount-1; D75ErrCount=D75ErrCount-1; D100ErrCount=D100ErrCount-1;
    
    
    %% Behavioural Analyses
    
    %% Individual Measures
    
    % Median RTs per ptt (row), per condition (column) % ~~~3rd dim nec?
    RTMed(PttCount,1,1)=median(T25(:,22));
    RTMed(PttCount,2,1)=median(T75(:,22));
    RTMed(PttCount,3,1)=median(T100(:,22));
    RTMed(PttCount,4,1)=median(D25(:,22));
    RTMed(PttCount,5,1)=median(D75(:,22));
    RTMed(PttCount,6,1)=median(D100(:,22));
    
    % RT Med inval 75% locs
    RTMed75InvalLoc(PttCount,3,1)=median(T75InvalLoc(:,22)); % inval
    RTMed75InvalLoc(PttCount,4,1)=median(D75InvalLoc(:,22));
    RTMed75InvalLoc(PttCount,1,1)=RTMed(PttCount,2,1); % val
    RTMed75InvalLoc(PttCount,2,1)=RTMed(PttCount,5,1);
    
    % Use Grand avg of Mean for each ptt and get difference to each
    % condition (to reduce SE)
%     SubjectMean(PttCount,1) = mean(RTMed(PttCount,:));
    
%     RTMedDiff(PttCount,:) = RTMed(PttCount,:) - SubjectMean(PttCount,1);
        
    % Mean RTs per ptt (row), per condition (column)
    RTM(PttCount,1,1)=mean(T25(:,22));
    RTM(PttCount,2,1)=mean(T75(:,22));
    RTM(PttCount,3,1)=mean(T100(:,22));
    RTM(PttCount,4,1)=mean(D25(:,22));
    RTM(PttCount,5,1)=mean(D75(:,22));
    RTM(PttCount,6,1)=mean(D100(:,22));
       
    %Standard Error ~~~ should use MAD or quartile for median
    SE(PttCount,1,1)=std(T25(:,22))/sqrt(length(T25(:,22))); %T25
    SE(PttCount,2,1)=std(T75(:,22))/sqrt(length(T75(:,1))); %T75
    SE(PttCount,3,1)=std(T100(:,22))/sqrt(length(T100(:,1))); %T100
    SE(PttCount,4,1)=std(D25(:,22))/sqrt(length(D25(:,1))); %D25
    SE(PttCount,5,1)=std(D75(:,22))/sqrt(length(D75(:,1))); %D75
    SE(PttCount,6,1)=std(D100(:,22))/sqrt(length(D100(:,1))); %D100
        
    % SE 75% inval Locs
    SE75InvalLoc(PttCount,1,1)=std(T75InvalLoc(:,22))/sqrt(length(T75InvalLoc(:,1)));
    SE75InvalLoc(PttCount,2,1)=std(D75InvalLoc(:,22))/sqrt(length(D75InvalLoc(:,1)));
    
    % MAD (larger than SE by factor 10!)
    RTmad(PttCount,1,1) = mad(T25(:,22),1);
    RTmad(PttCount,2,1) = mad(T25(:,22),1);
    RTmad(PttCount,3,1) = mad(T25(:,22),1);
    RTmad(PttCount,4,1) = mad(T25(:,22),1);
    RTmad(PttCount,5,1) = mad(T25(:,22),1);
    RTmad(PttCount,6,1) = mad(T25(:,22),1);
        
    % Accuracy
    Acc(PttCount,1,1)=(T25Count/(T25Count+T25ErrCount));
    Acc(PttCount,4,1)=(D25Count/(D25Count+D25ErrCount));
    Acc(PttCount,2,1)=(T75Count/(T75Count+T75ErrCount)); % NB: count in inval 75T loc to reflect true ratio
    Acc(PttCount,5,1)=(D75Count/(D75Count+D75ErrCount)); % NB: count in inval 75D loc to reflect true ratio
    Acc(PttCount,3,1)=(T100Count/(T100Count+T100ErrCount));
    Acc(PttCount,6,1)=(D100Count/(D100Count+D100ErrCount));
    
    Acc(PttCount,:,1) = Acc(PttCount,:,1)*100;% round(Acc(PttCount,:,1)*100); % Round accuracies to % integers
        
    Acc75InvalLoc(PttCount,3,1)=sum(T75InvalLoc(:,14))/length(T75InvalLoc(:,14)); % inval
    Acc75InvalLoc(PttCount,4,1)=sum(D75InvalLoc(:,14))/length(D75InvalLoc(:,14)); 
    Acc75InvalLoc(PttCount,1,1)=(T75Count/(T75Count+T75ErrCount)); % val
    Acc75InvalLoc(PttCount,2,1)=(T75Count/(T75Count+T75ErrCount));
    
    Acc75InvalLoc(PttCount,:,1) = Acc75InvalLoc(PttCount,:,1)*100;% round(Acc(PttCount,:,1)*100); % Round accuracies to % integers
    
    % Number of trials per condition per subject
    N(PttCount,1,1)=T25Count;
    N(PttCount,4,1)=D25Count;
    N(PttCount,2,1)=T75Count;
    N(PttCount,5,1)=D75Count;
    N(PttCount,3,1)=T100Count;
    N(PttCount,6,1)=D100Count;
        
    
    %% Individual RT Plots
    
    if indivPlots==1
        
        Fig = figure; % Open figure window
        
        % Make three subplots into Fig (RT, Acc, N)
        % RT
        subplot(3,1,1);
        plotMedRTT = plot(RTMed(PttCount,1:3,1)); % Plot RTs for T condition
        title(a{Ptt(PttCount)}); xlabel ('Condition'); ylabel('Median RT');
        set(plotMedRTT,'Color','red','LineWidth',2);
        set(gca, 'XTick',1:3, 'XTickLabel',{ '25' '75' '100' });
        hold on % hold on to current sbplot to put rest of info it
        plotRTmD = plot(RTMed(PttCount,4:6,1)); % RTs for Ds
        set(plotRTmD,'Color','blue','LineWidth',2);
        plotSE = errorbar(RTMed(PttCount,1:3,1), SE(PttCount,1:3,1), 'k', 'linestyle','none'); % SEs for T
        plotSE = errorbar(RTMed(PttCount,4:6,1), SE(PttCount,4:6,1), 'k', 'linestyle','none'); % SEs for D
        legend('T','D')
        hold off
        
        % Acc
        subplot(3,1,2);
        plotAcc = bar([Acc(PttCount,1:3,1)' Acc(PttCount,4:6,1)']);
        title('Accuracy'); xlabel ('Condition'); ylabel('Accuracy (%)'); ylim([75 100])
        set(plotAcc(1),'facecolor','red');
        set(plotAcc(2),'facecolor','blue');
        set(gca, 'XTick',1:3, 'XTickLabel',{'25' '75' '100' });
        
        % N
        subplot(3,1,3);
        plotN = bar([N(PttCount,1:3,1)' N(PttCount,4:6,1)']);
        title('Number of Trials'); xlabel ('Condition'); ylabel('Number of trials');
        set(plotN(1),'facecolor','red');
        set(plotN(2),'facecolor','blue');
        set(gca, 'XTick',1:3, 'XTickLabel',{ '25' '75' '100' });
        
        
    end
    
    %% Individual Scatterplots
    % to analyse RT decrease throughout block
    
    %     % ???
    %     for j=1:max(T25(:,2))
    %         data_trialn=find(data_T1(:,2) == j);
    %         mTrialnT1(j,PttCount)= median(data_T1(data_trialn,22));
    %     end
    %
    %     for j = 1:max(data_T2(:,2))
    %         data_trialn=find(data_T2(:,2) == j);
    %         mTrialnT2(j,PttCount)= median(data_T2(data_trialn,22));
    %     end
    %     for j = 1:max(data_T3(:,2))
    %         data_trialn=find(data_T3(:,2) == j);
    %         mTrialnT3(j,PttCount)= median(data_T3(data_trialn,22));
    %     end
    %     for j = 1:max(data_D1(:,2))
    %         data_trialn=find(data_D1(:,2) == j);
    %         mTrialnD1(j,PttCount)= median(data_D1(data_trialn,22));
    %     end
    %     for j = 1:max(data_D2(:,2))
    %         data_trialn=find(data_D2(:,2) == j);
    %         mTrialnD2(j,PttCount)= median(data_D2(data_trialn,22));
    %     end
    %     for j = 1:max(data_D3(:,2))
    %         data_trialn=find(data_D3(:,2) == j);
    %         mTrialnD3(j,PttCount)= median(data_D3(data_trialn,22));
    %     end
    %
    
    PttCount=PttCount+1;
    
end % end Ptt-loop

% For Total analyses, subtract 1 from the PttCount to reflect true number of ptts
PttCount=PttCount-1;

%% TOTAL group measures

Bhv.outputDat.RTMed = RTMed; % Put info into outputDat
Bhv.outputDat.RTmad = RTmad; % Put info into outputDat
Bhv.outputDat.RTM = RTM; % Put info into outputDat
Bhv.outputDat.SE = SE; % Put info into outputDat
Bhv.outputDat.Acc = Acc; % Put info into outputDat
Bhv.outputDat.N = N; % Put info into outputDat

% Total RT Mean of indiv Medians
% NB: 3rd dimension redundant, left over from previous script
TotRTm(1,1,1)=mean(RTMed(:,1));
TotRTm(1,2,1)=mean(RTMed(:,2));
TotRTm(1,3,1)=mean(RTMed(:,3));
TotRTm(1,4,1)=mean(RTMed(:,4));
TotRTm(1,5,1)=mean(RTMed(:,5));
TotRTm(1,6,1)=mean(RTMed(:,6));

% same for 75%InvalLocs
TotRTm75InvalLoc(1,2)=mean(RTMed75InvalLoc(:,3)); % inval
TotRTm75InvalLoc(2,2)=mean(RTMed75InvalLoc(:,4));
TotRTm75InvalLoc(1,1)=TotRTm(1,2,1); % val
TotRTm75InvalLoc(2,1)=TotRTm(1,5,1);


% Total RT Mean of Med indiv Differences
TotRTmDiff(1,1,1)=mean(RTMedDiff(:,1));
TotRTmDiff(1,2,1)=mean(RTMedDiff(:,2));
TotRTmDiff(1,3,1)=mean(RTMedDiff(:,3));
TotRTmDiff(1,4,1)=mean(RTMedDiff(:,4));
TotRTmDiff(1,5,1)=mean(RTMedDiff(:,5));
TotRTmDiff(1,6,1)=mean(RTMedDiff(:,6));

% % Total RTs for all subjects - Mean of indiv Medians 
% diff format for diff graph
% TotRTm(1,1)=mean(RTMed(:,1));
% TotRTm(2,1)=mean(RTMed(:,2));
% TotRTm(3,1)=mean(RTMed(:,3));
% TotRTm(1,2)=mean(RTMed(:,4));
% TotRTm(2,2)=mean(RTMed(:,5));
% TotRTm(3,2)=mean(RTMed(:,6));

Bhv.outputDat.TotRTm = TotRTm; % Put info into outputDat

% Total RT SE
TotSE(1,1)=std(RTMed(:,1))/sqrt(length(Ptt));
TotSE(1,2)=std(RTMed(:,2))/sqrt(length(Ptt));
TotSE(1,3)=std(RTMed(:,3))/sqrt(length(Ptt));
TotSE(1,4)=std(RTMed(:,4))/sqrt(length(Ptt));
TotSE(1,5)=std(RTMed(:,5))/sqrt(length(Ptt));
TotSE(1,6)=std(RTMed(:,6))/sqrt(length(Ptt));

% same for 75% invalLocs
TotSE75InvalLoc(1,2)=std(RTMed75InvalLoc(:,3))/sqrt(length(Ptt)); % inval
TotSE75InvalLoc(2,2)=std(RTMed75InvalLoc(:,4))/sqrt(length(Ptt));
TotSE75InvalLoc(1,1)=TotSE(1,2); % val
TotSE75InvalLoc(2,1)=TotSE(1,5);


% Total RT MADs
TotMAD(1,1) = mad(RTMed(:,1));
TotMAD(1,2) = mad(RTMed(:,2));
TotMAD(1,3) = mad(RTMed(:,3));
TotMAD(1,4) = mad(RTMed(:,4));
TotMAD(1,5) = mad(RTMed(:,5));
TotMAD(1,6) = mad(RTMed(:,6));


% Total RT SE Med Diff
TotSEMedDiff(1,1)=std(RTMedDiff(:,1))/sqrt(length(Ptt));
TotSEMedDiff(1,2)=std(RTMedDiff(:,2))/sqrt(length(Ptt));
TotSEMedDiff(1,3)=std(RTMedDiff(:,3))/sqrt(length(Ptt));
TotSEMedDiff(1,4)=std(RTMedDiff(:,4))/sqrt(length(Ptt));
TotSEMedDiff(1,5)=std(RTMedDiff(:,5))/sqrt(length(Ptt));
TotSEMedDiff(1,6)=std(RTMedDiff(:,6))/sqrt(length(Ptt));

% % Total RT SE
% diff format for diff graph
% TotSE(1,1)=std(RTMed(:,1))/sqrt(length(Ptt));
% TotSE(2,1)=std(RTMed(:,2))/sqrt(length(Ptt));
% TotSE(3,1)=std(RTMed(:,3))/sqrt(length(Ptt));
% TotSE(1,2)=std(RTMed(:,4))/sqrt(length(Ptt));
% TotSE(2,2)=std(RTMed(:,5))/sqrt(length(Ptt));
% TotSE(3,2)=std(RTMed(:,6))/sqrt(length(Ptt));

Bhv.outputDat.TotSE = TotSE; % Put info into outputDat

% Total Accuracy Mean
TotAccM(1,1)=mean(Acc(:,1));
TotAccM(2,1)=mean(Acc(:,2));
TotAccM(3,1)=mean(Acc(:,3));
TotAccM(1,2)=mean(Acc(:,4));
TotAccM(2,2)=mean(Acc(:,5));
TotAccM(3,2)=mean(Acc(:,6));

% same for 75%InvalLocs
TotAccM75InvalLoc(1,2)=mean(Acc75InvalLoc(:,3)); % val
TotAccM75InvalLoc(2,2)=mean(Acc75InvalLoc(:,4));
TotAccM75InvalLoc(1,1)=TotAccM(2,1); % inval
TotAccM75InvalLoc(2,1)=TotAccM(2,2);

Bhv.outputDat.TotAccM = TotAccM; % Put info into outputDat

% Total Accuracy SE
TotAccSE(1,1)=std(Acc(:,1))/sqrt(length(Ptt));
TotAccSE(2,1)=std(Acc(:,2))/sqrt(length(Ptt));
TotAccSE(3,1)=std(Acc(:,3))/sqrt(length(Ptt));
TotAccSE(1,2)=std(Acc(:,4))/sqrt(length(Ptt));
TotAccSE(2,2)=std(Acc(:,5))/sqrt(length(Ptt));
TotAccSE(3,2)=std(Acc(:,6))/sqrt(length(Ptt));

% same for 75% ivnalLocs
TotAccSE75InvalLoc(1,2)=std(Acc75InvalLoc(:,3))/sqrt(length(Ptt)); % inval
TotAccSE75InvalLoc(2,2)=std(Acc75InvalLoc(:,4))/sqrt(length(Ptt));
TotAccSE75InvalLoc(1,1)= TotAccSE(2,1); % val
TotAccSE75InvalLoc(2,1)= TotAccSE(2,2);


Bhv.outputDat.TotAccSE = TotAccSE; % Put info into outputDat

TotN(1,1) = sum(N(:,1));
TotN(1,2) = sum(N(:,2));
TotN(1,3) = sum(N(:,3));
TotN(1,4) = sum(N(:,4));
TotN(1,5) = sum(N(:,5));
TotN(1,6) = sum(N(:,6));

Bhv.outputDat.TotN = TotN; % Put info into outputDat

TotNProp = TotN(1,:)/length(Ptt);


%% within-subject error bar normalization

%subject means
for i=1:length(RTMed(:,1));
    SubjectMean(i,1) = mean(RTMed(i,:));
end

%grand average
GrndAvrg = mean(SubjectMean(:,1));

%normalization
for l=1:length(RTMed(:,1));
    for m=1:length(RTMed(1,:));
        normRT(l,m) = ...
            RTMed(l,m) - SubjectMean(l,1) + GrndAvrg;
    end
end

clear TOTnormRT
for n=1:length(normRT(1,:))
    TOTnormRT(1,n) = mean(normRT(:,1));
end

% SEs
for o=1:length(normRT(1,:))
    normRTSE(1,o) = std(normRT(o,:))/sqrt(length(normRT(:,1)));
end


figure
newplot = plot(TotRTm(1,1:3,1));
title('Group Mean of Median RTs'); xlabel('Spatial Predictability (%)'); ylabel('RT (s)');
set(gca, 'XTick',1:3, 'XTickLabel',{ '25' '75' '100' });% ylim([0.495 0.55])
set(newplot,'Color','red','LineWidth',2);
hold on
newplot2 = plot(TotRTm(1,4:6,1));
set(newplot2,'Color','blue','LineWidth',2); 
plotSE = errorbar(TotRTm(1,1:3,1), normRTSE(1,1:3,1), 'k', 'linestyle','none');
plotSE = errorbar(TotRTm(1,4:6,1), normRTSE(1,4:6,1), 'k', 'linestyle','none');
legend('T','D');
hold off

%% Total group figure

Fig=figure;

% Total Mean of Median RTs
% subplot(2,1,1);
plotRTmT = plot(TotRTm(1,1:3,1));
title('Group Mean of Median RTs'); xlabel('Spatial Predictability (%)'); ylabel('RT (s)');
set(gca, 'XTick',1:3, 'XTickLabel',{ '25' '75' '100' });% ylim([0.495 0.55])
set(plotRTmT,'Color','red','LineWidth',2);
hold on
plotRTmD = plot(TotRTm(1,4:6,1));
set(plotRTmD,'Color','blue','LineWidth',2); 
plotSE = errorbar(TotRTm(1,1:3,1), TotSE(1,1:3,1), 'k', 'linestyle','none');
plotSE = errorbar(TotRTm(1,4:6,1), TotSE(1,4:6,1), 'k', 'linestyle','none');
legend('T','D');
hold off

% Total Accuracy
Fig=figure;
% subplot(2,1,2);
plotTotAccM = barwitherr(TotAccSE(:,:,1), TotAccM(:,:,1)); % maybe leave 1 at 3rd dim
set(gca, 'XTick',1:3, 'XTickLabel',{ '25' '75' '100' }); ylim([89 100])
set(plotTotAccM(1),'facecolor','red'); set(plotTotAccM(2),'facecolor','blue');
title('Group Mean Accuracies'); xlabel('Spatial Predictability (%)'); ylabel('Accuracy (%)');
legend('Target','Distractor');

%% Alpha power analysis
clear meanAlpha
meanAlpha(1,1) = mean(spssDat(:,1));
meanAlpha(2,1) = mean(spssDat(:,2));
meanAlpha(3,1) = mean(spssDat(:,3));
meanAlpha(1,2) = mean(spssDat(:,4));
meanAlpha(2,2) = mean(spssDat(:,5));
meanAlpha(3,2) = mean(spssDat(:,6));

clear SEAlpha
SEAlpha(1,1) = std(spssDat(:,1))/sqrt(length(spssDat(:,1)));
SEAlpha(2,1) = std(spssDat(:,2))/sqrt(length(spssDat(:,1)));
SEAlpha(3,1) = std(spssDat(:,3))/sqrt(length(spssDat(:,1)));
SEAlpha(1,2) = std(spssDat(:,4))/sqrt(length(spssDat(:,1)));
SEAlpha(2,2) = std(spssDat(:,5))/sqrt(length(spssDat(:,1)));
SEAlpha(3,2) = std(spssDat(:,6))/sqrt(length(spssDat(:,1)));

% subject averages
subjMeanAlpha = [];
ct = 1;
i=1;
for i= 1:length(spssDat(:,1));
    subjMeanAlpha(ct,1) = mean(spssDat(i,:));
    ct=ct+1;
end

% grand average
grndAvg = mean(subjMeanAlpha(:,1));

% normalization to grand average

for subjLength=1:length(spssDat(:,1));
    for condLength=1:length(spssDat(1,:));
        normAlpha(subjLength,condLength) = ...
            spssDat(subjLength,condLength) - subjMeanAlpha(subjLength,1) + grndAvg;
    end
end

for i= 1:length(normAlpha(:,1));
    testsubjMeanAlpha(i,1) = mean(normAlpha(i,:));
end

clear NEWmeanAlpha
NEWmeanAlpha(1,1) = mean(normAlpha(:,1));
NEWmeanAlpha(2,1) = mean(normAlpha(:,2));
NEWmeanAlpha(3,1) = mean(normAlpha(:,3));
NEWmeanAlpha(1,2) = mean(normAlpha(:,4));
NEWmeanAlpha(2,2) = mean(normAlpha(:,5));
NEWmeanAlpha(3,2) = mean(normAlpha(:,6));

clear NEWSEAlpha
NEWSEAlpha(1,1) = std(normAlpha(:,1))/sqrt(length(normAlpha(:,1)));
NEWSEAlpha(2,1) = std(normAlpha(:,2))/sqrt(length(normAlpha(:,1)));
NEWSEAlpha(3,1) = std(normAlpha(:,3))/sqrt(length(normAlpha(:,1)));
NEWSEAlpha(1,2) = std(normAlpha(:,4))/sqrt(length(normAlpha(:,1)));
NEWSEAlpha(2,2) = std(normAlpha(:,5))/sqrt(length(normAlpha(:,1)));
NEWSEAlpha(3,2) = std(normAlpha(:,6))/sqrt(length(normAlpha(:,1)));


Fig=figure;
% subplot(2,1,2);
plotNEWRT = barwitherr(NEWSEAlpha(:,:), NEWmeanAlpha(:,:)); % maybe leave 1 at 3rd dim
set(gca, 'XTick',1:3, 'XTickLabel',{ '25' '75' '100' }); ylim([-1 1])
set(plotNEWRT(1),'facecolor','red'); set(plotNEWRT(2),'facecolor','blue');
title('Mean Alpha Lateralization (contra>ipsi) over 500 ms prestimulus period'); xlabel('Spatial Predictability (%)'); ylabel('Mean alpha lateralization (microV)');
legend('Target','Distractor');


%% 75% Val vs Inval Locs

% Total Mean of RTMed
Fig=figure;
% subplot(2,1,2);
plotTotRTMMed75InvalLoc = barwitherr(TotSE75InvalLoc(:,:,1), TotRTm75InvalLoc(:,:,1));
set(gca, 'XTick',1:2, 'XTickLabel',{ 'T75' 'D75' }); ylim([.5 .565])
set(plotTotRTMMed75InvalLoc(1),'facecolor','red'); set(plotTotRTMMed75InvalLoc(2),'facecolor','blue');
title('Group Mean of Median RTs in 75% for val vs inval Location'); xlabel('Stimulus Type'); ylabel('RT (s)');
legend('valid','invalid');


% Total Mean Acc
Fig=figure;
% subplot(2,1,2);
plotTotAcc75InvalLoc = barwitherr(TotAccSE75InvalLoc(:,:,1), TotAccM75InvalLoc(:,:,1));
set(gca, 'XTick',1:2, 'XTickLabel',{ 'T75' 'D75' }); ylim([90 95])
set(plotTotAcc75InvalLoc(1),'facecolor','red'); set(plotTotAcc75InvalLoc(2),'facecolor','blue');
title('Group Mean Accuracies in 75% for val vs inval Location'); xlabel('Stimulus Type'); ylabel('Accuracy (%)');
legend('valid','invalid');

% % Total Mean of Median RTs: Bar Chart
% 
% InvalvsValLoc(:,1) = RTT75all;
% InvalvsValLoc(:,3) = RTD75all;
% InvalvsValLoc(:,2) = RTT75;
% InvalvsValLoc(:,4) = RTD75;
% 
% InvalvsValLocSE(:,1) = std(InvalvsValLoc(:,1))/sqrt(length(Ptt));
% InvalvsValLocSE(:,3) = std(InvalvsValLoc(:,3))/sqrt(length(Ptt));
% InvalvsValLocSE(:,2) = std(InvalvsValLoc(:,2))/sqrt(length(Ptt));
% InvalvsValLocSE(:,4) = std(InvalvsValLoc(:,4))/sqrt(length(Ptt));
% 
% Fig=figure;
% plotInvalvsValLoc = barwitherr(InvalvsValLocSE(:,:), InvalvsValLoc(:,:));
% title('Effects of Invalid Location Removal'); xlabel('Condition'); ylabel('RT (s)');
% set(plotInvalvsValLoc(1),'facecolor','red'); set(plotInvalvsValLoc(2),'facecolor','blue');
% set(gca, 'XTick',1:3, 'XTickLabel',{'25' '75' '100' }); ylim([0.5 0.565])
% legend('T','D')

%% Paired Samples T-Tests (for quick overview)

if ttestOutput == 1;      
    % Stimulus
        % T25 vs T25 to ascertain randomization
        RTT25 = RTMed(:,1,1);
        RTD25 = RTMed(:,4,1);
        [hT25vD25, pT25vD25] = ttest(RTT25, RTD25) % h = H0 or H1; p = p-value

        % T75 vs D75
        RTT75 = RTMed(:,2,1);
        RTD75 = RTMed(:,5,1);
        [hT75vD75, pT75vD75] = ttest(RTT75, RTD75)   

        % T100 vs D100
        RTT100 = RTMed(:,3,1);
        RTD100 = RTMed(:,6,1);
        [hT100vD100, pT100vD100] = ttest(RTT100, RTD100)  
    
    % Spatial Predictability
        
        % T25 vs T75
        RTT25 = RTMed(:,1,1);
        RTT75 = RTMed(:,2,1);
        [hT25vT75, pT25vT75] = ttest(RTT25, RTT75)
    
        % T25 vs T100
        RTT25 = RTMed(:,1,1);
        RTT100 = RTMed(:,3,1);
        [hT25vT100, pT25vT100] = ttest(RTT25, RTT100)
        
        % T75 vs T100
        RTT75 = RTMed(:,2,1);
        RTT100 = RTMed(:,3,1);
        [hT75vT100, pT75vT100] = ttest(RTT75, RTT100)        

        % D25 vs D100
        RTD25 = RTMed(:,4,1);
        RTD100 = RTMed(:,6,1);
        [hD25vD100, pD25vD100] = ttest(RTD25, RTD100)    
                
        % D25 vs D75
        RTD25 = RTMed(:,4,1);
        RTD75 = RTMed(:,5,1);
        [hD25vD75, pD25vD75] = ttest(RTD25, RTD75) 

        % D75 vs D100
        RTD75 = RTMed(:,5,1);
        RTD100 = RTMed(:,6,1);
        [hD75vD100, pD75vD100] = ttest(RTD75, RTD100)
end

% %% Scatterplots across ptts
%
% for j = 1:length(mTrialnT1)
%     mmTrialnT1(j)=mean(mTrialnT1(j,:));
% end
% for j = 1:length(mTrialnT2)
%     mmTrialnT2(j)=mean(mTrialnT2(j,:));
% end
% for j = 1:length(mTrialnT3)
%     mmTrialnT3(j)=mean(mTrialnT3(j,:));
% end
% for j = 1:length(mTrialnD1)
%     mmTrialnD1(j)=mean(mTrialnD1(j,:));
% end
% for j = 1:length(mTrialnD2)
%     mmTrialnD2(j)=mean(mTrialnD2(j,:));
% end
% for j = 1:length(mTrialnD3)
%     mmTrialnD3(j)=mean(mTrialnD3(j,:));
% end
%
% Fig=figure;
% subplot(3,2,1)
% scatter(1:30,mmTrialnT1(1:30))
% title('T1')
% subplot(3,2,3)
% scatter(1:30,mmTrialnT2(1:30))
% title('T2')
% subplot(3,2,5)
% scatter(1:30,mmTrialnT3(1:30))
% title('T3')
% subplot(3,2,2)
% scatter(1:30,mmTrialnD1(1:30))
% title('D1')
% subplot(3,2,4)
% scatter(1:30,mmTrialnD2(1:30))
% title('D2')
% subplot(3,2,6)
% scatter(1:30,mmTrialnD3(1:30))
% title('D3')
%
%
%

% Save
% save(file, 'Bhv');   

