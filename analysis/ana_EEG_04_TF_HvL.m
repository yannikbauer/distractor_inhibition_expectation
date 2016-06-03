% Prestimulus Analysis with TFmaps in fieldtrip
% Alex von Lautz, Yannik Bauer OHBA, 2016

%%%%%%% Todo %%%%%%%

%%%%%%%      %%%%%%%
%% Set up
clear all;
dat_dir =  '/Volumes/Seagate/DistexpEEG/';
resultsdir = '/Users/otte/Desktop/DistExp/results';
subjects=[1:6,8,10:18,20:22];
channels='all';
subcounter=1;
TFwindow = [0 1];
induced=1;

%% Load conditions
load([dat_dir 'details.mat'])   
count=1;
conds=fieldnames(detail);
for i=1:(length(conds)-1)
    if ~any([regexpi(conds{i}, 'Bottom') regexpi(conds{i}, 'Top')])
        rlconds{count}=conds{i};
        count=count+1;
    end
end

%% Loop over subjects
tic
g=waitbar(0,'load');
for isub=subjects;
    
waitbar(isub/length(subjects), g, strcat('Subject Nr.', num2str(isub), ' time elapsed:', num2str(round(toc)), ' seconds'));
    
    subDir = fullfile(dat_dir, detail.sub_details(isub).in_fname{1});
    subDir_out=fullfile('/Volumes/Seagate/DistExpEEG/', detail.sub_details(isub).in_fname{1});
 if ~exist(subDir_out,'dir'), cmd = ['mkdir ' subDir_out]; unix(cmd); end
 
    %% load data
    %fname_in = ['ft' detail.sub_details(isub).in_fname{1} '_fixPost1000.mat'];
    %fname_in = ['pre1000post1000_' num2str(isub)];
    %load(fullfile(subDir,fname_in));
    %load(fullfile(subDir,'condListStimPost1000.mat'))
    %load(fullfile(subDir,'vbadTrials_fixPost1000.mat'))
    %load(fullfile(subDir,strcat(detail.sub_details(isub).in_fname{1}, 'Behav.mat')))
    load([subDir, '/tf_' num2str(isub)]);

    %% Preprocess data
%         cfg              = [];
%         cfg.demean       = 'yes';
%         cfg.detrend      = 'yes';
%         cfg.baselinewindow = [0 2];
%         dat=ft_preprocessing(cfg,dat);
    if ~exist('tfdat')
        cfg              = [];
        %cfg.trial        = ismember(condlist_stim, conds{regexp(Right)}
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 6:1:30;
        cfg.t_ftimwin    = 5./cfg.foi;%
        cfg.toi          = TFwindow(1):.01:TFwindow(2);
        cfg.keeptrials   = 'yes'; % classifiers operate on individual trials
        tfdat=ft_freqanalysis(cfg, dat);
        tfdat.trialinfo = condList_stim;
        tfdat.detail=detail;
        tfdat.vbadtrials=vbadTrials;
        tfdat.rlconds=rlconds;
        tfdat.behav=outputDat.data;

        dat.trialinfo = condList_stim;
        dat.detail=detail;
        dat.vbadtrials=vbadTrials;
        dat.rlconds=rlconds;
        dat.behav=outputDat.data;
        %delete([subDir_out 'tf_' num2str(isub) '.mat'])
        %delete([subDir_out 'pre1000post1000_' num2str(isub) '.mat'])
        %save([subDir_out '/tf_' num2str(isub)], 'tfdat')
        %save([subDir_out '/pre1000post1000_' num2str(isub)], 'dat')
    end

    %% Set up conditions
    
    for i=1:length(rlconds)
        cfg         = [];
        cfg.trials   = setdiff(find(ismember(tfdat.trialinfo, detail.(rlconds{i}))),tfdat.vbadtrials);
        C.(rlconds{i})=ft_selectdata(cfg,tfdat);
        cfg         = [];
        cfg.keeptrials='no';
        M.(rlconds{i})=ft_freqdescriptives(cfg,C.(rlconds{i}));
    end
    %% Subsample trials for tests/statistics
    totrials=[];
    for i=1:length(rlconds)
        totrials(i)=size(C.(rlconds{i}).powspctrm,1);      
    end
    ntrials=min(totrials);
    for i=1:length(rlconds)
        cfg=[];
        cfg.trials = datasample(1:size(C.(rlconds{i}).powspctrm,1), ntrials,'Replace',false);
        C.(rlconds{i}) = ft_selectdata(cfg, C.(rlconds{i}));
    end
    %% Subtract left trials per condition from right right trials
    for i=1:2:length(rlconds)
        cfg = [];
        cfg.parameter = 'powspctrm';
        %cfg.operation = '(x1-x2)/(x1+x2)'; %Normalizing here no good
        cfg.operation = 'subtract';
        MM{isub}.(rlconds{i})= ft_math(cfg, C.(rlconds{i}), C.(rlconds{i+1}));
    end
end

%% Get subject averages
for i=1:2:length(rlconds)
    cfg = [];
    for isub=subjects
        acMM{isub}.(rlconds{i})=ft_freqdescriptives(cfg, MM{isub}.(rlconds{i}));
    end
end
%% Grandaverage
for i=1:2:length(rlconds)
    cfg = [];
    for isub=subjects
        cMM{isub}=acMM{isub}.(rlconds{i});
    end
    GA.rlconds{i}=ft_freqgrandaverage(cfg, cMM{subjects});
end
%%
figure;
cfg = [];
cfg.layout    = 'easycapM11.mat'; % specify the layout file that should be used for plotting
cfg.colorbar = 'yes';
cfg.ylim=[7 13];
cfg.baseline=[0 0.4];
cfg.baselinetype='absolute';
ft_topoplotTFR(cfg, GA.rlconds{11})
%%
eeg_pleft = {'PO3','P3','P1','P5','PO7','O1','P7','TP7'}';
eeg_pright = {'PO4','P4','P2','P6','PO8','O2','P8','TP8'}';
figure;
titles={'T25', 'T75', 'T100', 'D25', 'D75', 'D100'};
for i=1:6
    subplot(2,3,i)
    cfg = [];
    cfg.layout    = 'easycapM11.mat'; % correct cap???
    cfg.colorbar = 'yes';
    cfg.ylim=[7 14];
    cfg.baseline=[0 0.4];
    cfg.baselinetype='absolute';
    cfg.zlim=[-.5 .5];
    cfg.channel={'PO3','P3','P1','P5','PO7','O1','P7','TP7','PO4','P4','P2','P6','PO8','O2','P8','TP8','Oz','POz'};
    ft_topoplotTFR(cfg, GA.rlconds{i*2-1})
    title(titles{i})
end

%% Do statistics
% Multivariate Regression?
% Classification (maybe with 4 locations), then crossvalidation with 10 folds?
% With PCA reduced data - feature selection?

%% Dependent samples Regresssion + cluster correction + Permutation test

cfg=[];
    cfg.method='template';
    cfg.template='/Users/otte/Documents/MATLAB/fieldtrip-20150610/template/neighbours/easycapM11_neighb.mat';
    cfg.layout='easycapM11.lay';
    %cfg.feedback='yes';
    cfg.channel={'PO3','P3','P1','P5','PO7','O1','P7','TP7','PO4','P4','P2','P6','PO8','O2','P8','TP8','Oz','POz'};
    neighbours=ft_prepare_neighbours(cfg, M.(rlconds{1}));    

cfg = [];
%%%% for cluster correction
cfg.neighbours=neighbours;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = .05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
%%%%
%cfg.channel          = '';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesregrT';
cfg.tail             = 0;
cfg.alpha            = .05;   % FWE correction threshold
cfg.numrandomization = 1000; %number of permutations
cfg.channel={'PO3','P3','P1','P5','PO7','O1','P7','TP7','PO4','P4','P2','P6','PO8','O2','P8','TP8','Oz','POz'};
count=1;
udesign=[];
idesign=[];

for isub=subjects
    udesign=[udesign count*ones(size(MM{isub}.(rlconds{1}).powspctrm,1)*3,1)'];
    count=count+1;
end
for isub=subjects
for icond=1:3
    idesign=[idesign icond*ones(size(MM{isub}.(rlconds{1}).powspctrm,1),1)'];
end
end
cfg.design(1,:)  = udesign;
cfg.design(2,:)  = idesign;
cfg.uvar     = 1;                  % row number of design that contains the labels of the UOs (subjects or trials)
cfg.ivar     = 2;                  % row number of the design that contains the labels of the conditions that must be
                                    % compared. The labels are the numbers 1 and 2. 
                                    
count=1;                                    
for isub=subjects   
        T25{count}=MM{isub}.(rlconds{1});
        T75{count}=MM{isub}.(rlconds{3});
        T100{count}=MM{isub}.(rlconds{5});
        D25{count}=MM{isub}.(rlconds{7});
        D75{count}=MM{isub}.(rlconds{9});
        D100{count}=MM{isub}.(rlconds{11});
        count=count+1;
end                                    

stat_T = ft_freqstatistics(cfg, T25{:}, T75{:}, T100{:});
stat_D = ft_freqstatistics(cfg, D25{1:19}, D75{:}, D100{:});

%% Plot results

figure;
cfg = [];
cfg.layout    = 'easycapM11.mat'; % specify the layout file that should be used for plotting
cfg.colorbar = 'yes';
cfg.parameter = 'stat';
cfg.ylim=[7 14];
cfg.baseline=[0 0.4];
%cfg.baselinetype='absolute';
ft_topoplotTFR(cfg, stat_D)

%% Classify high/low expectation SVM + cross validation
O.decoding_start=0;
O.timewins=1:50;
O.decode_steps=.02;
for isub=subjects
   g=waitbar(0,strcat('Decoding Subject Nr.', num2str(isub)));
for itime=O.timewins
    timewindow=[O.decoding_start+itime*O.decode_steps O.decoding_start+(itime+1)*O.decode_steps];

waitbar(itime/length(O.timewins), g, strcat('Decoding Subject Nr.', num2str(isub)));
    cfg=[];

cfg.method  = 'crossvalidate';
cfg.design  = [ones(size(MM{isub}.T25Right.powspctrm,1),1); 2*ones(size(MM{isub}.T100Right.powspctrm,1),1);]';
%cfg.design  = [ones(size(MM{isub}.D25Right.powspctrm,1),1); ones(size(MM{isub}.T25Right.powspctrm,1),1); 2*ones(size(MM{isub}.T100Right.powspctrm,1),1); 2*ones(size(MM{isub}.T75Right.powspctrm,1),1);]';
cfg.latency     = timewindow; % Timewindow of interest
%cfg.freq =2:9; %2:9 is 7-14 hz
cfg.channel={'PO3','P3','P1','P5','PO7','O1','P7','TP7','PO4','P4','P2','P6','PO8','O2','P8','TP8','Oz','POz'};
cfg.statistic = {'accuracy' 'binomial' 'contingency'};
cfg.nfolds = 10;
cfg.layout = 'easycapM11.mat';
%cfg.parameter = 'fourierspctrm';
%cfg.avgovertime = 'yes';
%cfg.avgoverfreq = 'yes';
stat_tf{isub}     = ft_freqstatistics(cfg, MM{isub}.T25Right, MM{isub}.T100Right);
%stat_tf{isub}     = ft_freqstatistics(cfg, MM{isub}.D25Right, MM{isub}.T25Right, MM{isub}.T100Right, MM{isub}.T75Right);
TFstats(isub,(itime*2-1))=stat_tf{1,isub}.statistic.accuracy;
TFstats(isub,(itime*2))=stat_tf{1,isub}.statistic.binomial;

end
delete(g);
end

randomnumber=randi(1000);
save(strcat(resultsdir, 'alpha_D25vsD100_10folds',num2str(randomnumber)),'TFstats', 'O')

t=O.decoding_start+O.decode_steps:O.decode_steps:O.decode_steps*length(O.timewins)+O.decoding_start;
 a=figure; plot(t,smooth(mean(TFstats(subjects,1:2:length(TFstats))),1));
 
 for i=1:2:length(TFstats)
  [H(i),P(i)]=ttest(TFstats(subjects,i)-.50);
   TFstd(i)=std(TFstats(subjects,i))./sqrt(3);
 end
 b=figure; plot(t,P(1:2:length(TFstats)));

 d=figure; shadedErrorBar(t,smooth(mean(TFstats(subjects,1:2:length(TFstats))),3), TFstd(1:2:length(TFstats)))

