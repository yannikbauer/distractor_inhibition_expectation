%% ana_EEG_02_epoch
% Second analysis script. Epochs the filtered and fieldtrip-converted EEG data files.
%
% >>>INPUTS>>>:
% NAME                                  TYPE       DESCRIPTION
% sXX_XX_filt.set                       struct     band-pass filtered EEG data file
% sXX_XX_filt.fdt                       struct     fieltrip-converted band-pass filtered EEG data file
% subjDetails.m                         struct     script creating subject list and bad channels
% 
% <<<OUTPUTS<<<:
% NAME                                  TYPE       DESCRIPTION
% sXX_XX_filt_epFixPost1000.set         struct     epoched EEG data
% sXX_XX_filt_epFixPost1000.fdt         struct     epoched EEG data in fieldtrip format
% sXX_XX_condList_epFixPost1000.mat     struct     list of valid trials w condition triggers
% sXX_XX_badTrials_epFixPost1000.mat    double     list of bad trials (artefacts, blinks) for exclusion
% 
% INFO:
% epoched two sets. One from 0 - 1 (Pre) and 1 - 2 (Post). Alex is Pre, Yannik Post.
% 
% @ Attention Group 2016, Oxford Centre for Human Brain Activity, Dept. of Psychiatry,
% Oxford University

%% Setup
clear all;
close all;
% clc;

% Directories
dirs.dat =  '/Volumes/Y/distractor_inhibition_expectation/EEG/data/eeg'; % input data
dirs.code = '/Users/Yannik/Documents/MATLAB/OX/distractor_inhibition_expectation/analysis/'; % analysis code
dirs.ft = '/Users/Yannik/Documents/MATLAB/Addins/fieldtrip-20160215'; % fieldtrip toolbox
addpath(genpath(dirs.dat), genpath(dirs.code), genpath(dirs.ft));
% dirs.eeglab = '/home/mnoonan/Documents/MATLAB/eeglab12_0_2_0b'; % add eeglab toolbox path
% addpath(genpath(dirs.eeglab));
% rmpath(genpath(['/home/mstokes/matlab/spm8/external/fieldtrip'])); % remove truncated spm8-fieltrip

% Load EEG channel configurations
try run(fullfile(dirs.code,'eegChanConfig.m')), % run eegChanConfig.m
catch, error('Script for generation of EEG channel configurations not found!');
end

% Epoching events using M/EEG triggers for trial start (fixation period)
% We varied ntrials per block, so if it is larger than 37 or very variable you should select 
% events after the behaviour has been loaded and base selEvents on max(ntrial)
selEvents = [2:2:80]; % only even triggers were sent as the acquisition device did not like odd nums
nEvents = length(selEvents);

% Epoching time window
epTime = [-.05 1.05]; % Choose epoch time relative to trigger (trial start)
epName = ['FixPost1000']; % epoch name for file name

% Subject Details (runs and updates external script)
try run(fullfile(dirs.code,'subjectDetails.m')), % run subjectDetails.m
catch, error('Script for generation of subject details not found!');
end
nSubs = length(subDetails);
doSubs = [1];

%% Loop through subject files
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Starting Epoching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
for iSub = doSubs    
    %% Epoch continuous EEG data
    % Directories & file IDs
    fid.in.eeg = [subDetails(iSub).fid{1} '_filt.set'];
    % include epoch info in output name (e.g. 1000ms post and pre fixation = '_filt_epFixPre1000Post1000.set')!
    fid.out.eeg = [subDetails(iSub).fid{1}, '_filt_ep', epName, '.set'];    
    dirs.sub = fullfile(dirs.dat, subDetails(iSub).fid{1});
    % Load experiment task output data
    load(fullfile(dirs.sub, [subDetails(iSub).fid{1} '_outputData.mat']));
    ALLEEG = eeglab;
    
    if  0%exist(fullfile(dirs.sub,fid.out.eeg),'file')
        fprintf('Epoching already done.\n')
        EEG = pop_loadset(fid.out.eeg, dirs.sub);
    else
        EEG = pop_loadset(fid.in.eeg, dirs.sub);  % Load continuous EEG-data
        EEG = pop_epoch(EEG, num2cell(selEvents), epTime); % epoch
        EEG.setname = sprintf('%s epoched',fid.out.eeg);
        pop_saveset(EEG, 'filename', fid.out.eeg, 'filepath', dirs.sub); % Save epoched data
    end
    
    %% Extract M/EEG trigger codes for valid epochs (creates subject condition list containing trial info)
    % Needed because not all epoched trials will go into final analysis
    % Trial start trigger (fixation period) containing trial number within block
    condList.fix = []; % Condition list from fixation triggers
    nEpochs = size(EEG.epoch,2);
    for iEpoch = 1:nEpochs
        condList.fix(iEpoch) = cell2mat(EEG.epoch(iEpoch).eventtype(1)); 
    end
    condList.fix = condList.fix';
    
    % Stimulus location trigger
    condList.stim = [];
    nEpochs = size(EEG.epoch,2);
    for iEpoch = 1:nEpochs
        condList(iEpoch) = cell2mat(EEG.epoch(iEpoch).eventtype(2));
    end
    condList.stim = condList.stim';
    
    fid.out.condList = [subDetails(iSub).fid{1}, '_condList_ep', epName, '.mat'];
    save(fullfile(dirs.sub, fid.out.condList), 'condList')
        
    %% Find bad trials (mainly artefacts and blinks)
    % Automatically makes the most diagnostic traces in addition to manual feedback
    fid.in.bad = fid.out.eeg;
    fid.out.bad = [subDetails(iSub).fid{1}, '_badTrials_ep', epName, '.mat'];
    if exist(fullfile(dirs.sub,fid.out.bad), 'file')
        fprintf('Previous badTrials found.\n')
        load(fullfile(dirs.sub, fid.out.bad)); % load badTrial list        
        EEG = pop_loadset(fid.in.bad, dirs.sub); % load epoched EEG
        tmpEEG = EEG;        
        eeg_data = EEG.data(eeg_chaninds(EEG,chan.eeg),:,:);
        x = reshape(eeg_data,[size(eeg_data,1) size(eeg_data,2)*size(eeg_data,3)]);
        x = zscore(x,[],2);
        y = reshape(x,size(eeg_data));
        
%         z=5;
%         ind = sum((y>z),2)>0;
%         ind = sum(ind)>size(eeg_data,1)/2;
%         ind = squeeze(sum(sum((y>z),1),2));
%         ind = squeeze(sum((mean(abs(y),1)>1),2));
        % tmpEEG.reject.rejmanual = ind>0;
        % tmpEEG.reject.rejmanual(1) = 1;
        % tmpEEG.reject.rejmanualE = zeros(tmpEEG.nbchan,tmpEEG.trials);
        
        VBT = zeros(1,length(tmpEEG.epoch));
        count = 1;
        for i = 1:max(badTrials)
            if i == badTrials(1,count)
                VBT(1,i) = 1;
                count = count+1;
            end
        end
        
        tmpEEG.reject.rejmanual = VBT';
        tmpEEG.reject.rejmanual(1) = 1;
        tmpEEG.reject.rejmanualE = zeros(tmpEEG.nbchan,tmpEEG.trials);
        
        %         tmpEEG.data(eeg_chaninds(tmpEEG,{'HEOG'}),:,:)        = [];
        %         tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'HEOG'})) = [];
        %
        %         tmpEEG.data(eeg_chaninds(tmpEEG,{'VEOG'}),:,:)        = [];
        %         tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'VEOG'})) = [];
        
                tmpEEG.data(eeg_chaninds(tmpEEG,{'dveog'}),:,:)        = [];
                tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'dveog'})) = [];
        
                tmpEEG.data(eeg_chaninds(tmpEEG,{'ddveog'}),:,:)        = [];
                tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'ddveog'})) = [];
        
        tmpEEG.nbchan = size(tmpEEG.data,1);
        pop_eegplot(tmpEEG,1,1,0);
        % pop_eegplot(tmpEEG,1,0,0);
        input('press any key to continue')
        badTrials = find(ALLEEG.reject.rejmanual);
        save(fullfile(dirs.sub,fid.out.bad),'badTrials')
        clear tmpEEG;
        
    else
        EEG = pop_loadset(fid.in, dirs.sub);
        tmpEEG = EEG;
        
        eeg_data = EEG.data(eeg_chaninds(EEG,chan.eeg),:,:);
        x = reshape(eeg_data,[size(eeg_data,1) size(eeg_data,2)*size(eeg_data,3)]);
        x = zscore(x,[],2);
        y = reshape(x,size(eeg_data));
        
        z = 5;
        idx = sum((y>z),2)>0;
        idx = sum(idx)>size(eeg_data,1)/2;
        idx = squeeze(sum(sum((y>z),1),2));
        %ind = squeeze(sum((mean(abs(y),1)>1),2));
        tmpEEG.reject.rejmanual = idx>0;
        tmpEEG.reject.rejmanual(1) = 1;
        tmpEEG.reject.rejmanualE = zeros(tmpEEG.nbchan,tmpEEG.trials);        
        
        %         tmpEEG.data(eeg_chaninds(tmpEEG,{'HEOG'}),:,:)        = [];
        %         tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'HEOG'})) = [];
        %
        %         tmpEEG.data(eeg_chaninds(tmpEEG,{'VEOG'}),:,:)        = [];
        %         tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'VEOG'})) = [];
        
                tmpEEG.data(eeg_chaninds(tmpEEG,{'dveog'}),:,:)        = [];
                tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'dveog'})) = [];
        
                tmpEEG.data(eeg_chaninds(tmpEEG,{'ddveog'}),:,:)        = [];
                tmpEEG.chanlocs(eeg_chaninds(tmpEEG,{'ddveog'})) = [];
              
        tmpEEG.nbchan = size(tmpEEG.data,1);
        pop_eegplot(tmpEEG,1,1,0);
        input('press any key to continue')
        badTrials = find(ALLEEG.reject.rejmanual);
        save(fullfile(dirs.sub,fid.out.bad),'badTrials')
        clear tmpEEG;
    end
    
end

