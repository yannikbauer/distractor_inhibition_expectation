%% ana_EEG_03_
% Third analysis script. Automatic and manual artefact removal. 
%
% >>>INPUTS>>>:
% NAME                                  TYPE       DESCRIPTION
% sXX_XX_filt_epFixPost1000.set         struct     epoched EEG data
% sXX_XX_filt_epFixPost1000.fdt         struct     epoched EEG data in fieldtrip format
% sXX_XX_badTrials_epFixPost1000.mat    double     list of bad trials (artefacts, blinks) for exclusion
% subjDetails.m                         struct     script creating subject list and bad channels
% condList.mat                          struct     script creating list of condition trigger codes
% <<<OUTPUTS<<<:
%  ...
% 
% @ Attention Group 2016, Oxford Centre for Human Brain Activity, Dept. of Psychiatry,
% Oxford University

%% Setup
clear all;
close all;
clc;

% Directories
dirs.dat =  '/Volumes/Y/distractor_inhibition_expectation/EEG/data/eeg'; % input data
dirs.code = '/Users/Yannik/Documents/MATLAB/OX/distractor_inhibition_expectation/analysis/';
addpath(dirs.dat, genpath(dirs.code));
if ~exist(dirs.dat, 'dir'), error('EEG data directory not found!'), end

dirs.eeglab = '/home/mnoonan/Documents/MATLAB/eeglab12_0_2_0b'; % add eeglab toolbox path
addpath(genpath(dirs.eeglab));
rmpath(genpath(['/home/mstokes/matlab/spm8/external/fieldtrip'])); % remove truncated spm8-fieltrip
dirs.fieldtrip = '/home/benmc/matlab/fieldtrip-20130507';
addpath(genpath(dirs.fieldtrip));

% Load EEG channel configurations
try run(fullfile(dirs.code,'eegChanConfig.m')), % run eegChanConfig.m
catch, error('Script for generation of EEG channel configurations not found!');
end

% Select trigger events
% currently epoching from fixation - if ntrials is larger than 37 or very variable you should select 
% events after the Behav. has been loaded and base on max(ntrial)
selecEvents = [2:2:80]; 
nEvents = length(selecEvents);

% Subject Details (runs and updates external script)
try run(fullfile(dirs.code,'subjectDetails.m')),
catch, error('Script for generation of subject details not found!');
end
nsubs=length(subDetails);
doSubs = [1:4];

%% Loop through subject files
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Starting Artefact Removal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
for iSub = doSubs    
    dirs.sub = fullfile(dirs.dat,subDetails(iSub).in_fname{1});
    fid.ft.in = [subDetails(iSub).fid{1} '_filt_epFixPost1000.set'];   
    fid.ft.out = [subDetails(iSub).fid{1} '_filt_epFixPost1000_ft.mat'];    
    
    if 0%exist(fullfile(dirs.sub,fid.out),'file')
        fprintf('convert trials to fieldtrip already done\n')
    else
        cfg = [];
        cfg.dataset = fullfile(dirs.sub, fid.ft.in);
        cfg.continuous = 'no';
        cfg.baselinewindow = [-.1 .0];
        cfg.detrend = 'no';
        dat = ft_preprocessing(cfg);
        save(fullfile(dirs.sub, fid.ft.out), 'dat');
    end
    
    % Find blinks in epochs (store vector)
    fid.ica.in = fid.ft.out;
    fid.ica.out = [subDetails(iSub).fid{1} '_filt_epFixPost1000_ft_ica.mat'];
    if 0%exist(fullfile(dirs.sub,fid.ica.out),'file')
        fprintf('ICA already done.\n')
        load(fullfile(dirs.sub,fid.ica.out));
    else
        % Load EEG file and subject specific bad trials
        load(fullfile(dirs.sub,fid.ica.in));
        load(fullfile(dirs.sub,[subDetails(iSub).fid{1},'_badTrials_epFixPost1000.mat']));
        
        allTrials = 1:size(dat.trial,2);       
        nTrials =  setdiff(allTrials,badTrials); % Exclude bad trials
        
        cfg = [];
        cfg.method = 'runica';
        cfg.channel = chan.eeg;
        cfg.trials = nTrials;
        comp = ft_componentanalysis(cfg, dat);
        comp.nTrials = nTrials;
        comp.badTrials = badTrials;
        save(fullfile(dirs.sub,fid.ica.out),'comp');
   end
    
end % separate file here for processing of upper part that takes time, 
% then bottom half requires input, so comment that out and do later.
    
    
%% Independent Component Analysis to remove artefacts
% requires manual input
    fid.icaClean.in = fid.ica.out;
    fid.icaClean.in2 = fid.ft.out;
    fid.out = [subDetails(iSub).fid{1} '_filt_epFixPost1000_ft_ica_clean.mat'];
    if 0%exist(fullfile(dirs.sub,fid.out),'file')
        fprintf('ica clean already done\n')
    else
        load(fullfile(dirs.sub,fid.icaClean.in));
        load(fullfile(dirs.sub,fid.icaClean.in2));
        load(fullfile(dirs.sub,[subDetails(iSub).fid{1},'_badTrials_epFixPost1000.mat']));
        
        allTrials = 1:size(dat.trial,2);
        nTrials =  setdiff(allTrials,vbadTrials);
        
        ndat = dat;
        time = {};
        trial = {};
        for iTrial = 1:length(nTrials)
            time{iTrial}=dat.time{nTrials(iTrial)}; 
            trial{iTrial}=dat.trial{nTrials(iTrial)};  
        end
        ndat.sampleinfo = dat.sampleinfo(nTrials,:);
        ndat.time=time;
        ndat.trial=trial;
        dat = ndat; clear ndat;
        
        nComp = size(comp.trial{1},1);
        tempR = zeros(nComp,1);
        comp_dat = zeros(nComp,length(comp.trial),length(comp.time{1}));
        for iTrial=1:length(comp.trial), comp_dat(:,iTrial,:) = comp.trial{iTrial}; end % put in the time in col 2
        
        % Blinks
        veog_dat = zeros(length(dat.trial),length(dat.time{1}));
        for iTrial=1:length(dat.trial), veog_dat(iTrial,:) = dat.trial{iTrial}(ismember(dat.label,'VEOG'),:); end
        for c=1:nComp
            tmpdat = squeeze(comp_dat(c,:,:));
            x=veog_dat(:); z=tmpdat(:);
            r = corrcoef(x,z);
            tempR(c) = r(2);
        end       
        [s si]=sort(abs(tempR),'descend');
        
        cfg = [];
        cfg.component          = si(1:20);
        cfg.layout = 'easycapM11.lay';
        cfg.comment = 'Blink';
        cfg.marker = 'labels';
        figure(2),ft_topoplotIC(cfg, comp);
        
        % Saccades
        heog_dat = zeros(length(dat.trial),length(dat.time{1}));
        for iTrial=1:length(dat.trial), heog_dat(iTrial,:) = dat.trial{iTrial}(ismember(dat.label,'HEOG'),:); end
        for c=1:nComp
            tmpdat = comp_dat(c,:,:);
            r = corrcoef(heog_dat(:),tmpdat(:));
            tempR(c) = r(2);
        end
        [s si]=sort(abs(tempR),'descend');
        
        cfg = [];
        cfg.component          = si(1:20);
        cfg.layout = 'easycapM11.lay';
        cfg.comment = 'Saccades';
        cfg.marker = 'labels';
        figure(3),ft_topoplotIC(cfg, comp);
        
        % Remove component
        input('press any key to continue'); 
        badComp = input('Which comps?','s');
        badComp = str2num(badComp);
        
        cfg = [];
        cfg.component = unique(badComp);
        dat = ft_rejectcomponent(cfg, comp);
        dat.ntrials = comp.ntrials;
        dat.vbadTrials = comp.vbadTrials;
        dat.bad_comp = badComp;
        dat.comp = comp;
        dat.heog = heog_dat;
        dat.veog = veog_dat;
        save(fullfile(dirs.sub,fid.out),'dat');
    end    
    
   