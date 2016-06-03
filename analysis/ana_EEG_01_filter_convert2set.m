%% ana_EEG_01_filter_convert2set
% First analysis script. Band-pass filters raw EEG data files and converts to .set.
%
% >>>INPUTS>>>:
% NAME                                  TYPE       DESCRIPTION
% sXX_XX.cnt                            rec file   raw continuous EEG data files
% subjDetails.m                         struct     subject list and bad channels
% 
% <<<OUTPUTS>>>:  
% NAME                                  TYPE       DESCRIPTION
% sXX_XX_filt.set                       struct     filtered EEG data file
% sXX_XX_filt.fdt                       struct     fieltrip-converted data file
%
% @ Attention Group 2016, Oxford Centre for Human Brain Activity, Dept. of Psychiatry,
% Oxford University

%% TODO
% avoid eeg-lab and do compoletely in ft

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

% Load EEG channel configuration
try run(fullfile(dirs.code,'eegChanConfig.m')), % run eegChanConfig.m
catch, error('Script for generation of EEG channel configurations not found!');
end

% Load Subject Details (runs and updates external script)
try run(fullfile(dirs.code,'subjectDetails.m')), % run subjectDetails.m
catch, error('Script for generation of subject details not found!');
end

nSubs = length(subDetails);
doSubs = [1:23];

%% Loop through subject files
fprintf('%%%%%%%%%%%% Starting Filtering and File Conversion %%%%%%%%%%%%\n');

for iSub = doSubs
    dirs.sub = fullfile(dirs.dat,subDetails(iSub).fid{1});
    close all
    ALLEEG = eeglab;
    
    fid.in = [subDetails(iSub).fid{1} '.cnt'];
    fid.out = [subDetails(iSub).fid{1} '_filt.set'];
    if  0%exist(fullfile(dirs.sub,fid.out),'file')
        fprintf('filter already done\n')
        EEG = pop_loadset(fid.out, dirs.sub);
    else
        %% Load raw data
        EEG = pop_loadcnt(fullfile(dirs.sub,fid.in),'dataformat','int32');
        subNum = subDetails(iSub).fid{1};
        EEG.filename = subDetails(iSub).fid{1};
        EEG.setname  = sprintf('S%s',subNum);
        
        % Re-name channels from chan.eegOrig to chan.eeg
        for iChan = 1:length(chan.eegOrig), 
            chanIdx = eeg_chaninds(EEG, chan.eegOrig{iChan});
            EEG.chanlocs(chanIdx).labels = chan.all{iChan};
        end
        
        % Subselect EEG/EOG channels        
        chaneeg = eeg_chaninds(EEG,chan.eeg);
        chaneog = eeg_chaninds(EEG,chan.eog);
        chanall = [chaneeg chaneog]; % cutting out EMG and ref; we should make sure all else is good
        chanref = EEG.data(eeg_chaninds(EEG,chan.ref(2)),:)/2; % because only M2 is saved
        
        EEG.data(chaneeg,:) = EEG.data(chaneeg,:) - repmat(chanref,length(chaneeg),1);        
        EEG.chanlocs = EEG.chanlocs(chanall);
        EEG.data     = EEG.data(chanall,:);
        EEG.nbchan = length(chanall);        
        EEG.setname = sprintf('S%s re-ref',subNum);
        
        % Interpolate bad channels 
        % I am still a bit unsure about this. Maybe OK to use sparingly?? ##YB: WHO is unsure?
        if isempty(subDetails(iSub).badChan)
            pop_eegplot(EEG); % just to bring up the viewer to look for bad channels
            error('badChan not done') % to break out of the script
            % get bad channels, complete in subject details, and re-run
        else
            if strcmp(subDetails(iSub).badChan,'none')
                fprintf('no bad channels\n')
            else
                badChan = eeg_chaninds(EEG,subDetails(iSub).badChan);
                chanlocfile = sprintf('%s/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp', dirs.eeglab);
                EEG = pop_chanedit(EEG, 'lookup', chanlocfile);
                EEG = pop_interp(EEG, badChan, 'spherical');
            end
        end
        
        % Resample to 250 Hz
        if EEG.srate ~= 250
            EEG = pop_resample( EEG, 250);
            EEG.setname = sprintf('S%s resampled',subNum);
            [ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
        end
        
        % Band-pass filter EEG channels
        veog = squeeze(EEG.data(eeg_chaninds(EEG,'VEOG'),:));
        % tmp filter for blink data
        hpfreq = .05; % low-pass frequency (Hz)
        lpfreq = 40; % low-pass frequency (Hz)
        chaneeg = eeg_chaninds(EEG,chan.eeg);
        chaneog = eeg_chaninds(EEG,chan.eog);
        
        eog_dat = EEG.data(chaneog,:);
        eog_dat = eegfilt(eog_dat,EEG.srate,hpfreq,0);
        eog_dat = eegfilt(eog_dat,EEG.srate,0,lpfreq);
        EEG.data(chaneog,:) = eog_dat;
        
        % Avoid low-freq filter artefacts
        hpfreq = .05; % low-pass frequency (Hz)
        lpfreq = 40; % low-pass frequency (Hz)
        
        data = EEG.data(chaneeg,:);
        data = eegfilt(data,EEG.srate,hpfreq,0); % avoid using two filters, use one bandpass
        data = eegfilt(data,EEG.srate,0,lpfreq);
        EEG.data(chaneeg,:) = data;
        
        % Add new useful channels
        % add a smoothed channel for EOG
        veog = gsmooth(veog, 40*EEG.srate/1000); % sd = 40ms 
        dveog = [0 diff(veog)];
        ddveog = [0 0 diff(diff(veog))];
        
        EEG.data(EEG.nbchan+1,:) = dveog';
        EEG.chanlocs(EEG.nbchan+1).labels = 'dveog';        
        EEG.data(EEG.nbchan+2,:) = ddveog';
        EEG.chanlocs(EEG.nbchan+2).labels = 'ddveog';
        EEG.nbchan = size(EEG.data,1);

        pop_saveset(EEG, 'filename', fid.out, 'filepath', dirs.sub);
    end   
end
