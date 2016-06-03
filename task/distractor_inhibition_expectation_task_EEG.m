%% distractor_inhibition_expectation_task_EEG
%
% Introduction: This task explores the neural mechanisms of experience-based distractor inhibition
% via the build-up of expectations arising from statistical regularities in stimulus presentations.
% It specifically investigates the role of oscillatory alpha power in this mechanism.

% Task: Visual Discrimination task on chequerboard-like stimulus target patches
% appearing in one of four quadrants. Participants have to indicate whether
% the chequerboard has many small or a few large squares. The likelihood of
% appearance in one location is manipulated independently for both targets
% and distracters at 25%, 75% and 100%, to test independent effects of
% target facilitation and distracter inhibition and their relation to each
% other.

% Script info: configured for EEG; includes main task + training task (for stimulus decoding model)

% @ Attention Group 2016, Oxford Centre for Human Brain Activity, Dept. of Psychiatry, Oxford University
% Authors: Yannik Bauer, Horst Alexander von Lautz

%% TO DO:

%% Administrative Setup
clear all
close all

echo off; % Prevents MATLAB from reprinting the source code when the program runs

isEEG = true; % true if on stim PC
isMEG = false;
isEyeTrack = false;

%% Random number generation & control
rng('default'); % resets seed to default (Mersenne Twister with seed=0) to control randomization
rng('shuffle'); % sets random seed to ensure random number generation
seed = rng; % saves state for later recall using rng(seed)

% For later use of escape-key to quit expt (see below)
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');

%% Input Subject Details
% Asks for input in prompt window
subjDetails.subjNumber = input('Subject number (e.g. 01): ','s'); % 's' = string input
if isempty(subjDetails.subjNumber)
   subjDetails.subjNumber = 'test'; % Gives Subj Name 'test' if SubjDetails.Name is empty
end
subjDetails.SeriesNum = input('Subject series number (e.g. 01): ','s');
subjDetails.Age = input('Subject Age: ','s');
subjDetails.Sex = input('Subject Sex (m/f): ','s');
subjDetails.Hand = input('Subject Handedness (r/l): ','s');

%% Directories, Files and Output data
% Path to the parallel port to EEG
addpath('C:\Documents and Settings\mstokes\My Documents\MATLAB\NTPort');

% Data storage directory
dirs.data = fullfile(pwd,'data');
if ~exist(dirs.data,'dir'), mkdir(dirs.data), end % Creates new data folder if folder does not yet exist

% File name setup in specified directory
timeNow = fix(clock); % Returns time for use in subjectFileDetails
timeNow = [num2str(timeNow(4)) '-' num2str(timeNow(5)) '-' num2str(timeNow(6))];
fid.subjDet = ['s', subjDetails.subjNumber, '_', subjDetails.SeriesNum, '__', date '_' timeNow];
fid.output = fullfile(dirs.data,[fid.subjDet,'.mat']);

% Create output data matrix (details on output matrix towards end)
output.data = []; % Creates empty sub-matrix in output into which later data are inserted
output.SubjDetails = subjDetails; % Creates subject details sub-file
output.seed = seed; % records seed to replicate trial conditions in case of data loss
output.train = []; % creates empty training data matrix

try
    %%  Psychtoolbox Setup
    % For explanation of use of try-loop see end.
    % Note that some of the administrative setup commands could have been set up in the first Administrative Setup.
    % Included here are the commands specific to PTB (rather than general matlab ones).
    
    AssertOpenGL; % Makes sure this is running on OpenGL Psychtoolbox
    
    % Set screen preferences
    oldSyncLevel = Screen('Preference', 'SkipSyncTests', 0); % old = 2; set to 1 to skip (NB: skipping gives bad timings!)
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 0);
    
    % Settings: % ##YB: change for debugging/expt
    % These are inserted here (rather than earlier) so that the user still
    % is in control of the script during the SubjDetails prompt
    HideCursor; % Hides cursor
    FlushEvents; % Remove events from the system event queue % ##YB: what events?
    ListenChar(2); % Suppresses output of keypresses in the command window.
    
    if isEEG||isMEG||isEyeTrack,
        IOPortfolder = 'C:\Users\Attention_Lab\Documents\MATLAB\IOPort';
        addpath(IOPortfolder);
        [portobject, portaddress] = OpenIOPort;
        
        triggerPulseWidth = 0.05; % Send trigger for 50 ms
        triggerHoldValue = 0; 
    end
    
    %% WINDOW SETUP
    screenid = max(Screen('Screens')); % Selects screen with maximum id for output window
    
    % Open a fullscreen, onscreen window with grey background. Enable 32bpc
    % floating point framebuffer via imaging pipeline on it, if this is possible
    % on your hardware while alpha-blending is enabled. Otherwise use a 16bpc
    % precision framebuffer together with alpha-blending. We need alpha-blending
    % here to implement the nice superposition of overlapping gabors. The expt will
    % abort if your graphics hardware is not capable of any of this (see ProceduralGaborDemo).
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    [win, winRect] = PsychImaging('OpenWindow', screenid, 128);
    
    Priority(MaxPriority(win)); % confer max priority to PTB-screen
    
    [winWidth, winHeight] = Screen('WindowSize', win); % Function to retrieve width and height of window
    [xCentre, yCentre] = WindowCenter(win); % Function to retrieve window centre
    
    % Alpha-blending
    % NB: will be turned on only for overlapping Gabor section, and back to normal afterwards
    % For overlapping Gabors (messes with Fixation circle etc.):
    %         Screen('BlendFunction', win, GL_ONE, GL_ONE);
    % Back to normal:
    %         Screen('BlendFunction', win, GL_ONE, GL_ZERO);
    
    % To ensure accurate consistent timings:
    ifi = Screen('GetFlipInterval', win, 100, 50e-6, 10);
    % All timings below will use following procedure:
    % flipTime = previousTime + (round(interval/ifi)-0.5*ifi)
    % where round(interval/ifi) gets rounded integer of #flipsPerInterval,
    % prepared half a flip beforehand (to ensure buffer is ready),
    % and multiplied by ifi to ouput the TIME
    % eg. tCueOnset = Screen('Flip', win, tStartFixation+(round(0.200/0.0167)-0.5)*0.0167)
    % = tStartFixation+11.5*0.0167 = tStartFixation+.192;
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL PARAMETER SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Color Setup
    black = BlackIndex(win); % Retrieves the CLUT code for black.
    white = WhiteIndex(win); % Retrieves the CLUT code for white.
    grey = (black+white)/2; % Computes the CLUT code for grey.
    darkGrey = (black+grey)/2;
    lightGrey = (white+grey)/2;
    %             blue = [0 153 255];
    %             orange = [255 102 0];
    %             red = [200 0 0];
    %             green = [0 181 0];
    %             brown = [138 138 0];
    
    %% Stimulus Parameters
    % Fixation Cross
    fixCrossDiameter = 11;
    
    % Gabor Patch
    width = 200;
    height = 200;
    TTilt = 0; % orientation angle in degrees (0-360)
    phase = 90; % phase of the gabors sine grating in degrees
    sc = 25; % spatial constant of the gaussian hull function of the gabor,
    % ie. the "sigma" value in the exponential function
    freq = [0.0345, 0.0275]; % spatial frequency in cycles per pixel & old=[0.0345, 0.0275]
    % new=old-0.001, old+0.001
    contrast = 14; % amplitude of your gabor in intensity units (old: 1.5)
    aspectratio = 1.0; % aspect ratio of the hull of the gabor
    
    % Four possible positions;
    stimPos = [200 -200; 200 200; -200 200; -200 -200]; % ensure it subtends 4° vis angle
    % Redefine Stimulus position relative to centre
    stimPos(:,1) = stimPos(:,1) - width./2 + xCentre;
    stimPos(:,2) = stimPos(:,2) - height./2 + yCentre;
    
    % Call CreateProceduralGabor (done here to load methods into memory;
    % not actually strictly necessary at this point)
    % NB: for more info, check "help CreateProceduralGabor" and "Screen('DrawTexture?')" or ProceduralGaborDemo
    [gaborTex, gaborrect] = CreateProceduralGabor(win, width, height, 0, [0 0 0 0], 0, 1);
    
    %% Timing Parameters (in seconds)
    startFixationDur = 1;
    stimulusDur = 0.200; % old=0.080
    ITI = 0.2; % inter-trial interval (ITI)
    
    %% Response Key/Button Allocation
    rng('shuffle'); % shuffle seed for random response allocation
    responseAllocation = randi(2);
    rng(seed); % immediately reset it to initial state to ensure same randomization
    
    %% Introduction Text Setup
    % Text Screen Parameters
    Screen('TextSize', win, 20);
    Screen('TextFont', win, 'Times New Roman');
    Screen('TextStyle', win, 0); % 0=normal,1=bold,2=italic,4=underline,8=outline,32=condense,64=extend
    
    % Introduction page 1
    introTxt1 = ['Welcome to the experiment!'...
        '\nYour TARGETS are round patches with a chequered pattern appearing in one of four possible locations. '...
        'Your task is to discriminate their square size and respond using the keyboard:'];
    DrawFormattedText(win, introTxt1, 'center', yCentre-80, black, 80, [], [], 1.2);
    
    switch responseAllocation
        case 1
            DrawFormattedText(win, '\n\n\n\n\nPress "c" for SMALL SQUARES and "m" for BIG SQUARES.',...
                'center', yCentre-80, black, 80, [], [], 1.2);
        case 2
            DrawFormattedText(win, '\n\n\n\n\nPress "c" for BIG SQUARES and "m" for SMALL SQUARES.',...
                'center', yCentre-80, black, 80, [], [], 1.2);
    end
    
    introTxt2 = ['You will get audio-feedback for each response:',...
        '\n HIGH tones for correct and LOW tones for incorrect responses'];
    DrawFormattedText(win, introTxt2, 'center', yCentre+100, black, 80, [], [], 1.2);    
    DrawFormattedText(win, 'Please press any key to continue reading the intro...', 'center', yCentre+300, black);
    
    % Present Two Sample Target Stimuli
    % NB: explanations on Gabors further below    
    Screen('BlendFunction', win, GL_ONE, GL_ONE);    
    TdstRect1 = OffsetRect(gaborrect, stimPos(4,1), stimPos(4,2));
    Screen('DrawTexture', win, gaborTex, [], TdstRect1, TTilt, [], [], [], [], kPsychDontDoRotation,...
        [phase, freq(1), sc, contrast, aspectratio, 0, 0, 0]);
    Screen('DrawTexture', win, gaborTex, [], TdstRect1, TTilt-90, [], [], [], [], kPsychDontDoRotation,...
        [phase, freq(1), sc, contrast, aspectratio, 0, 0, 0]);
    
    TdstRect2 = OffsetRect(gaborrect, stimPos(1,1), stimPos(1,2));
    Screen('DrawTexture', win, gaborTex, [], TdstRect2, TTilt, [], [], [], [], kPsychDontDoRotation,...
        [phase, freq(2), sc, contrast, aspectratio, 0, 0, 0]);
    Screen('DrawTexture', win, gaborTex, [], TdstRect2, TTilt-90, [], [], [], [], kPsychDontDoRotation,...
        [phase, freq(2), sc, contrast, aspectratio, 0, 0, 0]);    
    Screen('Flip', win);
    WaitSecs(.5);
    KbWait;
    
    % Introduction page 2
    introTxt3 = ['There will also be a DISTRACTER, a round patch of black-white stripes '...
        'appearing in another location. TRY TO IGNORE THE DISTRACTER.'];
    DrawFormattedText(win, introTxt3, 'center', yCentre-110, black, 80, [], [], 1.2);
    DrawFormattedText(win, 'Please press any key to continue reading the intro...', 'center', yCentre+300, black);
    
    % Present Sample Target and Distractor Stimuli
    TdstRect1 = OffsetRect(gaborrect, stimPos(1,1), stimPos(1,2));
    Screen('DrawTexture', win, gaborTex, [], TdstRect1, TTilt, [], [], [], [], kPsychDontDoRotation,...
        [phase, freq(2), sc, contrast, aspectratio, 0, 0, 0]);
    Screen('DrawTexture', win, gaborTex, [], TdstRect1, TTilt-90, [], [], [], [], kPsychDontDoRotation,...
        [phase, freq(2), sc, contrast, aspectratio, 0, 0, 0]);
    
    Dtilts = [.5:16].*180/16;
    Dtiltselect = randi(16);
    DTilt=Dtilts(Dtiltselect);
    
    DdstRect1 = OffsetRect(gaborrect, stimPos(2,1), stimPos(2,2));
    Screen('DrawTexture', win, gaborTex, [], DdstRect1, DTilt, [], [], [], [], kPsychDontDoRotation,...
        [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);
    Screen('DrawTexture', win, gaborTex, [], DdstRect1, DTilt, [], [], [], [], kPsychDontDoRotation,...
        [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);
    Screen('DrawTexture', win, gaborTex, [], DdstRect1, DTilt, [], [], [], [], kPsychDontDoRotation,...
        [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);
    Screen('DrawTexture', win, gaborTex, [], DdstRect1, DTilt, [], [], [], [], kPsychDontDoRotation,...
        [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);    
    Screen('BlendFunction', win, GL_ONE, GL_ZERO);
    
    % Fixation Cross
    Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
    Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);    
    Screen('Flip', win);
    WaitSecs(.5);
    KbWait;
    
    % Introduction page 3
    introTxt4 = ['Some tips:'...
        '\nAlways focus on the CENTRAL FIXATION CROSS.'...
        '\nRELAX and blink as little as possible.'...
        '\n(Otherwise, the trial will be invalid.)'...
        '\n\nAim for about 90% accuracy.'...
        '\n\n\n\n\n\nYou will get a BREAK about every 7 minutes.'];
    DrawFormattedText(win, introTxt4, 'center', yCentre-300, black, 80, [], [], 1.2);
    DrawFormattedText(win, 'Press any key to start the experiment. \nGOOD LUCK!', 'center', yCentre+250, black);
    
    % Fixation Cross
    Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
    Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);    
    Screen('Flip', win);
    WaitSecs(.5);
    KbWait;
    
    %% Condition Factorization Setup
    % This ensures combination of all factors (location, tilt...) for all conditions in specific blocks,
    % executed in the (trial-)for-loop created below
    
    % BLOCK CONDITIONS
    % here: 1[nDisplaySetSize] x 2[nConditions] x 3[nProbManip] x 2[blockReps] x4[nlocations]
    %       = 48 blocks = 58 mins = ca. 2.6 secs per trial)
    nDisplaySetSize = 1; % currently T+D only (but left here for future manipulations)
    nStimulusManip = 2; % T manipulation, D manipulation
    nProbability = 3; % Makes three probability conditions
    nblockReps = 2; % number of repetitions of whole set of conditions
    nlocations = 4; % number stimulus locations
    
    blockTypeCnt = 1; % sets blockType count to one;
    for iDisplaySetSize = 1:nDisplaySetSize % T or T+D
        for iCond = 1:nStimulusManip
            for iProbability = 1:nProbability
                for iRepeat = 1:nblockReps
                    for iLocation = 1:nlocations;
                        if iProbability == 1
                            locationProb = 0.25; % Probability of stimulus appearance in given location
                        elseif iProbability == 2
                            locationProb = 0.75;
                        elseif iProbability == 3
                            locationProb = 1;
                        end
                        blockType(blockTypeCnt,1) = iDisplaySetSize; % set size
                        blockType(blockTypeCnt,2) = iCond; % 2 conditions
                        blockType(blockTypeCnt,3) = iProbability;
                        blockType(blockTypeCnt,4) = locationProb;
                        blockType(blockTypeCnt,5) = iLocation;
                        blockTypeCnt = blockTypeCnt + 1;
                    end
                end
            end
        end
    end    
    blockType = blockType(randperm(size(blockType,1)),:); % randomize block order
    
    % TRIALS: This makes 30±2 trials within each block
    nFreq = length(freq); % = spatial frequency of T (big or small)
    trialReps = 1; % 1 repetition =~ 7 mins (for 144 trials)
    mBlocklength = 30; % mean length of one block
    SDBlocklength = 2; % Standard deviation of blocklength
    for blockLengthIteration = 1:length(blockType);
        % This is the function for mean±SD blocklength
        Blocklength(blockLengthIteration) = mBlocklength+round(randn*SDBlocklength);
    end
    
    trialTypeCnt = 1; % sets trialType counter for all trials to one;
    for iBlock = 1:length(blockType)
        for  iTrialType = 1:Blocklength(iBlock)
            if blockType(iBlock,5) == 1 % Makes location matrices
                locations = [1 2 3 4];
            elseif blockType(iBlock,5) == 2
                locations = [2 1 3 4];
            elseif blockType(iBlock,5) == 3
                locations = [3 1 2 4];
            elseif blockType(iBlock,5) == 4
                locations = [4 1 2 3];
            end
            
            % Set actual location to fixed location with given probability
            value = rand; % random number between 0 and 1
            if value < blockType(iBlock,4) % If the random is smaller than the location probability
                locProb=1; % use fixed location (= first in matrix)
            else
                locProb=randi(3)+1; % use one random location of the three others
            end
            
            % Set locations of controlled stimulus and other stimulus
            if blockType(iBlock,2) == 1
                trialType(trialTypeCnt,1) = locations(locProb);
                % To ensure that D is once in each of the 3 locations where T is not:
                other_locations = [1:4];
                other_locations = other_locations(locations(locProb)~=other_locations);
                trialType(trialTypeCnt,2) = other_locations(randi(3));
            elseif blockType(iBlock,2) == 2
                trialType(trialTypeCnt,2) = locations(locProb);
                % To ensure that D is once in each of the 3 locations where T is not:
                other_locations = [1:4];
                other_locations = other_locations(locations(locProb)~=other_locations);
                trialType(trialTypeCnt,1) = other_locations(randi(3));
            end
            
            trialType(trialTypeCnt,3) = freq(randi(2));
            trialType(trialTypeCnt,4) = blockType(iBlock,2); % T or D manipulation
            trialType(trialTypeCnt,5) = blockType(iBlock,3); % probability manipulation type
            trialType(trialTypeCnt,6) = blockType(iBlock,4); % locationProbability
            trialType(trialTypeCnt,7) = blockType(iBlock,5); % location
            trialType(trialTypeCnt,8) = iBlock; % Blocknumber
            trialType(trialTypeCnt,9) = iTrialType; % Trialnumber in block
            trialType(trialTypeCnt,10) = value; % value of i
            
            trialTypeCnt = trialTypeCnt+1;
        end
    end
    % trialType = trialType(randperm(size(trialType,1)),:); % randomize trial order
    
    %% Start Experiment
    tStartExpt = GetSecs;
    
    % Set all counters to 1 (NB: do NOT put into trial-loop, otherwise blockCounter will be reset to 1)
    trialCnt = 1; % trial counter
    trainTrialCnt = 1; % training trial counter
    
    %% Start trials
    for iTrial = 1:length(trialType)
        %% Give error message when accuracy is low
        if trialCnt < 20
            % Checking every 5 trials that for the last 4 not more than 2 mistakes were made.
            if mod(trialCnt,5) == 0 && sum(output.data((trialCnt-4):(trialCnt-1),14))<3
                
                errorTxt = ['You have been making mistakes lately. Please ensure that you respond correctly:'...
                    '\nYou may take more time to respond or take a short break if you want'...
                    '\n\n\nPress any key to continue'];
                DrawFormattedText(win, errorTxt, 'center', yCentre-200, black);
                % Fixation Cross
                Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                
                Screen('Flip', win);
                WaitSecs(1);
                KbWait;
            end
        end
        
        if trialCnt > 20
            % Checking every 5 trials that for the last 20 not more than 4 mistakes were made.
            if mod(trialCnt,10) == 0 && sum(output.data((trialCnt-20):(trialCnt-1),14))<15
                
                errorTxt = ['You have been making mistakes lately. Please ensure that you respond correctly:'...
                    '\nYou may take more time to respond or take a short break if you want'...
                    '\n\n\nPress any key to continue'];
                DrawFormattedText(win, errorTxt, 'center', yCentre-200, black);
                % Fixation Cross
                Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                
                Screen('Flip', win);
                WaitSecs(1);
                KbWait;
            end
        end
        
        %% Take a break every x trials & Current Scores
        % About 8 mins - every 7 blocks. (6 breaks)
        
        if trialCnt > 1
            %with mod(y,x) specify with 'x' after how many blocks breaks will happen
            if mod(trialType(trialCnt,8),7) == 0 && trialType(trialCnt,9) == 1
                
                DrawFormattedText(win, 'Break Time! \nPress any key to continue', 'center', yCentre-200, black);
                % Fixation Cross
                Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                
                % Participant Current Score Counter: Accuracy and RTs                
                % calculate Accuracy
                accScore = sum(output.data(:,14))/ length(output.data(:,14));
                % calculate mean RTs
                meanRTScore = mean(output.data(:, 22));
                % Display Text & Vars
                DrawFormattedText(win, ['Your current scores:'...
                    '\n\n Percent Correct: ' num2str(round(accScore*100)) '%'...
                    '\n Mean Reaction Time: ' num2str(meanRTScore) 's'], 'center', yCentre+200, black);
                
                Screen('Flip', win);
                WaitSecs(3);
                KbWait;
                
                %% Training setup
                trainCnt = 1;
                nLocs = 4;
                nTilts = 16;
                nReps = 2;
                for iRep = 1:nReps
                    for iLoc = 1:nLocs
                        for iTilt = 1:nTilts
                            trainType(trainCnt,1) = iLoc;
                            trainType(trainCnt,2) = iTilt;
                            trainCnt = trainCnt+1;
                        end
                    end
                end
                TrainTypeRnd = trainType(randperm(size(trainType,1)),:);
                
                %% Training intro                
                DrawFormattedText(win, ['After breaks, we ask you to do a small and simple intermission task.'...
                    '\n\n It includes only one stimulus.'...
                    'Please indicate on which SIDE of the screen the stimulus appears.'...
                    '\n\n Press C for LEFT, and M for RIGHT.'], 'center', yCentre-300, black);
                DrawFormattedText(win, 'Press any key to continue', 'center', yCentre+200, black);
                
                Screen('Flip', win);
                WaitSecs(3);
                KbWait;
                
                for iTrain = 1:length(TrainTypeRnd);% start training task loop
                    
                    %% M/EEG Trigger Codes (sent to recording computer)
                    meegTrigStartTrain = 98; % M/EEG trigger code: training trial start
                    meegTrigLocTrain = 244+2*TrainTypeRnd(iTrain,1); % M/EEG trigger code: training stimulys location
                    
                    %% Training presentation
                    % lightGrey cross marks new trial
                    Screen('DrawLine',win,lightGrey,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                    Screen('DrawLine',win,lightGrey,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                    tTrainOffset = Screen('Flip',win);
                    
                    % black cross appears 1 second before stimuli
                    Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                    Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                    
                    % TIME: Set and get start fixation cross onset time
                    tTrainStartFixation = Screen('Flip', win, tTrainOffset+(ITI+0.1)); % NO RAND => 300 ms ITI
                                                            
                    % Send M/EEG trigger code
                    if isEEG
                        io64(portobject, portaddress, meegTrigStartTrain); % This sends the trigger
                        WaitSecs(triggerPulseWidth);
                        io64(portobject, portaddress, triggerHoldValue); % This sets the trigger channel back to 0
                    end
                    
                    %% Present Stimuli
                    Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                    Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                    
                    % Enable alpha-blending for overlapping Gabors 
                    % (disable afterwards, as it messes up fixationCircle etc.)
                    Screen('BlendFunction', win, GL_ONE, GL_ONE);
                    
                    % Distractor Gabor
                    DxPos = stimPos(TrainTypeRnd(iTrain,1),1);
                    DyPos = stimPos(TrainTypeRnd(iTrain,1),2);
                    DdstRect = OffsetRect(gaborrect, DxPos, DyPos);
                    
                    % Draws distractor Gabor with random tilt and spatial frequency that is
                    % half between the two Target spat freqs
                    % NB: two Gabors to equal contrast intensity of Target Gabor
                    Dtilts = [.5:16].*180/16;
                    DTilt = Dtilts(TrainTypeRnd(iTrain,2));
                    Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
                        kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast,...
                        aspectratio, 0, 0, 0]);
                    Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
                        kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast,...
                        aspectratio, 0, 0, 0]);
                    Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
                        kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast,...
                        aspectratio, 0, 0, 0]);
                    Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
                        kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast,...
                        aspectratio, 0, 0, 0]);
                    
                    % TIME: Set & get Stimulus Onset Time
                    % ##AvL: Change of timing, in script before times of cues here
                    tTrainStimOnset = Screen('Flip', win, ...
                        tTrainStartFixation+(round(startFixationDur/ifi)-0.5)*ifi);
                    
                    % Send second M/EEG trigger code
                    if isEEG
                        io64(portobject, portaddress, meegTrigLocTrain); % This sends the trigger
                        WaitSecs(triggerPulseWidth);
                        io64(portobject, portaddress, triggerHoldValue); % This sets the trigger channel back to 0
                    end
                    
                    % Disable alpha-blending
                    Screen('BlendFunction', win, GL_ONE, GL_ZERO);
                    
                    %%  P Response & Audio-Feedback
                    % ##YB: get response before offset of stimuli?                    
                    % Present fixation circle
                    Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                    Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                    
                    % TIME: Set & get response period onset
                    tTrainResponseOnset = Screen('Flip',win, tTrainStimOnset+(round(stimulusDur/ifi)-0.5)*ifi);
                    
                    % Get participant responses & give audio feedback
                    keyIsDown = 0; % set keyboard buttons to 0
                    while ~keyIsDown
                        [keyIsDown, tResponse, keyCode] = KbCheck;                                                
                        if keyIsDown
                            % Press 'escape' to quit experiment
                            if keyCode(escapeKey);
                                Screen('CloseAll');
                                ShowCursor;
                                ListenChar(0);
                                error('Aborted using escape key.')
                                
                                % If no escape is pressed, proceed as usual
                            else
                                RTTrainOnset = tResponse - tTrainStartFixation; % Gets RT from Trial Onset
                                RTTrainStimulusOnset = tResponse - tTrainStimOnset; % Gets RT from Stimulus Onset
                                
                                % Feedback according to random response key/button allocation
                                switch TrainTypeRnd(iTrain,1) % location
                                    case {3 4} % left
                                        if strcmp(KbName(keyCode),'c') % if 'c' is pressed;
                                            kbResponse = 1; % assign num value rather than char for easier output ana
                                            trainResponseAccuracy = 1; % Correct
                                            beep = MakeBeep(900,.05);
                                            Snd('Open');
                                            Snd('Play',beep);
                                        elseif strcmp(KbName(keyCode),'m') % elseif 'm' is pressed
                                            kbResponse = 2;
                                            trainResponseAccuracy = 0; % Incorrect
                                            beep = MakeBeep(500,.05);
                                            Snd('Open');
                                            Snd('Play',beep);
                                        else keyIsDown = 0; % keep waiting for 'c' or 'm' response
                                        end
                                        
                                    case  {1 2} % right
                                        if strcmp(KbName(keyCode),'m')
                                            kbResponse = 2;
                                            trainResponseAccuracy = 1; % Correct
                                            beep = MakeBeep(900,.05);
                                            Snd('Open');
                                            Snd('Play',beep);
                                        elseif strcmp(KbName(keyCode),'c')
                                            kbResponse = 1;
                                            trainResponseAccuracy = 0; % Incorrect
                                            beep = MakeBeep(500,.05);
                                            Snd('Open');
                                            Snd('Play',beep);
                                        else keyIsDown = 0; % keep waiting for 'c' or 'm' response
                                        end
                                end % switch TrainTypeRnd(iTrain,1) % location                                
                            end % if keyCode(escapeKey)
                        end % if keyIsDown
                    end % while ~keyIsDown
                    
                    output.train(trainTrialCnt,1) = TrainTypeRnd(iTrain,1); % location
                    output.train(trainTrialCnt,2) = TrainTypeRnd(iTrain,2); % tilt
                    output.train(trainTrialCnt,3) = kbResponse; % Participant Response ('c' = 1; 'm' = 2)
                    output.train(trainTrialCnt,4) = trainResponseAccuracy; % response Accuracy
                    output.train(trainTrialCnt,5) = tTrainOffset; % End of last trial = beginning of ITI period
                    output.train(trainTrialCnt,6) = tTrainStartFixation; % Fixation Start
                    output.train(trainTrialCnt,7) = tTrainStartFixation-tTrainOffset; % measured ITI from previous trial
                    output.train(trainTrialCnt,8) = tTrainStimOnset; % Stimulus Onset
                    output.train(trainTrialCnt,9) = tTrainStimOnset-tTrainStartFixation; % measured SOA
                    output.train(trainTrialCnt,10) = tTrainResponseOnset; % Response Period Onset
                    output.train(trainTrialCnt,11) = tTrainResponseOnset-tTrainStimOnset; % measured stimulusDur
                    output.train(trainTrialCnt,12) = RTTrainStimulusOnset; % RT from stimulus onset
                    output.train(trainTrialCnt,13) = RTTrainOnset; % RT from trial onset (just double-checking timings)
                    output.train(trainTrialCnt,14) = meegTrigStartTrain; % M/EEG trigger code: training trial start
                    output.train(trainTrialCnt,15) = meegTrigLocTrain; % M/EEG trigger code: training stimulus location
                                        
                    save(fid.output, 'output'); % Save trial output data into file
                    trainTrialCnt = trainTrialCnt+1;                    
                end % for iTrain = 1:length(TrainTypeRnd);% start training task loop
                
                %% Back to task
                DrawFormattedText(win, 'Now back to the main task. Remember: ', 'center', yCentre-200, black);
                
                switch responseAllocation
                    case 1
                        DrawFormattedText(win, 'Press "c" for SMALL SQUARES and "m" for BIG SQUARES.', ...
                            'center', yCentre-80, black, 80, [], [], 1.2);
                    case 2
                        DrawFormattedText(win, 'Press "c" for BIG SQUARES and "m" for SMALL SQUARES.', ...
                            'center', yCentre-80, black, 80, [], [], 1.2);
                end
                
                DrawFormattedText(win, 'Press any key to continue', 'center', yCentre+200, black);
                % Fixation Cross
                Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
                Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
                
                Screen('Flip', win);
                WaitSecs(3);
                KbWait;                
            end
        end
        
        %% M/EEG Trigger Codes
        % M/EEG Trigger for Trial Start Fixation (codes trial # within block)
        meegTrigStart = trialType(iTrial,9)*2; % NB: x2 because EEG-amp doesn't like odd #s
        
        % M/EEG Trigger for Stimulus Presentation - codes location of T & D
        if trialType(trialCnt,4) == 1 && trialType(trialCnt,5) == 1; %T25
            if trialType(trialCnt,1) == 1 %T in loc 1
                if trialType(trialCnt,2)==2 % D in loc 2
                    meegTrigLoc=102; 
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=104;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=106;
                end
            elseif trialType(trialCnt,1) == 2% T in loc 2
                if trialType(trialCnt,2)==1 % D in loc 1
                    meegTrigLoc=108;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=110;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=112;
                end
            elseif trialType(trialCnt,1) == 3
                if trialType(trialCnt,2)==1
                    meegTrigLoc=114;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=116;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=118;
                end
            elseif trialType(trialCnt,1) == 4
                if trialType(trialCnt,2)==1
                    meegTrigLoc=120;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=122;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=124;
                end
            end
        elseif trialType(trialCnt,4) == 1 && trialType(trialCnt,5)==2;%T75
            if trialType(trialCnt,1) == 1 % T in loc 1
                if trialType(trialCnt,2)==2 % D in loc 2
                    meegTrigLoc=126;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=128;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=130;
                end
            elseif trialType(trialCnt,1) == 2
                if trialType(trialCnt,2)==1
                    meegTrigLoc=132;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=134;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=136;
                end
            elseif trialType(trialCnt,1) == 3
                if trialType(trialCnt,2)==1
                    meegTrigLoc=138;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=140;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=142;
                end
            elseif trialType(trialCnt,1) == 4
                if trialType(trialCnt,2)==1
                    meegTrigLoc=144;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=146;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=148;
                end
            end
        elseif trialType(trialCnt,4) == 1 && trialType(trialCnt,5)==3;%T100
            if trialType(trialCnt,1) == 1 % T in loc 1
                if trialType(trialCnt,2)==2 % D in loc 2
                    meegTrigLoc=150;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=152;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=154;
                end
            elseif trialType(trialCnt,1) == 2
                if trialType(trialCnt,2)==1
                    meegTrigLoc=156;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=158;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=160;
                end
            elseif trialType(trialCnt,1) == 3
                if trialType(trialCnt,2)==1
                    meegTrigLoc=162;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=164;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=166;
                end
            elseif trialType(trialCnt,1) == 4
                if trialType(trialCnt,2)==1
                    meegTrigLoc=168;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=170;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=172;
                end
            end
        elseif trialType(trialCnt,4) == 2 && trialType(trialCnt,5)==1;%D25
            if trialType(trialCnt,1) == 1 % T in loc 1 (i.e. same coding principle as above)
                if trialType(trialCnt,2)==2 % D in loc 2
                    meegTrigLoc=174;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=176;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=178;
                end
            elseif trialType(trialCnt,1) == 2
                if trialType(trialCnt,2)==1
                    meegTrigLoc=180;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=182;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=184;
                end
            elseif trialType(trialCnt,1) == 3
                if trialType(trialCnt,2)==1
                    meegTrigLoc=186;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=188;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=190;
                end
            elseif trialType(trialCnt,1) == 4
                if trialType(trialCnt,2)==1
                    meegTrigLoc=192;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=194;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=196;
                end
            end
        elseif trialType(trialCnt,4) == 2 && trialType(trialCnt,5)==2;%D75
            if trialType(trialCnt,1) == 1 % T in loc 1
                if trialType(trialCnt,2)==2 % D in loc 2
                    meegTrigLoc=198;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=200;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=202;
                end
            elseif trialType(trialCnt,1) == 2
                if trialType(trialCnt,2)==1
                    meegTrigLoc=204;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=206;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=208;
                end
            elseif trialType(trialCnt,1) == 3
                if trialType(trialCnt,2)==1
                    meegTrigLoc=210;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=212;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=214;
                end
            elseif trialType(trialCnt,1) == 4
                if trialType(trialCnt,2)==1
                    meegTrigLoc=216;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=218;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=220;
                end
            end
        elseif trialType(trialCnt,4) == 2 && trialType(trialCnt,5)==3;%D100
            if trialType(trialCnt,1) == 1 % T in loc 1
                if trialType(trialCnt,2)==2 % D in loc 2
                    meegTrigLoc=222;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=224;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=226;
                end
            elseif trialType(trialCnt,1) == 2
                if trialType(trialCnt,2)==1
                    meegTrigLoc=228;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=230;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=232;
                end
            elseif trialType(trialCnt,1) == 3
                if trialType(trialCnt,2)==1
                    meegTrigLoc=234;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=236;
                elseif trialType(trialCnt,2)==4
                    meegTrigLoc=238;
                end
            elseif trialType(trialCnt,1) == 4
                if trialType(trialCnt,2)==1
                    meegTrigLoc=240;
                elseif trialType(trialCnt,2)==2
                    meegTrigLoc=242;
                elseif trialType(trialCnt,2)==3
                    meegTrigLoc=244;
                end
            end
        end
        
        %% Present ITI fixation cross and then trial start fixation cross
        % NB: inserted here rather than end of trial loop to enable recording of ITI
        
        % lightGrey cross marks new trial
        Screen('DrawLine',win,lightGrey,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
        Screen('DrawLine',win,lightGrey,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);        
        tTrialOffset = Screen('Flip',win);
        
        % black cross appears 1 second before stimuli
        Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
        Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
        
        % TIME: Set and get start fixation cross onset time
        tStartFixation = Screen('Flip', win, tTrialOffset+(ITI+randi(4)*0.1)); % For randomized ITI: 300-600 ms
        
        %Send M/EEG trigger code
        if isEEG
            io64(portobject, portaddress, meegTrigStart); % This sends the trigger
            WaitSecs(triggerPulseWidth);
            io64(portobject, portaddress, triggerHoldValue); % This sets the trigger channel back to 0
        end
        
        %% Present Stimuli
        
        Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
        Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
        
        % Enable alpha-blending for overlapping Gabors (disable afterwards, as it messes up fixationCircle etc.)
        Screen('BlendFunction', win, GL_ONE, GL_ONE);
        
        % Stimulus Position
        TxPos = stimPos(trialType(iTrial,1),1);
        TyPos = stimPos(trialType(iTrial,1),2);
        TdstRect = OffsetRect(gaborrect, TxPos, TyPos);
        
        % Draws two overlapping Gabors to create plaid stimulus with spatial frequency
        % specified in (trialType(trialLoop,3))
        Screen('DrawTexture', win, gaborTex, [], TdstRect, TTilt, [], [], [], [],...
            kPsychDontDoRotation, [phase, (trialType(iTrial,3)), sc, contrast, aspectratio, 0, 0, 0]);
        Screen('DrawTexture', win, gaborTex, [], TdstRect, TTilt-90, [], [], [], [],...
            kPsychDontDoRotation, [phase, (trialType(iTrial,3)), sc, contrast, aspectratio, 0, 0, 0]);
        
        % Distractor Gabor
        DxPos = stimPos(trialType(iTrial,2),1);
        DyPos = stimPos(trialType(iTrial,2),2);
        DdstRect = OffsetRect(gaborrect, DxPos, DyPos);
        
        % Draws distractor Gabor with random tilt and spatial frequency that is half between
        % the two Target spat freqs
        % NB: two Gabors to equal contrast intensity of Target Gabor
        Dtilts = [.5:16].*180/16;
        Dtiltselect = randi(16);
        DTilt=Dtilts(Dtiltselect);
        Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
            kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);
        Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
            kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);
        Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
            kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);
        Screen('DrawTexture', win, gaborTex, [], DdstRect, DTilt, [], [], [], [],...
            kPsychDontDoRotation, [phase, (freq(1)-(freq(1)-freq(2))/2), sc, contrast, aspectratio, 0, 0, 0]);
        
        % TIME: Set & get Stimulus Onset Time
        % ##AvL: Change of timing, in script before times of cues here
        tStimOnset = Screen('Flip', win, tStartFixation+(round(startFixationDur/ifi)-0.5)*ifi);
        
        % Send M/EEG location trigger
        if isEEG
            io64(portobject, portaddress, meegTrigLoc); % This sends the trigger
            WaitSecs(triggerPulseWidth);
            io64(portobject, portaddress, triggerHoldValue); % This sets the trigger channel back to 0
        end
        
        % Disable alpha-blending
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        
        %%  Participant Response & Audio-Feedback
        % ##YB: get response before offset of stimuli?        
        % Present fixation circle
        Screen('DrawLine',win,black,xCentre,yCentre-fixCrossDiameter,xCentre,yCentre+fixCrossDiameter,2.0);
        Screen('DrawLine',win,black,xCentre-fixCrossDiameter,yCentre,xCentre+fixCrossDiameter,yCentre,2.0);
        
        % TIME: Set & get response period onset
        tResponseOnset = Screen('Flip',win, tStimOnset+(round(stimulusDur/ifi)-0.5)*ifi);        
        
        % Get participant responses & give audio feedback
        keyIsDown = 0; % set keyboard buttons to 0
        while ~keyIsDown
            [keyIsDown, tResponse, keyCode] = KbCheck;            
            if keyIsDown
                % Press 'escape' to quit experiment
                if keyCode(escapeKey);
                    Screen('CloseAll');
                    ShowCursor;
                    ListenChar(0);
                    error('Aborted using escape key.')                                       
                else % If no escape is pressed, proceed as usual
                    RTTrialOnset = tResponse - tStartFixation; % Gets RT from Trial Onset
                    RTStimOnset = tResponse - tStimOnset; % Gets RT from Stimulus Onset
                    
                    % Feedback according to random response key/button allocation
                    switch responseAllocation                        
                        case 1
                            switch trialType(iTrial,3) % Spatial Frequency
                                case freq(1) % Small Pattern
                                    if strcmp(KbName(keyCode),'c') % if 'c' is pressed;
                                        kbResponse = 1; % assign num value rather than char for easier output analysis
                                        responseAccuracy = 1; % Correct
                                        beep = MakeBeep(900,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    elseif strcmp(KbName(keyCode),'m') % elseif 'm' is pressed
                                        kbResponse = 2;
                                        responseAccuracy = 0; % Incorrect
                                        beep = MakeBeep(500,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    else keyIsDown = 0; % keep waiting for 'c' or 'm' response
                                    end                                    
                                case freq(2) % Big Pattern
                                    if strcmp(KbName(keyCode),'m')
                                        kbResponse = 2;
                                        responseAccuracy = 1; % Correct
                                        beep = MakeBeep(900,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    elseif strcmp(KbName(keyCode),'c')
                                        kbResponse = 1;
                                        responseAccuracy = 0; % Incorrect
                                        beep = MakeBeep(500,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    else keyIsDown = 0; % keep waiting for 'c' or 'm' response
                                    end
                            end                            
                        case 2
                            switch trialType(iTrial,3) % Spatial Frequency
                                case freq(1) % Small Pattern
                                    if strcmp(KbName(keyCode),'c') % 'c'-press
                                        kbResponse = 1;
                                        responseAccuracy = 0; % Incorrect
                                        beep = MakeBeep(500,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    elseif strcmp(KbName(keyCode),'m')
                                        kbResponse = 2;
                                        responseAccuracy = 1; % Correct
                                        beep = MakeBeep(900,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    else keyIsDown = 0; % keep waiting for 'c' or 'm' response
                                    end                                    
                                case freq(2) % Big Pattern
                                    if strcmp(KbName(keyCode),'m') % 'm'-press
                                        kbResponse = 2;
                                        responseAccuracy = 0; % Incorrect
                                        beep = MakeBeep(500,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    elseif strcmp(KbName(keyCode),'c')
                                        kbResponse = 1;
                                        responseAccuracy = 1; % Correct
                                        beep = MakeBeep(900,.05);
                                        Snd('Open');
                                        Snd('Play',beep);
                                    else keyIsDown = 0; % keep waiting for 'c' or 'm' response
                                    end
                            end % switch trialType(iTrial,3) % Spatial Frequency
                    end % switch responseAllocation
                end % if keyCode(escapeKey);
            end % if keyIsDown
        end % while ~keyIsDown        
        tEndExpt = GetSecs;
        
        %%  Fill Output Matrix        
        % Trial & Block info
        output.data(trialCnt,1)  = trialType(iTrial,8);  output.hdr(1,1) = {'block number (1-30)'};
        output.data(trialCnt,2)  = trialType(iTrial,9);  output.hdr(1,2) = {'trial number within block (20±4)'};
        output.data(trialCnt,3)  = trialType(iTrial,5);  output.hdr(1,3) = {'probability Manipulation'};
        output.data(trialCnt,4)  = trialType(iTrial,6);  output.hdr(1,4) = {'probability'};
        output.data(trialCnt,5)  = trialType(iTrial,4);  output.hdr(1,5) = {'expected Stimulus manipulation (T or D)'};
        output.data(trialCnt,6)  = trialType(iTrial,7);  output.hdr(1,6) = {'expected location'};
        output.data(trialCnt,7)  = trialType(iTrial,1);  output.hdr(1,7) = {'actual T position'};
        output.data(trialCnt,8)  = trialType(iTrial,2);  output.hdr(1,8) = {'actual D position'};
        output.data(trialCnt,9)  = DTilt;                output.hdr(1,9) = {'random D Tilt (0-180)'};
        output.data(trialCnt,10) = trialType(iTrial,10); output.hdr(1,10) = {'value (testing if script works as supposed)'};
        output.data(trialCnt,11) = trialType(iTrial,3);  output.hdr(1,11) = {'T Spatial Frequency: 1 = small; 2 = big'};
        output.data(trialCnt,12) = responseAllocation;   output.hdr(1,12) = {'resp key alloc: 1=[L=SMALL,R=BIG square];2=opposite'};
        % Measurements: responses & timings
        output.data(trialCnt,13) = kbResponse;           output.hdr(1,13) = {'participant Response ("c" = 1; "m" = 2)'};
        output.data(trialCnt,14) = responseAccuracy;     output.hdr(1,14) = {'response accuracy'};
        output.data(trialCnt,15) = tTrialOffset;         output.hdr(1,15) = {'end of last trial = beginning of ITI period'};
        output.data(trialCnt,16) = tStartFixation;       output.hdr(1,16) = {'fixation Start'};
        output.data(trialCnt,17) = tStartFixation-tTrialOffset; output.hdr(1,17) = {'measured ITI from previous trial'};
        output.data(trialCnt,18) = tStimOnset;           output.hdr(1,18) = {'stimulus Onset'};
        output.data(trialCnt,19) = tStimOnset-tStartFixation;   output.hdr(1,19) = {'measured SOA'};
        output.data(trialCnt,20) = tResponseOnset;       output.hdr(1,20) = {'response period onset'};
        output.data(trialCnt,21) = tResponseOnset-tStimOnset;   output.hdr(1,21) = {'measured stimulusDur'};
        output.data(trialCnt,22) = RTStimOnset;          output.hdr(1,22) = {'RT from stimulus onset'};
        output.data(trialCnt,23) = RTTrialOnset;         output.hdr(1,23) = {'RT from trial onset (just double-check timings)'};
        output.data(trialCnt,24) = tEndExpt-tStartExpt;  output.hdr(1,24) = {'Total expt duration'};
        output.data(trialCnt,25) = meegTrigStart;        output.hdr(1,25) = {'M/EEG trigger: trial start'};
        output.data(trialCnt,26) = meegTrigLoc;          output.hdr(1,26) = {'M/EEG trigger: stimuli locations'};
        
        %%  Save Output Before Next Trial                
        save(fid.output, 'output'); % Save trial output data into file                
        trialCnt = trialCnt+1; % Add 1 to trial counter
        
    end % trialLoop              
    
    %% Set up End of Experiment    
    % Set End Screen
    Screen('TextSize', win, 50);
    DrawFormattedText(win, 'Thank you!\n You are great!','center', yCentre-80, black);
    Screen('Flip', win);
    
    % Wait for Keyboard response and close all screens
    Waitsecs(2);
    KbWait;
    Screen('CloseAll');
    
    % Close Parallel Port
    if isEEG||isMEG||isEyeTrack, CloseParPort;
    end
    
    % Restore general settings
    Priority(0);
    FlushEvents;
    ShowCursor;
    ListenChar(0);
    
    % Restore preferences
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Screen('Preference', 'SkipSyncTests', oldSyncLevel);
        
catch % catch for error-handling    
    % ---------------- Error Handling & Clean-up ----------------
    % If there is an error in our code, we will end up here.
    % The try-catch block ensures that Screen will restore the display and return us
    % to the MATLAB prompt even if there is an error in our code.  Without this try-catch
    % block, Screen could still have control of the display when MATLAB throws an error, in
    % which case the user will not see the MATLAB prompt.
    
    % Save output & rethrow error
    save(fid.output, 'output');
    psychrethrow(psychlasterror); % Rethrow the error again so the user sees the error description.
    
    % Close Parallel Port
    if isEEG||isMEG||isEyeTrack, CloseParPort;
    end
    
    % Restore Settings
    Screen('CloseAll');
    
    Priority(0);
    FlushEvents;
    ShowCursor;
    ListenChar(0);
    
    % Restore preferences
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Screen('Preference', 'SkipSyncTests', oldSyncLevel);        
end % end of try-catch loop
