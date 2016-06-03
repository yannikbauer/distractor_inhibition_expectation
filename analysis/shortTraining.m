% Part of the alpha Quadrant Encoding script

%close all
clc;
isub=1;
if(0)
       setfile = [pwd '/es05_01_filt_fixTrain.set'];
       cfg.dataset = setfile; 
       cfg.continuous = 'no';
       data = ft_preprocessing(cfg)
       save ft_es05_01_filt_fixTrain.mat data cfg
    end
    ftfile = '/home/alexvl/Documents/Randomblock/s05_01/ft_es05_01_filt_fixTrain.mat';
    load(ftfile,'data');
%%
    %load(filename);
    %waitbar(0,hbar,sprintf('preprocessing S%02d',sublist(isub)));
    
    % Preprocessing Settings
    base_single = false; %remove prestimulus baseline, separately for each sensor and trial? only for ERPs
    use_norm    = false; %normalize data point-by-point? this will remove the ERP and ERP-related variance, so probably only good for regression
    norm_single = true;  %only active if use_norm=true

    chanlist = 1:60;
    nchans   = length(chanlist);

    fsample = data.fsample;

    dtime = fix([-0.300 +1.996]*fsample); %time window of interest, relative to cue onset
    dbase = fix([-0.200 +0.000]*fsample); %baseline period, relative to cue onset
    dstep = 1;

    ftype = 'gaussian';      % smooth the data with a gaussian filter (for regression analysis, not ERP/TF)
    ftype = [];              % don't filter the data
    fsize = 0.032*fsample;   % 32 ms = 8 samples at 250 Hz
    timelist = (dtime(1):dstep:dtime(2))/fsample;
    ntimes = length(timelist);

    % re-cut epochs around probe onset and possibly remove baseline
    %waitbar(0,hbar,sprintf('processing S%02d...',sublist(isub)));
    
    % Filter out bad trials
%     ntrials = setdiff(data.trial, vBadTrials);
%     
    
    ntrials = length(data.trial);
    datanew = data;
    dataerp_elem = nan(ntrials,nchans,ntimes);
    for itrial = 1:ntrials
        if ~isempty(ftype)
            data.trial{itrial} = filtfast(data.trial{itrial},2,[],ftype,fsize);
        end
        ionset = find(data.time{itrial} <= 0,1,'last');
        itime  = ionset+(dtime(1):dstep:dtime(2));
        itime  = min(itime,size(data.trial{itrial},2));
        dataerp_elem(itrial,:,:) = data.trial{itrial}(chanlist,itime);
        if base_single
            ibase = ionset+(dbase(1):dbase(2));
            dataerp_elem(itrial,:,:) = squeeze(dataerp_elem(itrial,:,:))-repmat(mean(data.trial{itrial}(:,ibase),2),[1,ntimes]);
        end
        datanew.trial{itrial} = squeeze(dataerp_elem(itrial,:,:));
        datanew.time{itrial}  = timelist;
    end
    data = datanew; 
    data.label = data.label(chanlist);
    clear datanew dataerp_elem
    
    % Add zero-padding
    add_padding = 1; %use zero-padding to avoid TF edge artifacts?
    if add_padding
        pad = ceil([1.000 1.000]*fsample)/fsample; %1  'PO3' 'PO4' 'P3' 'P4'sec should be plenty for alpha (8 Hz = 125 ms, at 5 cycles per wavelet = 750 ms)
        nchans = size(data.trial{1},1);
        for itrial = 1:ntrials
            data.trial{itrial} = [zeros(nchans,pad(1)*fsample),data.trial{itrial},zeros(nchans,pad(2)*fsample)];
            data.time{itrial}  = linspace(data.time{itrial}(1)-pad(1),data.time{itrial}(end)+pad(2),size(data.trial{itrial},2));
        end
    end
    %
    
    remove_erp = false; 
    if remove_erp % remove phase-locked activity for TF analysis? Probably not a great idea...
        cc    = [ 0 0 +1; ...% Retro-Cued Left
                  1 0 +1; ...% Retro-Cued Right
                  0 1 +1; ...% Pre-Cued Left
                  1 1 +1; ...% Pre-Cued Right
                  0 0 -1; ...% Neutral, Probed Left
                  1 0 -1];   % Neutral, Probed Right

        for ct   = 1:size(cc,1)
            icue = [];
            if cc(ct,3) == 1
                icue = (usable & csid == cc(ct,1) & cued == cc(ct,2) & cues == cc(ct,3));
            else
                icue = (usable & psid == cc(ct,1) & cued == cc(ct,2) & cues == cc(ct,3));
            end
            icue    = reshape(find(icue),1,[]);
            dataerp = [];
            if remove_erp
                dataerp = mean(cat(3,data.trial{icue}),3);
                for itrial = icue
                    data.trial{itrial} = data.trial{itrial}-dataerp;
                end
            end
        end
        clear dataerp
    end
    
    flip_data = false;
    if flip_data
        ichanorig = [01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35];
        ichanflip = [02 01 07 06 05 04 03 12 11 10 09 08 17 16 15 14 13 22 21 20 19 18 27 26 25 24 23 32 31 30 29 28 35 34 33];
       for itrial = 1:ntrials
           if csid(itrial) == 1
              data.trial{itrial} = data.trial{itrial}(ichanflip,:); 
           end
       end
    end    
    %waitbar(0,hbar,sprintf('processing S%02d',sublist(isub)));
    
    tfmeth = 'hanning'; % time-frequency analysis method (multitapering or hanning or wavelet)
    chanoi = {};
    tons = -0.300;
    toff = +2.000;
    timeoi = data.time{1}(data.time{1} >= tons & data.time{1} <= toff); % time samples of interest (s)
    dstep  = 1;
    timeoi = timeoi(1:dstep:end);
    triallen = timeoi(end)-timeoi(1);%window length in sec
    freqoi = [5:1:20];
    
    switch lower(tfmeth)
        
        case 'multitapering' % time-frequency analysis via multitapering (k = 2*tw*fw-1)
            cfg            = [];
            cfg.method     = 'mtmconvol';
            cfg.output     = 'pow';
            cfg.keeptrials = 'yes';
            cfg.feedback   = 'none';
            cfg.channel    = chanoi;
            cfg.foi        = freqoi;
            cfg.toi        = timeoi;
            cfg.taper      = 'dpss';
            cfg.t_ftimwin  = 5./cfg.foi; % time window (s)
            cfg.tapsmofrq  = 0.25*cfg.foi; % frequency window (Hz)
            
        case 'hanning' % time-frequency analysis via Hanning tapers (fw = 1/tw)
            cfg            = [];
            cfg.method     = 'mtmconvol';
            cfg.output     = 'pow';
            cfg.keeptrials = 'yes';
            cfg.feedback   = 'none';
            cfg.channel    = chanoi;
            cfg.foi        = freqoi;
            cfg.toi        = timeoi;
            cfg.taper      = 'hanning';
            cfg.pad        = 7;%make sure this is an integer larger than the max (padded) trial length, in secs
            cfg.t_ftimwin  = 5./cfg.foi; % time window (s)
            
        case 'wavelet' % time-frequency analysis via Morlet wavelets (fw = f/w and tw = 1/fw)
            cfg            = [];
            cfg.method     = 'wavelet';
            cfg.output     = 'pow'; % fourier to get complex values, pow to get power
            cfg.keeptrials = 'yes';
            cfg.feedback   = 'none';
            cfg.channel    = chanoi;
            %cfg.channel    = 'all';
            cfg.foi        = freqoi;
            cfg.toi        = timeoi;
            cfg.width      = 5; % wavelet width
            
        otherwise
            error('invalid time-frequency analysis method!');
            
    end
    tfr = ft_freqanalysis(cfg,data); %clear data    
    eeg_chans = tfr.label;
    log_transform = 1;
    if log_transform
        tfr.powspctrm = 10*log10(tfr.powspctrm);
    end   
    %%  Train on Neutral Probe Location, Test on Cued Location
    
    load('/home/alexvl/Documents/Randomblock/s05_01/s05_01Behav.mat')
    prbloc = outputDat.train(:,1);
    
    if(1)
        tic
        fit_behavior = false;
        add_offset   = true;
        demean       = false;
        do_neutral   = false;
        tfrcurr = tfr;
        ichans       = ismember(tfrcurr.label,eeg_chans);
        %cues: -1=neutral, 1=cued, cued: 0=retro, 1=pre, vldt: 1=valid, 0=invalid
        ul = find(prbloc == 4);
        ll = find(prbloc == 3);
        lr = find(prbloc == 2);
        ur = find(prbloc == 1);
        
        totvec = cat(1,ul,ll,lr,ur);
        ind    = cat(1,ones(length(ul),1)*1,ones(length(ll),1)*2,ones(length(lr),1)*3,ones(length(ur),1)*4);
        
        foi = find(tfrcurr.freq>=8&tfrcurr.freq<=8); %train at all freq
        toi = find(tfrcurr.time>1.050&tfrcurr.time<1.5);
        tfrcurr.powspctrm(isnan(tfrcurr.powspctrm)) = 0;
        
        dat1=squeeze(mean(mean(tfrcurr.powspctrm(ul,ichans,foi,toi),3),4));
        dat2=squeeze(mean(mean(tfrcurr.powspctrm(ll,ichans,foi,toi),3),4));
        dat3=squeeze(mean(mean(tfrcurr.powspctrm(lr,ichans,foi,toi),3),4));
        dat4=squeeze(mean(mean(tfrcurr.powspctrm(ur,ichans,foi,toi),3),4));
        
        D = [dat1; dat2; dat3; dat4]';
        if demean
            Dmu = mean(D,2);
            D  = D - repmat(Dmu,1,size(D,2));
        end
        
        C = [[ones(size(dat1,1),1)*1 ones(size(dat1,1),1)*0 ones(size(dat1,1),1)*0 ones(size(dat1,1),1)*0];...
            [ones(size(dat2,1),1)*0 ones(size(dat2,1),1)*1 ones(size(dat2,1),1)*0 ones(size(dat2,1),1)*0];...
            [ones(size(dat3,1),1)*0 ones(size(dat3,1),1)*0 ones(size(dat3,1),1)*1 ones(size(dat3,1),1)*0];...
            [ones(size(dat4,1),1)*0 ones(size(dat4,1),1)*0 ones(size(dat4,1),1)*0 ones(size(dat4,1),1)*1]]';
        if add_offset
            C(5,:) = 1;
        end
        
        %test on all trials
        ulTest = find(prbloc == 4);
        llTest = find(prbloc == 3);
        lrTest = find(prbloc == 2);
        urTest = find(prbloc == 1);
        
        indTest = [ones(size(ulTest,1),1)*1;ones(size(llTest,1),1)*2;ones(size(lrTest,1),1)*3;ones(size(urTest,1),1)*4];        
        itrl   = [ulTest;llTest;lrTest;urTest];
        ntest  = length(itrl);
        nfreq  = length(find(tfrcurr.freq>=5&tfrcurr.freq<=20));
        ntime  = length(tfrcurr.time);
        fprintf('\nRetrocues:\n');
        C2 = [];
        for itrial = 1:ntest
            fprintf('.');
            itest  = itrl(itrial);
            itrain = setdiff(1:length(totvec),find(itest==totvec));
            D1     = D(:,itrain);
            C1     = C(:,itrain);
            
            W = D1*C1'*pinv(C1*C1');
            pinvW = pinv(W'*W)*W'; %compute here for speed
            
            for f=1:nfreq        
                %D2 = squeeze(mean(mean(tfrcurr.powspctrm(itest,ichans,f,t),3),4))';
                D2 = permute(tfrcurr.powspctrm(itest,ichans,f,:),[2 4 1 3]);
                if demean
                    D2  = D2 - repmat(Dmu,1,size(D2,2));
                end
                C2(:,itrial,f,:)  = pinvW*D2;
            end
        end
        fprintf('\n');
        resT1(isub,:,:,1,:) = permute([mean(C2(1,indTest==1,:,:),2) mean(C2(2,indTest==2,:,:),2) mean(C2(3,indTest==3,:,:),2) mean(C2(4,indTest==4,:,:),2)],[3 4 2 1]);
        resT2(isub,:,:,1,:) = permute([mean(C2(2,indTest==1,:,:),2) mean(C2(1,indTest==2,:,:),2) mean(C2(4,indTest==3,:,:),2) mean(C2(3,indTest==4,:,:),2)],[3 4 2 1]); % within hemifield
        resT3(isub,:,:,1,:) = permute([mean(C2(4,indTest==1,:,:),2) mean(C2(3,indTest==2,:,:),2) mean(C2(2,indTest==3,:,:),2) mean(C2(1,indTest==4,:,:),2)],[3 4 2 1]); % within upper/lower
        resT4(isub,:,:,1,:) = permute([mean(C2(3,indTest==1,:,:),2) mean(C2(4,indTest==2,:,:),2) mean(C2(1,indTest==3,:,:),2) mean(C2(2,indTest==4,:,:),2)],[3 4 2 1]); % diagonal       
        fprintf('\nDone.\n');
    end
    toc
    %%
    close all
    figure
    clims = [-2 2];
    subplot(2,2,1)
    imagesc(tfr.time,tfr.freq,squeeze(mean(resT1(1,:,:,:,:),5)),clims), axis xy, colorbar
    
    subplot(2,2,3)
    imagesc(tfr.time,tfr.freq,squeeze(mean(resT2(1,:,:,:,:),5)),clims), axis xy, colorbar
    
    subplot(2,2,2)
    imagesc(tfr.time,tfr.freq,squeeze(mean(resT3(1,:,:,:,:),5)),clims), axis xy, colorbar
    
    subplot(2,2,4)
    imagesc(tfr.time,tfr.freq,squeeze(mean(resT4(1,:,:,:,:),5)),clims), axis xy, colorbar