%% EEG channel configuration
% Script to generate EEG channel configuration used in all pre-processing analysis scripts.

chan.eegOrig = {'FP1','FPZ','FP2','AF7','AF3','AF4','AF8','F7','F5','F3','F1','FZ','F2',...
    'F4','F6','F8','FT7','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT8','T7','C5',...
    'C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4',...
    'CP6','TP8','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4',...
    'PO8','O2','O1','OZ'}';

chan.eeg = {'Fp1','Fpz','Fp2','AF7','AF3','AF4','AF8','F7','F5','F3','F1','Fz','F2',...
    'F4','F6','F8','FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8','T7','C5',...
    'C3','C1','Cz','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPz','CP2','CP4',...
    'CP6','TP8','P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4',...
    'PO8','O2','O1','Oz'}';

chan.eog = {'HEOG','VEOG'};       
chan.ref = {'M1','M2'};
