%% Condition Trigger Codes Translation for EEG-experiment
% To be used in conjunction with actual trigger code for each trial from condList file

% TARGETS    
% Target location codes
    % T25  (quadrants going clockwise from Top Right)
    cond.T.p25.R = 102:112; % Target Right, 25 percent
    cond.T.p25.TR = 102:106; % Target Top Right, 25 percent
    cond.T.p25.BR = 108:112; % Target Bottom Right, 25 percent
    cond.T.p25.L = 114:124;
    cond.T.p25.BL = 114:118;
    cond.T.p25.TL = 120:124;
    
    % T75
    cond.T.p75.R = 126:136;
    cond.T.p75.TR = 126:130;
    cond.T.p75.BR = 132:136;
    cond.T.p75.L = 138:148;
    cond.T.p75.BL = 138:142;
    cond.T.p75.TL = 142:148;

    % T100
    cond.T.p100.R = 150:160;
    cond.T.p100.TR = 150:154;
    cond.T.p100.BR = 156:160;
    cond.T.p100.L = 162:172;
    cond.T.p100.BL = 162:166;  
    cond.T.p100.TL = 168:172;         

% DISTRACTORS
% Distracter locations 
% NB: Codes for D location follow different order than in T conditions!
    % D25 % (quadrants going clockwise from Top Right)
    cond.D.p25.R = [174, 180, 186, 188, 192, 194];
    cond.D.p25.TR = [180, 186, 192];
    cond.D.p25.BR = [174, 188, 194];
    cond.D.p25.L = [176, 178, 182, 184, 190, 196];
    cond.D.p25.BL = [176, 182, 196];
    cond.D.p25.TL = [178, 184, 190];
              
    % D75
    cond.D.p75.R = [198, 204, 210, 212, 216, 218];
    cond.D.p75.TR = [204, 210, 216];
    cond.D.p75.BR = [198, 212, 218];
    cond.D.p75.L = [200, 202, 206, 208, 214, 220];
    cond.D.p75.BL = [200, 206, 220];
    cond.D.p75.TL = [202, 208, 214];

    % D100
    cond.D.p100.R = [222, 228, 234, 236, 240, 242];
    cond.D.p100.TR = [228, 234, 240];
    cond.D.p100.BR = [222, 236, 242];
    cond.D.p100.L = [224, 226, 230, 232, 238, 244];
    cond.D.p100.BL = [224, 230, 244];  
    cond.D.p100.TL = [226, 232, 238]; 
    
    