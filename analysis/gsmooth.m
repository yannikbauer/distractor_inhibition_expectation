function dat = gsmooth(dat,sd)

% smooth data (dat), where sd the standard deviation of the
% Gaussian function. The units depend on sample rate. I.e., if 1kHz, then the sd
% is in msec

edge = sd*10;
b=normpdf(-1*edge:edge,0,sd);
c = conv(dat,b);
dat = c((edge+1):end-edge);