function [d,z,p] = euclidFunc(dat,lab,nSims)
% calculate euclidean  distance for dat permuting over the first dimension (trials)
% dat [trials, features, samples, samples, samples, etc...]
% e.g. dat = rand(150,20,30,500); % random data, 150 trials, 20 channels, 30 freqs, 500 time points
% lab is a vector of binary condition values [1 or 2]
% nSim is the requests number of permutations
% output DistEEG is the median-centred distance, zDistEEG is the normalised distance
% written by Mark Stokes, 2014
%%

nSamples = size(dat,3);
sDist = zeros(nSims,nSamples);% variable for saving simulated data
for sim=1:nSims
    fprintf(['Doing ' num2str(sim) '\n'])
    
    tempLab = lab(randperm(size(lab,1)));
    
    P1 = squeeze(mean(dat(tempLab==1,:,:),1));
    P2 = squeeze(mean(dat(tempLab==2,:,:),1));
    
    sDist(sim,:) = sqrt(sum((P1-P2).^2));
end

P1 = squeeze(mean(dat(lab==1,:,:),1));
P2 = squeeze(mean(dat(lab==2,:,:),1));
d = sqrt(sum((P1-P2).^2));


z = zeros(nSamples,1);
for n=1:nSamples
    pVal = sum(sDist(:,n) >= d )./nSims;
    if pVal == 1
        pVal = (nSims - 1)/nSims;
    elseif pVal == 0
        pVal = 1/nSims;
    end
    p(n) = pVal;
    z(n) = norminv(1-pVal,0,1);
end

d = d - median(sDist);