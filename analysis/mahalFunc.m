function [d,z,p] = mahalFunc(dat,lab,nSims)
% calculate Mahalanobis  distance for dat permuting over the first dimension (trials)
% dat [trials, features, samples]
% e.g. dat = rand(150,20,500); % random data, 150 trials, 20 channels, 500 time points
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
    for n=1:nSamples
        P1 = mean(dat(tempLab==1,:,n),1);
        P2 = mean(dat(tempLab==2,:,n),1);
        sigma = cov(squeeze(dat(tempLab==1|tempLab==2,:,n)));
        sDist(sim,n) = (P1-P2)*pinv(sigma)*(P1-P2)';
    end
end

z = zeros(nSamples,1);
p = zeros(nSamples,1);
for n=1:nSamples
    P1 = mean(dat(lab==1,:,n),1);
    P2 = mean(dat(lab==2,:,n),1);
    sigma = cov(squeeze(dat(lab==1|lab==2,:,n)));
    d = (P1-P2)*pinv(sigma)*(P1-P2)';
    
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