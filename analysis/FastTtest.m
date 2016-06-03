function tval = FastTtest(x)

mx = squeeze(mean(x));
sx = squeeze(std(x));
tval = mx./(sx./sqrt(size(x,1)));