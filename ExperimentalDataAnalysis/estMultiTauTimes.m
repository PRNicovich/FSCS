function nTimePoints = estMultiTauTimes(nStart, nEnd, nSub)

nC = abs(nStart-nEnd)+1;
nS = nSub;

j = 1:nC*nS;

tau = 2.^floor((j-1)/nS);
tau = cumsum(tau);

nTimePoints = numel(tau);
