nC = 39; % N cascades
nS = 4; % N points per cascade

% Values for nC, nS are educated guesses for what values need to be there
% to span the data.  Better may be to have a 'while' loop so it at least
% gets to cover 1 -> TraceLength nanoseconds

dT = 1e9; % Range of times through data, in nanoseconds

j = 1:nC*nS;

tau = 2.^floor((j-1)/nS);
tau = cumsum(tau);

ExpLength = 2e9; 
k = max(j)+1;
while  ((max(tau) < ExpLength) && (k < 500))
    
    tauAppend = 2.^floor((k-1)*ones(1,nS)/nS);
    tau = [tau max(tau)+cumsum(tauAppend)];
    k = k+1;
    
end

nCs = ceil(k/nC);

%% 
[chan, AbsTime, MacroTime, MicroTime] = pt3Import('D:\MATLAB\FSCS\FCS_point_correlator/topfluorPE_2_1_1_1.pt3');

AbsTime = AbsTime(chan < 15);
MacroTime = MacroTime(chan < 15);
MicroTime = MicroTime(chan < 15);
chan = chan(chan < 15);

%%
w = zeros(numel(AbsTime), 2);
w(chan == 1, 1) = 1;
w(chan == 2, 2) = 1;
y = (AbsTime);

% w = ones(sum(chan==1), 1);
% y = AbsTime(chan==1);

[a, b] = multiTauWeighted(y, w, 5, 29, 5);

b = b/1e9;

figure(1)
semilogx(b, squeeze(a(:,:,1)));
if size(a, 3) > 1
hold on
semilogx(b, squeeze(a(:,:,2)), 'r');
semilogx(b, squeeze(a(:,:,3)), 'g');
hold off
end


% %%
% validPhotons = subChanArr[subChanArr < 3 ]
% % Creates boolean for photon events in either channel.
% num = zeros(size(validPhotons, 1), 2);
% num(:,1) = (np.array([np.array(validPhotons) ==self.ch_present[0]])).astype(np.int32)
% 		if self.numOfCH ==2:
% 			num(:,2) = (np.array([np.array(validPhotons) ==self.ch_present[1]])).astype(np.int32)

