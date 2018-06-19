% filepath = 'D:\Dropbox\Proposals\FSCS\DataFiles\Figure3\DOPC-50uMchol-6ch-2.pt3';
filepath = 'D:\Dropbox\Proposals\FSCS\DataFiles\Figure3\DOPC-1250uMchol-6ch-1.pt3';
refLoSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\LUVsB_PSM_Chol_NR12S_2_spectrum.txt';
refLdSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\SLBs_DOPC_Chol_NR12S_3_spectrum.txt';

SamplingRate = 1e6; % Hz
NCorrSplits = 10; % Number of equally-spaced curves to split each correlation into



% Load Lo, Ld ref spectra
LoData = dlmread(refLoSpectrum, '\t', 4, 0); 
Lo = LoData(:,2);
LdData = dlmread(refLdSpectrum, '\t', 4, 0);
Ld = LdData(:,2);

% Load data file
[chan, AbsTime, macroTime, microTime] = pt3Import(filepath, 'Parallel');

ch = microTime(chan < 15);

% Generate Filter coefficients and spectral weighting for this spectrum
% References are Lo and Ld defined above
spec = histc(ch, 0:5);
spec = spec/sum(spec(:));

f = ([Lo, Ld]'*((diag(spec)^-1))*[Lo, Ld])^-1*([Lo, Ld]'*(diag(spec)^-1));
% f = ([Lo(:,2), Ld(:,2)]'*((diag(spec)^-1))*[Lo(:,2), Ld(:,2)])^-1*([Lo(:,2), Ld(:,2)]'*(diag(spec)^-1));

% From LabVIEW implementation:
% dFiltersC = ((dPatternsN / diag(dFFileN)) * dPatternsN') \ dPatternsN / diag(dFFileN);  %own calculation
% Matlab translation of above: 
% dFiltersC = (([Lo(:,2), Ld(:,2)]' / diag(spec)) * [Lo(:,2), Ld(:,2)]) \ [Lo(:,2), Ld(:,2)]' / diag(spec);  %own calculation
% Values tested identical for Figure 3C data

% Given filters, generate spectral weighting coefficients and comparison
% sum spectrum
c1 = sum(spec'.*f(1,:));
c2 = sum(spec'.*f(2,:));

rS = zeros(size(Lo, 1), 3);

rS(:,1) = Lo*c1;
rS(:,2) = Ld*c2;
sumVal = sum(rS(:,2));
rS(:,3) = sum(rS,2);

figure(2)
plot(LoData(:,1), spec, 'r');
hold on
plot(LoData(:,1), rS(:,1), 'k:');
plot(LoData(:,1), rS(:,2), 'k--');
plot(LoData(:,1), rS(:,3), 'k');
hold off
legend({'Calc', 'Lo', 'Ld', 'Ref'});
xlabel('Wavelength (nm)'); ylabel('Normalized Intensity');


%% Generate traces w/ photon weights, given filters above

mT = AbsTime(macroTime > 0);
uChan = unique(ch);
bV = min(mT):(1e9/SamplingRate):max(mT); % nanoseconds

% Split each channel into 10 equal time bins
nSplit = floor(numel(bV)/NCorrSplits);

Mb = zeros(NCorrSplits, nSplit, 'uint32');

[MbHere, bins] = histc(mT, bV);

wHere = zeros(numel(MbHere), 2);

% Case where nPhotons in bin == 1
wHere(MbHere == 1, :) = f(:,1 + ch(ismember(bins, find(MbHere == 1))))';

% Case where nPhotons in bin == 2
% Assume that these cases are where two photons arrive consecutively 
% Then can decompose 2 photon case into a pair of one photon cases with
% alternating logical vector
altLog = [false(floor(numel(bins)/2), 1), true(floor(numel(bins)/2), 1)]';
altLog = altLog(:);
if mod(numel(bins), 2) == 1 % If odd, need to add one more 'false' to the end
    altLog = [altLog; 0];
end

wHere(MbHere == 2, :) = f(:,1 + ch(altLog & ismember(bins, find((MbHere == 2)))))' + f(:,1 + ch(~altLog & ismember(bins, find((MbHere == 2)))))';

% Special cases should cover nearly all (> 99.5%) using (faster) vectorized code
% Remainder covered by for loop here
remBins = find(MbHere > 2);
for k = 1:numel(remBins)
    
    wHere(remBins(k), :) = sum(f(:,1 + ch(bins == remBins(k))), 2)';
 
end


C = mat2cell(wHere(1:(NCorrSplits*nSplit), 1), repmat(nSplit, NCorrSplits, 1), 1);    
chan1 = cell2mat(C')';
C = mat2cell(wHere(1:(NCorrSplits*nSplit), 2), repmat(nSplit, NCorrSplits, 1), 1);  
chan2 = cell2mat(C')';

%% Calc autocorrs, cross-corr, cluster and average

% corrlimit = 0.1*size(Mb,2) ;
% [block_cur{j}, block_err{j}, t_block{j}] = blocking_analysis_ACF(chan1, 0.1, 1, corrlimit) ; 

[Tc1, ~, ~] = compute_CCFwWeights(chan1, chan1, 1/SamplingRate, 2);

[Tc2, ~, ~] = compute_CCFwWeights(chan2, chan2, 1/SamplingRate, 2);

[Tcc, t, crossBinned] = compute_CCFwWeights(chan1, chan2, 1/SamplingRate, 2);

%%
NclusterTraces = 10; % Number of time traces to carry forward.  N closest-clustered ones based on distance from mean will be carried forward.

[Tm1, Tsd1] = clusterAndAverageCorrelations(Tc1, NclusterTraces, 15:25);
[Tm2, Tsd2] = clusterAndAverageCorrelations(Tc2, NclusterTraces, 15:25);
[Tmc, Tsdc] = clusterAndAverageCorrelations(Tcc, NclusterTraces, 15:25);

Tm1 = Tm1 + 1;
Tm2 = Tm2 + 1;
Tmc = Tmc + 1;

%% Plot correlation results

figure(4)
clf(4)
semilogx(t, Tm1, 'bo', 'markersize', 4);
hold on
semilogx(t, Tm2, 'ro', 'markersize', 4);
semilogx(t, Tmc, 'go', 'markersize', 4);
hold off
xlabel('Time (s)');
ylabel('Correlation');
legend({'Lo', 'Ld', 'Cross-correlation'});

%% Add previoiusly-determined correlation results to plot

% refFile = 'D:\Dropbox\Proposals\FSCS\refscs_new\Figure 3\correlation-DOPC-50uMChol-6ch-1.txt';
refFile = 'D:\Dropbox\Proposals\FSCS\refscs_new\Figure 3\correlation-DOPC-1250uMchol-6ch-1.txt';

refCorrs = dlmread(refFile, '\t', 1, 0);
figure(4)
hold on
semilogx(refCorrs(:,1), refCorrs(:,2), 'b-');
semilogx(refCorrs(:,1), refCorrs(:,4), 'r-');
semilogx(refCorrs(:,1), refCorrs(:,6), 'g-');
hold off

% %%
% M = fft(chan1);
% N = conj(fft(chan1));
% 
% Q = ifft(M.*N)/(mean(chan1(:))*mean(chan1(:))*numel(chan1(:)));
% 
% 
% %%
% % Generate Lo autocorr, Ld autocorr, Lo + Ld cross-corr given filter
% % weights, channel spectra
% 
% ordVect = 1:size(Lo, 1);
% ordVect = repmat(ordVect, NCorrSplits, 1);
% holdOrd = ordVect(:);
% % circOrd = reshape(ordVect', [], 1);
% 
% Tmean = zeros(1, numel(t)); % Should be a way to get this value without having to run an auto-correlation first, but it works.
% Tcount = 0;
% 
% for k = 1:size(Lo, 1); % Iterate over number of channels
%                        % Each cycle shifts which channels are
%                        % cross-correlated with one another.
%     [Tcorr, t, crossBinned] = compute_CCFwWeights(Mb, circshift(Mb, [NCorrSplits*(k-1) 0]), ...
%         1/SamplingRate, 4, [], f(1,holdOrd), f(1,circshift(holdOrd, NCorrSplits*[0 (k-1)])));
% %     Tmean = Tmean + sum(Tcorr);
%     Tcount = Tcount + numel(holdOrd); 
%     
%     Tcorr = reshape(Tcorr, [size(Lo, 1) NCorrSplits size(Tcorr, 2)]);
%     
%     % Take NclusterTraces traces forward - filters out weird time points
%     % Pass to function to cluster, average, take stddev
%     [TcorrMean, TcorrStdDev] = clusterAndAverageTraces(Tcorr, NclusterTraces);
%     
%     
%       
%     Tmean = Tmean + TcorrMean + 1;
%     
% end
% 
% Tmean = Tmean/Tcount;
% 
% %% 
% % Calculate autocorrelations for all 6 channels
% 
% mT = AbsTime(macroTime > 0);
% uChan = unique(ch);
% bV = min(mT):(1e9/SamplingRate):max(mT); % nanoseconds
% 
% % Split each channel into 10 equal time bins
% nSplit = floor(numel(bV)/NCorrSplits);
% 
% Mb = zeros(numel(uChan)*NCorrSplits, nSplit, 'uint32');
% 
% cnt = 1;
% cntEnd = NCorrSplits;
% 
% for k = 1:numel(uChan);
%     
%     MbHere = histc(mT(ch == uChan(k)), bV);
%     C = mat2cell(MbHere(1:(NCorrSplits*nSplit)), repmat(nSplit, NCorrSplits, 1), 1);    
%     
%     Mb(cnt:cntEnd,:) = cell2mat(C')';
%     
%     cnt = cnt + NCorrSplits;
%     cntEnd = cntEnd + NCorrSplits;
% 
% end
% 
% j = 1;
% corrlimit = 0.1*size(Mb,2) ;
% [block_cur{j}, block_err{j}, t_block{j}] = blocking_analysis_ACF(Mb, 0.1, 1, corrlimit) ; 
% block_min_ind{j} = block_conver_cri41(block_cur{j}, block_err{j}) ;
%     block_min_ind_s{j} = block_min_ind{j} ;
%     for k = 1:numel(block_min_ind_s{j})
%        if block_min_ind_s{j}(k) == 0
%           block_min_ind_s{j}(k) = numel(t_block{j}) ;
%        else
%        end
%     end 
% [Tcorr, t, crossBinned] = compute_ACF(Mb, 1/SamplingRate, 4);
% Tcorr = Tcorr + 1;
% t = t/1e9;
% 
% figure(2)
% semilogx(t, Tcorr);
% 
% %%
% 
% 
% 
% 
% % TmeanLo = Tmean/Tcount + 1;
% 
% % Tmean = zeros(1, numel(t)); % Should be a way to get this value without having to run an auto-correlation first, but it works.
% Tcount = 0;
% for k = 1:size(Lo, 1);
%     [Tcorr, t, crossBinned] = compute_CCFwWeights(Mb(ordVect,:), Mb(circshift(ordVect, [0 NCorrSplits*(k-1)]),:), ...
%         1/SamplingRate, 4, [], f(2,ordVect), f(2,circshift(ordVect, [0 NCorrSplits*(k-1)])));
% %     Tmean = Tmean + sum(Tcorr);
%     Tcount = Tcount + numel(ordVect);
% end
% 
% % TmeanLd = Tmean/Tcount + 1;
% 
% % Tmean = zeros(1, numel(t)); % Should be a way to get this value without having to run an auto-correlation first, but it works.
% Tcount = 0;
% for k = 1:size(Lo, 1);
%     [Tcorr, t, crossBinned] = compute_CCFwWeights(Mb(ordVect,:), Mb(circshift(ordVect, [0 NCorrSplits*(k-1)]),:), ...
%         1/SamplingRate, 4, [], f(1,ordVect), f(2,circshift(ordVect, [0 NCorrSplits*(k-1)])));
% %     Tmean = Tmean + sum(Tcorr);
%     Tcount = Tcount + numel(ordVect);
% end
% 
% % TmeanCross = Tmean/Tcount + 1;
% 
% t = t/1e9;
% % Make plots
% 
% figure(1)
% semilogx(t, [TmeanLo; TmeanLd; TmeanCross]);
% 
% %% 
% % Wohland et al (Biophysical Journal 80(6) 2987–2999) provides a means
% % To average correlation curves and determine standard deviation
% 
% % Can use this to figure out some selection criteria for measured curves,
% % reject those that we don't want to keep.
% % 
% % Code prior generates 10 curves, which are then averaged.  So for each
% % time sub-set, do the cross-correlation analysis.  That yields 10 curves.
% % Then figure out way to reject some of them, average those together, and
% % figure the standard deviation from that.  Then you end up with your curve
% % + standard deviation for each signal.  
% 
% %% Labview call to FLCS64.dll
% 
% % loadlibrary('FLCS64.dll'); % Needs FLCS64.dll, FLCS64.h in same folder, pwd at that folder
% % 
% % fileName = filepath;
% % pdTimes = mT;
% % pdCorr = 1;
% % pdFilter = f;
% % uiNoRoutes = 6;
% % uiNoChannels = 6;
% % uiNoCompLocal = 2;
% % iNoIntervals = 10;
% % pdTime = 0;
% % piNoLagTimesFinal = 180;
% % dAveInt = 0;
% % iNorm = 1;
% % 
% % calllib('FLCS64','correlation', fileName, pdTimes, pdCorr, pdFilter, uiNoRoutes, ...
% %     uiNoChannels, uiNoCompLocal, iNoIntervals, pdTime, piNoLagTimesFinal, dAveInt, iNorm);
% % 
% % ^ Crashes MATLAB!  
% % 
% % % int correlation(char *sFileName, double *pdTimes, double *pdCorr, double *pdFilter, unsigned int uiNoRoutes,...
% % %     unsigned int uiNoChannels, unsigned int uiNoCompLocal, int iNoIntervals,...
% % %     double *pdTime, int *piNoLagTimesFinal, double *dAveInt, int iNorm);
% % 
% % unloadlibrary('FLCS64');
% 
% 
% 

















