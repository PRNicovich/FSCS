topFolder = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper';

filepath = 'D:\Dropbox\Proposals\FSCS\DataFiles\Figure3\DOPC-50uMchol-6ch-2.pt3';
% filepath = 'D:\Dropbox\Proposals\FSCS\DataFiles\Figure3\DOPC-1250uMchol-6ch-1.pt3';
refLoSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\LUVsB_PSM_Chol_NR12S_2_spectrum.txt';
refLdSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\SLBs_DOPC_Chol_NR12S_3_spectrum.txt';

NCorrSplits = 10; % Number of equally-spaced curves to split each correlation into

nCStart = 7;
nCEnd = 29;
nSub = 10;

averagingRange = 15:25; % Range of points to use as G_inf for std averaging

smoothFactor = 2; % Product of nSub for # points to smooth over
Nkeep = 6;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Lo, Ld ref spectra
LoData = dlmread(refLoSpectrum, '\t', 4, 0); 
Lo = LoData(:,2);
LdData = dlmread(refLdSpectrum, '\t', 4, 0);
Ld = LdData(:,2);

% Start loop here
fprintf(1, 'Loading file %s...\n', filepath);
% Load data file
[chan, AbsTime, macroTime, microTime] = pt3Import(filepath, 'Parallel');

% Filter out any overflow events
AbsTime = AbsTime(chan < 15);
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
plot(LoData(:,1), rS(:,3), 'r');
hold on
plot(LoData(:,1), rS(:,1), 'k:');
plot(LoData(:,1), rS(:,2), 'k--');
plot(LoData(:,1), spec, 'k');
hold off
legend({'Calc', 'Lo', 'Ld', 'Ref'});
xlabel('Wavelength (nm)'); ylabel('Normalized Intensity');

%% 
% Generate weight vectors
% End up dropping last nCut photons to get matrix to evenly divide into
% NCorrSplits columns

nSplit = floor(numel(AbsTime)/NCorrSplits);
nCut = mod(numel(AbsTime), NCorrSplits);
y = AbsTime(1:(end-nCut));

w = zeros(numel(y), 2);
w(:, 1) = f(1, ch(1:(end-nCut))+1);
w(:, 2) = f(2, ch(1:(end-nCut))+1);

% Split data into N channels
C = mat2cell(y(1:(NCorrSplits*nSplit), 1), repmat(nSplit, NCorrSplits, 1), 1);    
y = cell2mat(C');
y = y - repmat(min(y), size(y, 1), 1);
C1 = mat2cell(w(1:(NCorrSplits*nSplit), 1), repmat(nSplit, NCorrSplits, 1), 1);  
C2 = mat2cell(w(1:(NCorrSplits*nSplit), 2), repmat(nSplit, NCorrSplits, 1), 1);  
w = cat(3, cell2mat(C1'), cell2mat(C2'));

%%

% w = ones(sum(chan==1), 1);
% y = AbsTime(chan==1);

tPs = estMultiTauTimes(nCStart, nCEnd, nSub);

corrRes = zeros(tPs, NCorrSplits, 3); % timepoints x sub-samples x [ACF_1 ACF_2 CCF]

% Run correlation once, unweighted, single channel, to check for bad time
% segments.
fprintf(1, 'Running initial correlations...\n');
for k = 1:NCorrSplits
    [corrRes(:,k,1), corrTime] = multiTauWeighted(y(:,k), ones(size(y, 1), 1), nCStart, nCEnd, nSub);
end

% Find which time points best fit to mean curve.  Here going to take Nkeep
% closest ones and ditch the rest
keepPoints = find(timepointsClosestToMean(corrRes(:,:,1), Nkeep));

fprintf(1, 'Calculating weighted correlations...\n');
% Run again with the parts retained, this time with weights applied to get
% spectral channels
for kP = 1:Nkeep
    k = keepPoints(kP);
    [corrRes(:,k,:), corrTime] = multiTauWeighted(y(:,k), squeeze(w(:,k,:)), nCStart, nCEnd, nSub);
end

corrTime = corrTime/1e9;

%%
% Average data.  Should be no discarding points at this stage

[corrAvg, corrStd] = clusterAndAverageCorrelations(corrRes, Nkeep, averagingRange);

%% Smooth resulting Avg curve w/ running average
corrAvgSmooth = zeros(size(corrAvg, 1), size(corrAvg, 2));
for k = 1:size(corrAvg, 2);
    corrAvgSmooth(:,k) = smooth(corrAvg(:,k),nSub*smoothFactor);
end

%% Plot it up

figure(1)
semilogx(corrTime, corrAvg(:,1), 'o', 'MarkerSize', 4);
hold on
semilogx(corrTime, corrAvgSmooth(:,1), 'b');
if size(corrAvg, 2) > 1
semilogx(corrTime, corrAvg(:,2), 'ro', 'MarkerSize', 4);
semilogx(corrTime, corrAvgSmooth(:,2), 'r');
semilogx(corrTime, corrAvg(:,3), 'go', 'MarkerSize', 4);
semilogx(corrTime, corrAvgSmooth(:,3), 'g');
hold off
end

%% Output
fprintf(1, 'Output results to file...\n');
[pth, fname] = fileparts(filepath);

% Save each ACF, CCF as individual tab-delimited files
% Want out: 
% Metadata:
% NCorrSplits 
% nCStart = 7;
% nCEnd = 29;
% nSub = 10;
% Nkeep
% Data:
% Timepoints
% Not doing correlation 
% 
sfname = strcat(fullfile(pth, fname), '_correlationResults.txt');
fID = fopen(sfname, 'w+');
fprintf(fID, '# FSCS Correlation Results\r\n');
fprintf(fID, '# File: %s\r\n', filepath);
fprintf(fID, '# Lo Reference: %s\r\n', refLoSpectrum);
fprintf(fID, '# Ld Reference: %s\r\n', refLdSpectrum);
fprintf(fID, '# Lo Filters: %.3d\t%.3d\t%.3d\t%.3d\t%.3d\t%.3d\r\n', f(1, :));
fprintf(fID, '# Ld Filters: %.3d\t%.3d\t%.3d\t%.3d\t%.3d\t%.3d\r\n', f(2, :));
fprintf(fID, '# NCorrSplits: %.d\r\n', NCorrSplits);
fprintf(fID, '# nCStart / nCEnd: %.d\t%.d\r\n', [nCStart, nCEnd]);
fprintf(fID, '# nSub: %.d\r\n', nSub);
fprintf(fID, '# Nkeep: %.d\r\n', Nkeep);
fprintf(fID, '# Averaging Range: %.d\t:\t%.d\r\n', [averagingRange(1), averagingRange(end)]);
fprintf(fID, '# Smoothing Factor: %.d\r\n', smoothFactor);
fprintf(fID, '# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\r\n');
fprintf(fID, 'Time(s)\tACF_Lo\tACF_LoSmooth\tACF_Ld\tACF_LdSmooth\tCCF\tCCFSmooth\r\n');
for k = 1:numel(corrTime)
   
    fprintf(fID, '%.6d\t%.4d\t%.4d\t%.4d\t%.4d\t%.4d\t%.4d\r\n', ...
        [corrTime(k), corrAvg(k,1), corrAvgSmooth(k,1), corrAvg(k,2), corrAvgSmooth(k,2), corrAvg(k,3), corrAvgSmooth(k,3)]);
    
end
fclose(fID);

% Spectral outputs as separate file
% Reference spectra files
% File here
% Filter coefficients 
% Plot data
sfname = strcat(fullfile(pth, fname), '_spectralResults.txt');
fID = fopen(sfname, 'w+');
fprintf(fID, '# Spectral Filter Results\r\n');
fprintf(fID, '# File: %s\r\n', filepath);
fprintf(fID, '# Lo Reference: %s\r\n', refLoSpectrum);
fprintf(fID, '# Ld Reference: %s\r\n', refLdSpectrum);
fprintf(fID, '# Lo Filters: %.3d\t%.3d\t%.3d\t%.3d\t%.3d\t%.3d\r\n', f(1, :));
fprintf(fID, '# Ld Filters: %.3d\t%.3d\t%.3d\t%.3d\t%.3d\t%.3d\r\n', f(2, :));
fprintf(fID, '# Lo coefficient: %.4d\r\n', c1);
fprintf(fID, '# Ld coefficient: %.4d\r\n', c2);
fprintf(fID, '# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\r\n');
fprintf(fID, 'Wavelength(nm)\tData\tLoScaled\tLdScaled\tSum\r\n');
for k = 1:size(rS, 1);
    fprintf(fID, '%.4d\t%.4d\t%.4d\t%.4d\t%.4d\r\n', [LoData(k,1) spec(k) rS(k,:)]);
end
fclose(fID);

fprintf(1, 'Done!\n');
% End loop



