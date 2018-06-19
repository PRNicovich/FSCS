% Draw figures for FSCS paper.
% Do after FSCS Analysis Loop, moving folders around

% Need to define colors, markers for the many different samples used here
SumFilters = [52, 73, 94]/255;

LoVesicle7PSM_3Chol = [46, 204, 90]/255;
LdSLB_4DOPC_1Chol = [142, 68, 173]/255;
Ld_Vesicles = [242, 88, 173]/255;

MonophasicVes = [0, 188, 156]/255;
BiphasicVes = [241, 200, 0]/255;

HomogeneousLUVs = [230, 126, 34]/255;
HeterogeneousLUVs = [45, 200, 219]/255;

LiveCells = [41, 128, 185]/255;
Treated_7kc_mbck = [231, 76, 95]/255;
Treated_LatB = [243, 156, 18]/255;

LoMarker = 'o';
LdMarker = 's';
CrossMarker = 'd';
SumMarker = '^';

%%
% Figure 1A
% Spectra of Lo, Ld standard samples
refLoSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\LUVsB_PSM_Chol_NR12S_2_spectrum.txt';
refLdSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\SLBs_DOPC_Chol_NR12S_3_spectrum.txt';
refLdVesicleSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\NR12S-LUVs-DOPC-6ch-3_spectrum.txt';

% Load Lo, Ld ref spectra
LoData = dlmread(refLoSpectrum, '\t', 4, 0); 
Lo = LoData(:,2);
LdData = dlmread(refLdSpectrum, '\t', 4, 0);
Ld = LdData(:,2);
LdVes = dlmread(refLdVesicleSpectrum, '\t', 4, 0);
LdV = LdData(:,2);

figHand = figure(1);
clf(figHand);
set(figHand, 'color', [1 1 1], 'PaperPositionMode','auto');
ax = axes('parent', figHand, 'FontSize', 10);
hold on

% Loop through available data
plot(ax, LoData(:,1), Lo, 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
plot(ax, LdData(:,1), Ld, 'color', LdSLB_4DOPC_1Chol, 'lineWidth', 2);
% plot(ax, LdVes(:,1), LdV, '--', 'color', Ld_Vesicles, 'lineWidth', 2);

hold off
xlabel('Wavelength (nm)', 'Fontsize', 12); ylabel('Normalized Intensity (a.u.)', 'Fontsize', 12);
set(ax, 'xlim', [552.5 673.5], 'box', 'off', 'ylim', [0.08 0.32], 'Fontsize', 12, 'position', [0.1690    0.1629    0.7360    0.7766]);
pnow = get(gcf, 'position');
set(figHand, 'position', [pnow(1:2), 380 400], 'PaperPositionMode', 'auto');
legend({'L_o Vesicles', 'L_d SLB', 'L_d Vesicles'});

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1A.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1A.eps', '-depsc', '-r300');

%%
% Figure 1B
% Lo autocorrelation


dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 1\pt3\all data\PSM-Chol LUV\LUVsB_PSM_Chol_NR12S_2_correlationResults.txt';
   data = dlmread(dataFile, '\t', 14, 0);
   ampPeak = mean(data(30:53,3));
   ampEnd = mean(data(172:185, 3));
   
   ampHalf = ((ampPeak - ampEnd)/2) + ampEnd;
   
   peakGuess = [find(data(:,3) > ampHalf, 1, 'last'), find(data(:,3) < ampHalf, 1, 'first')];
   
   t_half = interp1(data(peakGuess,3), data(peakGuess,1), ampHalf); 
   t_halfAmp = interp1(data(peakGuess,1), data(peakGuess,3), t_half); 
   
   % Fit to standard correlation curve
   [fitLine, fits1C] = fitCorrelation(data(80:end,[1 2]), '3DT');
   fitLine = correlation3DTriplet(fits1C, data(:,1));
   diffTime = tauTothalfTripletIncluded(fits1C(2), fits1C(3), fits1C(4));
%    diffTime = fits1C(2);
   diffAmp = interp1(data(:,1), fitLine, diffTime); 

figHand = figure(2);
clf(2)
set(2, 'color', [1 1 1]);
semilogx(data(:,1), data(:,2), 'o', 'color', LoVesicle7PSM_3Chol, 'markersize', 6);
hold on

semilogx(data(:,1), data(:,3), '-', 'color', DarkenHue(LoVesicle7PSM_3Chol), 'LineWidth', 2);
% semilogx(data(:,1), fitLine, '--', 'color', ChangeSaturation(DarkenHue(LoVesicle7PSM_3Chol, .3), .1), 'LineWidth', 2);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7], 'linewidth', 1.5);
plot([data(1,1) t_half], [t_halfAmp, t_halfAmp], 'k--', 'linewidth', 1);
plot([t_half t_half], [t_halfAmp, .9], 'k--', 'linewidth', 1);
% plot([diffTime diffTime], [diffAmp, .9], 'k--');
hold off
% [legHand, hobj] = legend({'Autocorrelation', 'Smoothed', 'Fit'}, 'fontsize', 12);
% for k = 5
%     hobj(k).MarkerSize = 10; 
%     hobj(k).LineWidth = 1.5;
% end
axis([data(1,1), data(1,end), 0.98 1.85]);
pnow = get(gcf, 'position');
set(gcf, 'position', [pnow(1:2), 550 400], 'PaperPositionMode', 'auto');
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087]);
set(gca, 'Xtick', logspace(-6, 0, 7), 'xlim', [1e-5 1e0]);
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1B.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1B.eps', '-depsc', '-r300');

%%
% Figure 1C
% Ld autocorrelation, vesicles + SLBs
dataFileVes = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 1\pt3\all data\DOPC LUV\NR12S-LUVs-DOPC-6ch-3_correlationResults.txt';
   dataV = dlmread(dataFileVes, '\t', 14, 0);
   ampPeakV = mean(dataV(30:53,3)); % 1e-5 to 5e-5seconds
   ampEndV = mean(dataV(172:185, 3)); % 0.2-0.5 seconds
   
   ampHalfV = ((ampPeakV - ampEndV)/2) + ampEndV;
   
   peakGuessV = [find(dataV(:,3) > ampHalfV, 1, 'last'), find(dataV(:,3) < ampHalfV, 1, 'first')];
   
   t_halfV = interp1(dataV(peakGuessV,3), dataV(peakGuessV,1), ampHalfV); 
   t_halfAmpV = interp1(dataV(peakGuessV,1), dataV(peakGuessV,3), t_halfV); 
   
   % Fit to standard correlation curve
   [fitLineV, fits1CV] = fitCorrelation(dataV(80:end,[1 2]), '2DT');
   fitLineV = correlation2D(fits1CV, dataV(:,1));
%    diffTimeV = tauTothalfTripletIncluded(fits1CV(2), fits1CV(3), fits1CV(4));
    diffTimeV = fits1CV(2);
   diffAmpV = interp1(dataV(:,1), fitLineV, diffTimeV); 
   
 dataFileSLB = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 1\pt3\all data\DOPC SLB\SLBs_DOPC_Chol_NR12S_3_correlationResults.txt';  
   dataSLB = dlmread(dataFileSLB, '\t', 14, 0);
   ampPeakS = mean(dataSLB(30:53,3));
   ampEndS = mean(dataSLB(172:185, 3));
   
   ampHalfS = ((ampPeakS - ampEndS)/2) + ampEndS;
   
   peakGuessS = [find(dataSLB(:,3) > ampHalfS, 1, 'last'), find(dataSLB(:,3) < ampHalfS, 1, 'first')];
   
   t_halfS = interp1(dataSLB(peakGuessS,3), dataSLB(peakGuessS,1), ampHalfS); 
   t_halfAmpS = interp1(dataSLB(peakGuessS,1), dataSLB(peakGuessS,3), t_halfS); 
   
   % Fit to standard correlation curve
   [fitLineS, fits1CS] = fitCorrelation(dataSLB(:,[1 2]), '2DT');
%    fitLineS = correlation2D(fits1CS, dataSLB(:,1));

%    diffTimeS = tauTothalfTripletIncluded(fits1CS(2), fits1CS(3), fits1CS(4));
   
   diffAmpS = interp1(dataSLB(:,1), fitLineS, diffTimeS); 
   
figHand = figure(3);
clf(figHand)
set(figHand, 'color', [1 1 1]);
semilogx(dataV(:,1), dataV(:,2), 'o', 'color', Ld_Vesicles, 'markersize', 6);
hold on

semilogx(dataV(:,1), dataV(:,3), '-', 'color', DarkenHue(Ld_Vesicles), 'LineWidth', 2);
semilogx(dataV(:,1), fitLineV, '--', 'color', ChangeSaturation(DarkenHue(Ld_Vesicles, .3), .1), 'LineWidth', 2);
plot(dataV(:,1), ones(numel(dataV(:,1)), 1), ':', 'color', [0.7 0.7 0.7], 'linewidth', 1.5);
plot([dataV(1,1) t_halfV], [t_halfAmpV, t_halfAmpV], 'k--', 'linewidth', 1);
plot([t_halfV t_halfV], [t_halfAmpV, .9], 'k--', 'linewidth', 1);
% plot([diffTimeV diffTimeV], [diffAmpV, .9], 'k--');
hold off
% [legHand, hobj] = legend({'Autocorrelation', 'Smoothed', 'Fit'}, 'fontsize', 12);
% for k = 5
%     hobj(k).MarkerSize = 10; 
%     hobj(k).LineWidth = 1.5;
% end
axis([data(1,1), data(1,end), 0.98 1.7]);
pnow = get(gcf, 'position');
set(gcf, 'position', [pnow(1:2), 550 400], 'PaperPositionMode', 'auto');
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087]);
set(gca, 'Xtick', logspace(-6, 0, 7), 'xlim', [1e-5 1e0]);
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1C.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1C.eps', '-depsc', '-r300');

figHand = figure(4);
clf(figHand)
set(figHand, 'color', [1 1 1]);
semilogx(dataSLB(:,1), dataSLB(:,2), 'o', 'color', LdSLB_4DOPC_1Chol, 'markersize', 6);
hold on

semilogx(dataSLB(:,1), dataSLB(:,3), '-', 'color', DarkenHue(LdSLB_4DOPC_1Chol), 'LineWidth', 2);
semilogx(dataSLB(:,1), fitLineS, '--', 'color', ChangeSaturation(DarkenHue(LdSLB_4DOPC_1Chol, .3), .1), 'LineWidth', 2);
plot(dataSLB(:,1), ones(numel(dataSLB(:,1)), 1), ':', 'color', [0.7 0.7 0.7], 'linewidth', 1.5);
plot([dataSLB(1,1) t_halfS], [t_halfAmpS, t_halfAmpS], 'k--', 'linewidth', 1);
plot([t_halfS t_halfS], [t_halfAmpS, .9], 'k--', 'linewidth', 1);
% plot([diffTimeS diffTimeS], [diffAmpS, .9], 'k--');
hold off
% [legHand, hobj] = legend({'Autocorrelation', 'Smoothed', 'Fit'}, 'fontsize', 12);
% for k = 5
%     hobj(k).MarkerSize = 10; 
%     hobj(k).LineWidth = 1.5;
% end
axis([data(1,1), data(1,end), 0.997 1.026]);
pnow = get(gcf, 'position');
set(gcf, 'position', [pnow(1:2), 550 400], 'PaperPositionMode', 'auto');
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087]);
set(gca, 'Xtick', logspace(-6, 0, 7), 'ytick', [.99:.01:1.08]);
set(gca, 'Xtick', logspace(-6, 0, 7), 'xlim', [1e-5 1e0]);
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1D.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig1\Fig1D.eps', '-depsc', '-r300');

%% Figure 2A
% Spectral data in monophasic, diphasic vesicles

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 2\pt3\all data\Mono-SLB_DOPC+LUV_PSM\NR12S-SLBs-DOPC-LUVs-PSM-Chol-6ch-1_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

figHand = figure(10);
clf(figHand);
set(figHand, 'color', [1 1 1]);
pnow = get(gcf, 'position');
set(gcf, 'position', [100, 200, 275 700], 'PaperPositionMode', 'auto');
ax1 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.2600    0.7263    0.7046]);
plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', MonophasicVes, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(MonophasicVes));
set(ax1, 'ylim', [0 0.22]);

hold off
legHand = legend({'L_o Reference', 'L_d Reference', 'Data', 'Calculated'});
set(legHand, 'position', [0.4234    0.4083    0.4909    0.1369]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xticklabels', []);


% Second axes
ax2 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.0771    0.7263    0.1600]);

set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))]);
plot(ax2, data(:,1), f(1,:), 'color', MonophasicVes, 'Marker', LoMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold on
plot(ax2, data(:,1), f(2,:), 'color', MonophasicVes, 'Marker', LdMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
plot(ax2, data(:,1), sum(f), 'color', SumFilters, 'Marker', SumMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold off
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'ylim', [-2.5 3.5], 'ytick', [-3:2:3]);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12); ylabel(ax2, 'Filter Coefficient', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2A.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2A.eps', '-depsc', '-r300');

%% Figure 2B
% Spectral data in monophasic, diphasic vesicles

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 2\pt3\all data\Biphasic\LUVs-DOPC-PSM-CHol-NR12S-6ch-1_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

figHand = figure(10);
clf(figHand);
set(figHand, 'color', [1 1 1]);
pnow = get(gcf, 'position');
set(gcf, 'position', [100, 200, 275 700], 'PaperPositionMode', 'auto');
ax1 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.2600    0.7263    0.7046]);
plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', BiphasicVes, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(BiphasicVes));

hold off
legHand = legend({'L_o Reference', 'L_d Reference', 'Data', 'Calculated'});
set(legHand, 'position', [0.4416    0.4040    0.4909    0.1369]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xticklabels', []);
set(ax1, 'ylim', [0 0.21]);

% Second axes
ax2 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.0771    0.7263    0.1600]);

set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))]);
plot(ax2, data(:,1), f(1,:), 'color', BiphasicVes, 'Marker', LoMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold on
plot(ax2, data(:,1), f(2,:), 'color', BiphasicVes, 'Marker', LdMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
plot(ax2, data(:,1), sum(f), 'color', SumFilters, 'Marker', SumMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold off
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'ylim', [-2.5 3.5], 'ytick', [-3:2:3]);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12); ylabel(ax2, 'Filter Coefficient', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2B.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2B.eps', '-depsc', '-r300');

%% Fig 2C 
% Monophasic vesicle correlation

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 2\pt3\all data\Mono-SLB_DOPC+LUV+PSM\NR12S-SLBs-DOPC-LUVs-PSM-Chol-6ch-1_correlationResults.txt';
data = dlmread(dataFile, '\t', 14, 0);

figHand = figure(13);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
s1 = zeros(6, 1);
s1(1) = semilogx(data(:,1), data(:,2), 'marker', LoMarker, 'color', MonophasicVes, 'linestyle', 'none', 'markersize', 5);
hold on
s1(3) = semilogx(data(:,1), data(:,4), 'marker', LdMarker, 'color', ChangeSaturation(MonophasicVes, 0.4), 'linestyle', 'none', 'markersize', 5);
s1(5) = semilogx(data(:,1), data(:,6), 'marker', CrossMarker, 'color', DarkenHue(MonophasicVes, 0.6), 'linestyle', 'none', 'markersize', 5);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7]);
s1(2) = semilogx(data(:,1), data(:,3), 'color', DarkenHue(MonophasicVes), 'linewidth', 1.5);
s1(4) = semilogx(data(:,1), data(:,5), 'color', DarkenHue(ChangeSaturation(MonophasicVes, 0.4)), 'linewidth', 1.5);
s1(6) = semilogx(data(:,1), data(:,7), 'color', DarkenHue(DarkenHue(MonophasicVes, 0.6)), 'linewidth', 1.5);

hold off
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
set(gca, 'ylim', [.4 3.4])
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
[legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
for k = 5:2:9
    hobj(k).MarkerSize = 10; 
    hobj(k).LineWidth = 1.5;
end

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2C.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2C.eps', '-depsc', '-r300');

%% Fig 2D
% Biphasic vesicle correlation

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 2\pt3\all data\Biphasic\LUVs-DOPC-PSM-CHol-NR12S-6ch-1_correlationResults.txt';
data = dlmread(dataFile, '\t', 14, 0);

figHand = figure(14);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
s1 = zeros(6, 1);
s1(1) = semilogx(data(:,1), data(:,2), 'marker', LoMarker, 'color', BiphasicVes, 'linestyle', 'none', 'markersize', 5);
hold on
s1(3) = semilogx(data(:,1), data(:,4), 'marker', LdMarker, 'color', ChangeSaturation(BiphasicVes, 0.4), 'linestyle', 'none', 'markersize', 5);
s1(5) = semilogx(data(:,1), data(:,6), 'marker', CrossMarker, 'color', DarkenHue(BiphasicVes, 0.6), 'linestyle', 'none', 'markersize', 5);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7]);
s1(2) = semilogx(data(:,1), data(:,3), 'color', DarkenHue(BiphasicVes), 'linewidth', 1.5);
s1(4) = semilogx(data(:,1), data(:,5), 'color', DarkenHue(ChangeSaturation(BiphasicVes, 0.4)), 'linewidth', 1.5);
s1(6) = semilogx(data(:,1), data(:,7), 'color', DarkenHue(DarkenHue(BiphasicVes, 0.6)), 'linewidth', 1.5);

hold off
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
set(gca, 'ylim', [0.93 3.5])
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
[legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
for k = 5:2:9
    hobj(k).MarkerSize = 10; 
    hobj(k).LineWidth = 1.5;
end

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2D.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig2\Fig2D.eps', '-depsc', '-r300');
%% Figure 3A
% Homogeneous vs heterogeneous mixtures of vesicles

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 3\pt3\all data\Homogeneous\DOPC-50uMchol-6ch-2_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

figHand = figure(10);
clf(figHand);
set(figHand, 'color', [1 1 1]);
pnow = get(gcf, 'position');
set(gcf, 'position', [100, 200, 275 700], 'PaperPositionMode', 'auto');
ax1 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.2600    0.7263    0.7046]);
plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', HomogeneousLUVs, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(HomogeneousLUVs));

hold off
legHand = legend({'L_o Reference', 'L_d Reference', 'Data', 'Calculated'});
set(legHand, 'position', [0.4234    0.4083    0.4909    0.1369]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xticklabels', []);


% Second axes
ax2 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.0771    0.7263    0.1600]);

set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))]);
plot(ax2, data(:,1), f(1,:), 'color', HomogeneousLUVs, 'Marker', LoMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold on
plot(ax2, data(:,1), f(2,:), 'color', HomogeneousLUVs, 'Marker', LdMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
plot(ax2, data(:,1), sum(f), 'color', SumFilters, 'Marker', SumMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold off
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'ylim', [-2.5 3.5], 'ytick', [-3:2:3]);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12); ylabel(ax2, 'Filter Coefficient', 'fontsize', 12);
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3A.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3A.eps', '-depsc', '-r300');

%% Figure 3B
% Homogeneous vs heterogeneous mixtures of vesicles

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 3\pt3\all data\Heterogeneous\DOPC-1250uMchol-6ch-1_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

figHand = figure(11);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [300, 200, 275 700], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
ax1 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.2600    0.7263    0.7046]);
plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', HeterogeneousLUVs, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(HeterogeneousLUVs));

hold off
legHand = legend({'L_o Reference', 'L_d Reference', 'Data', 'Calculated'});
set(legHand, 'position', [0.4403    0.3983    0.4934    0.1369]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xticklabels', []);


% Second axes
ax2 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.0771    0.7263    0.1600]);

set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))]);
plot(ax2, data(:,1), f(1,:), 'color', HeterogeneousLUVs, 'Marker', LoMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold on
plot(ax2, data(:,1), f(2,:), 'color', HeterogeneousLUVs, 'Marker', LdMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
plot(ax2, data(:,1), sum(f), 'color', SumFilters, 'Marker', SumMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold off
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'ylim', [-2.5 3.5], 'ytick', [-3:2:3]);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12); ylabel(ax2, 'Filter Coefficient', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3B.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3B.eps', '-depsc', '-r300');

%% Fig 3C 
% Homogeneous vesicles correlation

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 3\pt3\all data\Homogeneous\DOPC-50uMchol-6ch-2_correlationResults.txt';
data = dlmread(dataFile, '\t', 14, 0);

figHand = figure(13);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
s1 = zeros(6, 1);
s1(1) = semilogx(data(:,1), data(:,2), 'marker', LoMarker, 'color', HomogeneousLUVs, 'linestyle', 'none', 'markersize', 5);
hold on
s1(3) = semilogx(data(:,1), data(:,4), 'marker', LdMarker, 'color', ChangeSaturation(HomogeneousLUVs, 0.4), 'linestyle', 'none', 'markersize', 5);
s1(5) = semilogx(data(:,1), data(:,6), 'marker', CrossMarker, 'color', DarkenHue(HomogeneousLUVs, 0.6), 'linestyle', 'none', 'markersize', 5);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7]);
s1(2) = semilogx(data(:,1), data(:,3), 'color', DarkenHue(HomogeneousLUVs), 'linewidth', 1.5);
s1(4) = semilogx(data(:,1), data(:,5), 'color', DarkenHue(ChangeSaturation(HomogeneousLUVs, 0.4)), 'linewidth', 1.5);
s1(6) = semilogx(data(:,1), data(:,7), 'color', DarkenHue(DarkenHue(HomogeneousLUVs, 0.6)), 'linewidth', 1.5);

hold off
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-4 1e0]);
set(gca, 'ylim', [.9 15])
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
[legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
for k = 5:2:9
    hobj(k).MarkerSize = 10; 
    hobj(k).LineWidth = 1.5;
end

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3C.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3C.eps', '-depsc', '-r300');

%% Fig 3D
% Heterogeneous vesicles correlation

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 3\pt3\all data\Heterogeneous\DOPC-1250uMchol-6ch-1_correlationResults.txt';
data = dlmread(dataFile, '\t', 14, 0);

figHand = figure(14);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
s1 = zeros(6, 1);
s1(1) = semilogx(data(:,1), data(:,2), 'marker', LoMarker, 'color', HeterogeneousLUVs, 'linestyle', 'none', 'markersize', 5);
hold on
s1(3) = semilogx(data(:,1), data(:,4), 'marker', LdMarker, 'color', ChangeSaturation(HeterogeneousLUVs, 0.4), 'linestyle', 'none', 'markersize', 5);
s1(5) = semilogx(data(:,1), data(:,6), 'marker', CrossMarker, 'color', DarkenHue(HeterogeneousLUVs, 0.6), 'linestyle', 'none', 'markersize', 5);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7]);
s1(2) = semilogx(data(:,1), data(:,3), 'color', DarkenHue(HeterogeneousLUVs), 'linewidth', 1.5);
s1(4) = semilogx(data(:,1), data(:,5), 'color', DarkenHue(ChangeSaturation(HeterogeneousLUVs, 0.4)), 'linewidth', 1.5);
s1(6) = semilogx(data(:,1), data(:,7), 'color', DarkenHue(DarkenHue(HeterogeneousLUVs, 0.6)), 'linewidth', 1.5);

hold off
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-4 1e0]);
set(gca, 'ylim', [0.93 2.42])
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
[legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
for k = 5:2:9
    hobj(k).MarkerSize = 10; 
    hobj(k).LineWidth = 1.5;
end

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3D.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig3\Fig3D.eps', '-depsc', '-r300');

%% Figure 4A
% Live cells and 3 different treatments
% Live cell spectra, filters

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-live-NR12S-6ch-4_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

figHand = figure(11);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [300, 200, 275 700], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
ax1 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.2600    0.7263    0.7046]);
plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', LiveCells, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(LiveCells));
set(ax1, 'ytick', 0:.05:.3, 'ylim', [0 .27]);

hold off
legHand = legend({'L_o Reference', 'L_d Reference', 'Data', 'Calculated'});
set(legHand, 'position', [0.4258    0.7954    0.4934    0.1369]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xticklabels', []);


% Second axes
ax2 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.0771    0.7263    0.1600]);

set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))]);
plot(ax2, data(:,1), f(1,:), 'color', LiveCells, 'Marker', LoMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold on
plot(ax2, data(:,1), f(2,:), 'color', LiveCells, 'Marker', LdMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
plot(ax2, data(:,1), sum(f), 'color', SumFilters, 'Marker', SumMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold off
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'ylim', [-2.5 3.5], 'ytick', [-3:2:3]);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12); ylabel(ax2, 'Filter Coefficient', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4A.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4A.eps', '-depsc', '-r300');

%% Figure 4B
% Live cells and 3 different treatments
% 7KC spectra, filters

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-7KC-liveB-NR12S-6ch-4_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

figHand = figure(11);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [300, 200, 275 700], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
ax1 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.2600    0.7263    0.7046]);
plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', Treated_7kc_mbck, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(Treated_7kc_mbck));

hold off
legHand = legend({'L_o Reference', 'L_d Reference', 'Data', 'Calculated'});
set(legHand, 'position', [0.2731    0.2740    0.4934    0.1369]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xticklabels', []);
set(ax1, 'ylim', [0 0.2], 'ytick', 0:.05:.3);

% Second axes
ax2 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.0771    0.7263    0.1600]);

set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))]);
plot(ax2, data(:,1), f(1,:), 'color', Treated_7kc_mbck, 'Marker', LoMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold on
plot(ax2, data(:,1), f(2,:), 'color', Treated_7kc_mbck, 'Marker', LdMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
plot(ax2, data(:,1), sum(f), 'color', SumFilters, 'Marker', SumMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold off
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'ylim', [-2.5 3.5], 'ytick', [-3:2:3]);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12); ylabel(ax2, 'Filter Coefficient', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4B.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4B.eps', '-depsc', '-r300');

%% Figure 4C
% Live cells and 3 different treatments
% LatB spectra, filters

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-LatB-NR12S-6ch-4_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

figHand = figure(11);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [300, 200, 275 700], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
ax1 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.2600    0.7263    0.7046]);
plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', Treated_LatB, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(Treated_LatB));

hold off
legHand = legend({'L_o Reference', 'L_d Reference', 'Data', 'Calculated'});
set(legHand, 'position', [0.4331    0.2726    0.4934    0.1369]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xticklabels', []);
set(ax1, 'ylim', [0 0.22], 'ytick', 0:.05:.3);

% Second axes
ax2 = axes('parent', figHand, 'fontsize', 12, 'position', [0.2336    0.0771    0.7263    0.1600]);

set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))]);
plot(ax2, data(:,1), f(1,:), 'color', Treated_LatB, 'Marker', LoMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold on
plot(ax2, data(:,1), f(2,:), 'color', Treated_LatB, 'Marker', LdMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
plot(ax2, data(:,1), sum(f), 'color', SumFilters, 'Marker', SumMarker, 'linewidth', 1.5, 'markersize', 8, 'linestyle', '-');
hold off
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'ylim', [-2.5 3.5], 'ytick', [-3:2:3]);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12); ylabel(ax2, 'Filter Coefficient', 'fontsize', 12);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4C.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4C.eps', '-depsc', '-r300');

%% Fig 4 A - C

% Alternative used for smaller figures, separate filter coefficients
figHand = figure(11);
clf(figHand);

pnow = get(figHand, 'position');
set(figHand, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-live-NR12S-6ch-4_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

ax1 = axes('parent', figHand, 'position', [0.14    0.1450    0.26    0.8146], 'fontsize', 12);

plot(ax1, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax1, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax1, data(:,1), data(:,2), 'color', LiveCells, 'linewidth', 2);
plot(ax1, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(LiveCells));
set(ax1, 'ytick', 0:.05:.3, 'ylim', [0 .27]);
ylabel(ax1, 'Normalized Intensity', 'fontsize', 12);
set(ax1, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xtick', 555:50:655);
set(ax1, 'ylim', [0 0.205], 'ytick', 0:.05:.3);
xlabel(ax1, 'Wavelength (nm)', 'fontsize', 12);

legHand = legend({'L_o Ref', 'L_d Ref', 'Data', 'Calc'}, 'fontsize', 8);
set(legHand, 'position', [0.1932    0.1608    0.1618    0.1883]);

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-7KC-liveB-NR12S-6ch-4_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

ax2 = axes('parent', figHand, 'position', [0.43   0.1450    0.26   0.8146], 'fontsize', 12);

plot(ax2, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax2, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax2, data(:,1), data(:,2), 'color', Treated_7kc_mbck, 'linewidth', 2);
plot(ax2, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(Treated_7kc_mbck));
set(ax2, 'ytick', 0:.05:.3, 'ylim', [0 .27]);
% ylabel(ax2, 'Normalized Intensity', 'fontsize', 12);
set(ax2, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xtick', 555:50:655);
set(ax2, 'ylim', [0 0.205], 'ytick', 0:.05:.3, 'yticklabels', []);
xlabel(ax2, 'Wavelength (nm)', 'fontsize', 12);

legHand = legend({'L_o Ref', 'L_d Ref', 'Data', 'Calc'}, 'fontsize', 8);
set(legHand, 'position', [0.4441    0.1558    0.1618    0.1883]);

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-LatB-NR12S-6ch-4_spectralResults.txt';
% Read in filter data on lines 5 and 6
fID = fopen(dataFile, 'r');
C = textscan(fID, '%s','delimiter', '\n');
fclose(fID);

filtText = textscan(C{1}{5}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(1, :) = horzcat(filtText{2:end});
filtText = textscan(C{1}{6}, '%s %f %f %f %f %f %f', 'delimiter', '\t');
filtText = filtText([1 7 2:6]);
filtText([1 2]) = textscan(filtText{1}{1}, '%s %f', 'delimiter', ':');
f(2, :) = horzcat(filtText{2:end});

data = dlmread(dataFile, '\t', 10, 0);

ax3 = axes('parent', figHand, 'position', [0.72   0.1450    0.26    0.8146], 'fontsize', 12);

plot(ax3, data(:,1), data(:,3), 'color', LoVesicle7PSM_3Chol, 'linewidth', 2);
hold on
plot(ax3, data(:,1), data(:,4), 'color', LdSLB_4DOPC_1Chol, 'linewidth', 2);
plot(ax3, data(:,1), data(:,2), 'color', Treated_LatB, 'linewidth', 2);
plot(ax3, data(:,1), data(:,5), '-', 'linewidth', 2, 'color', DarkenHue(Treated_LatB));
set(ax3, 'ytick', 0:.05:.3, 'ylim', [0 .27]);
% ylabel(ax2, 'Normalized Intensity', 'fontsize', 12);
set(ax3, 'fontsize', 12, 'xlim', [min(data(:,1)), max(data(:,1))], 'xtick', 555:50:655);
set(ax3, 'ylim', [0 0.205], 'ytick', 0:.05:.3, 'yticklabels', []);
xlabel(ax3, 'Wavelength (nm)', 'fontsize', 12);

legHand = legend({'L_o Ref', 'L_d Ref', 'Data', 'Calc'}, 'fontsize', 8);
set(legHand, 'position', [0.7768    0.1608    0.1618    0.1883]);

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4A-C.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4A-C.eps', '-depsc', '-r300');

%% Fig 4D
% Live cell correlation

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-live-NR12S-6ch-4_correlationResults.txt';
data = dlmread(dataFile, '\t', 14, 0);

   t_half = zeros(3, 1);
   t_halfAmp = zeros(3,1);
   for m = 1:2
   ampPeak = mean(data(53:63,3+(m-1))); % 5e-5 to 1e-4
   ampEnd = mean(data(185:195, 3+(m-1))); % .5 to 1
   
   ampHalf = ((ampPeak - ampEnd)/2) + ampEnd;
   
   peakGuess = [find(data(:,3+2*(m-1)) > ampHalf, 1, 'last'), find(data(:,3+2*(m-1)) < ampHalf, 1, 'first')];
   


    t_half(m) = interp1(data(peakGuess,3+2*(m-1)), data(peakGuess,1), ampHalf); 
    t_halfAmp(m) = interp1(data(peakGuess,1), data(peakGuess,3+2*(m-1)), t_half(m)); 
    
   end


figHand = figure(14);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
s1 = zeros(6, 1);
s1(1) = semilogx(data(:,1), data(:,2), 'marker', LoMarker, 'color', LiveCells, 'linestyle', 'none', 'markersize', 5);
hold on
s1(3) = semilogx(data(:,1), data(:,4), 'marker', LdMarker, 'color', ChangeSaturation(LiveCells, 0.4), 'linestyle', 'none', 'markersize', 5);
s1(5) = semilogx(data(:,1), data(:,6), 'marker', CrossMarker, 'color', DarkenHue(LiveCells, 0.6), 'linestyle', 'none', 'markersize', 5);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7]);
s1(2) = semilogx(data(:,1), data(:,3), 'color', DarkenHue(LiveCells), 'linewidth', 1.5);
s1(4) = semilogx(data(:,1), data(:,5), 'color', DarkenHue(ChangeSaturation(LiveCells, 0.4)), 'linewidth', 1.5);
s1(6) = semilogx(data(:,1), data(:,7), 'color', DarkenHue(DarkenHue(LiveCells, 0.6)), 'linewidth', 1.5);

plot([data(1,1) t_half(1)], [t_halfAmp(1), t_halfAmp(1)], '--', 'linewidth', 1, 'color', DarkenHue(LiveCells));
plot([t_half(1) t_half(1)], [t_halfAmp(1), .9], '--', 'linewidth', 1, 'color', DarkenHue(LiveCells));
plot([data(1,1) t_half(2)], [t_halfAmp(2), t_halfAmp(2)], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(LiveCells, 0.4)));
plot([t_half(2) t_half(2)], [t_halfAmp(2), .9], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(LiveCells, 0.4)));
% plot([data(1,1) t_half(3)], [t_halfAmp(3), t_halfAmp(3)], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));
% plot([t_half(3) t_half(3)], [t_halfAmp(3), .9], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));

hold off
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
set(gca, 'ylim', [0.9 1.72])
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
[legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
for k = 5:2:9
    hobj(k).MarkerSize = 10; 
    hobj(k).LineWidth = 1.5;
end

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4D.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4D.eps', '-depsc', '-r300');

%% Fig 4E
% Live cell correlation

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-7KC-liveB-NR12S-6ch-4_correlationResults.txt';
data = dlmread(dataFile, '\t', 14, 0);

   t_half = zeros(3, 1);
   t_halfAmp = zeros(3,1);
   for m = 1:2
   ampPeak = mean(data(53:63,3+(m-1))); % 5e-5 to 1e-4
   ampEnd = mean(data(185:195, 3+(m-1))); % .5 to 1
   
   ampHalf = ((ampPeak - ampEnd)/2) + ampEnd;
   
   pkVect = 1:numel(data(:,1));
   pkVect = (pkVect > 63) & (pkVect < 185);
   peakGuess = [find(pkVect' & (data(:,3+2*(m-1)) > ampHalf), 1, 'last'), find(pkVect' & (data(:,3+2*(m-1)) < ampHalf), 1, 'first')];
   


    t_half(m) = interp1(data(peakGuess,3+2*(m-1)), data(peakGuess,1), ampHalf); 
    t_halfAmp(m) = interp1(data(peakGuess,1), data(peakGuess,3+2*(m-1)), t_half(m)); 
    
   end

figHand = figure(15);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
s1 = zeros(6, 1);
s1(1) = semilogx(data(:,1), data(:,2), 'marker', LoMarker, 'color', Treated_7kc_mbck, 'linestyle', 'none', 'markersize', 5);
hold on
s1(3) = semilogx(data(:,1), data(:,4), 'marker', LdMarker, 'color', ChangeSaturation(Treated_7kc_mbck, 0.4), 'linestyle', 'none', 'markersize', 5);
s1(5) = semilogx(data(:,1), data(:,6), 'marker', CrossMarker, 'color', DarkenHue(Treated_7kc_mbck, 0.6), 'linestyle', 'none', 'markersize', 5);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7]);
s1(2) = semilogx(data(:,1), data(:,3), 'color', DarkenHue(Treated_7kc_mbck), 'linewidth', 1.5);
s1(4) = semilogx(data(:,1), data(:,5), 'color', DarkenHue(ChangeSaturation(Treated_7kc_mbck, 0.4)), 'linewidth', 1.5);
s1(6) = semilogx(data(:,1), data(:,7), 'color', DarkenHue(DarkenHue(Treated_7kc_mbck, 0.6)), 'linewidth', 1.5);

plot([data(1,1) t_half(1)], [t_halfAmp(1), t_halfAmp(1)], '--', 'linewidth', 1, 'color', DarkenHue(Treated_7kc_mbck));
plot([t_half(1) t_half(1)], [t_halfAmp(1), .9], '--', 'linewidth', 1, 'color', DarkenHue(Treated_7kc_mbck));
plot([data(1,1) t_half(2)], [t_halfAmp(2), t_halfAmp(2)], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_7kc_mbck, 0.4)));
plot([t_half(2) t_half(2)], [t_halfAmp(2), .9], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_7kc_mbck, 0.4)));

hold off
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
set(gca, 'ylim', [0.98 1.112])
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
[legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
for k = 5:2:9
    hobj(k).MarkerSize = 10; 
    hobj(k).LineWidth = 1.5;
end

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4E.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4E.eps', '-depsc', '-r300');

%% Fig 4E
% Live cell correlation

dataFile = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\COS7-LatB-NR12S-6ch-4_correlationResults.txt';
data = dlmread(dataFile, '\t', 14, 0);

   t_half = zeros(3, 1);
   t_halfAmp = zeros(3,1);
   for m = 1:3
   ampPeak = mean(data(53:63,3+(m-1))); % 5e-5 to 1e-4
   ampEnd = mean(data(185:195, 3+(m-1))); % .5 to 1
   
   ampHalf = ((ampPeak - ampEnd)/2) + ampEnd;
   
   peakGuess = [find(data(:,3+2*(m-1)) > ampHalf, 1, 'last'), find(data(:,3+2*(m-1)) < ampHalf, 1, 'first')];
   


    t_half(m) = interp1(data(peakGuess,3+2*(m-1)), data(peakGuess,1), ampHalf); 
    t_halfAmp(m) = interp1(data(peakGuess,1), data(peakGuess,3+2*(m-1)), t_half(m)); 
    
   end

figHand = figure(16);
clf(figHand);
pnow = get(gcf, 'position');
set(gcf, 'position', [400, 200, 550, 400], 'PaperPositionMode', 'auto');
set(figHand, 'color', [1 1 1]);
s1 = zeros(6, 1);
s1(1) = semilogx(data(:,1), data(:,2), 'marker', LoMarker, 'color', Treated_LatB, 'linestyle', 'none', 'markersize', 5);
hold on
s1(3) = semilogx(data(:,1), data(:,4), 'marker', LdMarker, 'color', ChangeSaturation(Treated_LatB, 0.4), 'linestyle', 'none', 'markersize', 5);
s1(5) = semilogx(data(:,1), data(:,6), 'marker', CrossMarker, 'color', DarkenHue(Treated_LatB, 0.6), 'linestyle', 'none', 'markersize', 5);
plot(data(:,1), ones(numel(data(:,1)), 1), ':', 'color', [0.7 0.7 0.7]);
s1(2) = semilogx(data(:,1), data(:,3), 'color', DarkenHue(Treated_LatB), 'linewidth', 1.5);
s1(4) = semilogx(data(:,1), data(:,5), 'color', DarkenHue(ChangeSaturation(Treated_LatB, 0.4)), 'linewidth', 1.5);
s1(6) = semilogx(data(:,1), data(:,7), 'color', DarkenHue(DarkenHue(Treated_LatB, 0.6)), 'linewidth', 1.5);

plot([data(1,1) t_half(1)], [t_halfAmp(1), t_halfAmp(1)], '--', 'linewidth', 1, 'color', DarkenHue(Treated_7kc_mbck));
plot([t_half(1) t_half(1)], [t_halfAmp(1), .9], '--', 'linewidth', 1, 'color', DarkenHue(Treated_7kc_mbck));
plot([data(1,1) t_half(2)], [t_halfAmp(2), t_halfAmp(2)], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_7kc_mbck, 0.4)));
plot([t_half(2) t_half(2)], [t_halfAmp(2), .9], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_7kc_mbck, 0.4)));

hold off
set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
set(gca, 'ylim', [0.98 1.5])
xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
[legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
for k = 5:2:9
    hobj(k).MarkerSize = 10; 
    hobj(k).LineWidth = 1.5;
end

print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4F.png', '-dpng', '-r300');
print(figHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4F.eps', '-depsc', '-r300');