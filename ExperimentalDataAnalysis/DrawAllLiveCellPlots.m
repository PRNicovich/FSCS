% Draw all of the plots from the live cell case so they can all be compared
% before making sweeping conclusions
% Also pull out t_half, amplitudes for each autocorrelation
% Peak time in cross-correlations if possible


topFolder = '\\149.171.80.222\users\Joanna Kwiatek\FSCS paper\Figure 4\pt3\all data\';
files = cellstr(ls(strcat(topFolder, '*.txt')));
files(vertcat(cellfun(@isempty, (strfind(files, '_correlationResults.txt'))))) = [];

LatBFiles = find(~cellfun(@isempty, strfind(files, 'LatB')));
KCFiles = find(~cellfun(@isempty, strfind(files, '7KC')));
LiveFiles = find(cellfun(@isempty, (strfind(files, 'LatB'))) &cellfun(@isempty, (strfind(files, '7KC'))));

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
% Live cell correlation

thalfLive = zeros(numel(LiveFiles), 3);
thalfAmpLive = zeros(numel(LiveFiles), 3);
peakAmpLive = zeros(numel(LiveFiles), 3);

close all

for k = 1:numel(LiveFiles);

    dataFile = fullfile(topFolder, files(LiveFiles(k,:)));
    data = dlmread(dataFile{1}, '\t', 14, 0);

       t_half = zeros(3, 1);
       t_halfAmp = zeros(3,1);
       t_halfPeak = zeros(3,1);
       for m = 1:2
       ampPeak = mean(data(53:63,3+(m-1))); % 5e-5 to 1e-4
       t_halfPeak(m) = ampPeak;
       ampEnd = mean(data(185:195, 3+(m-1))); % .5 to 1

       ampHalf = ((ampPeak - ampEnd)/2) + ampEnd;

       peakGuess = [find(data(1:end-5,3+2*(m-1)) > ampHalf, 1, 'last'), find(data(1:end-5,3+2*(m-1)) < ampHalf, 1, 'first')]; % Force this to not pick one of the last 5 points 



        t_half(m) = interp1(data(peakGuess,3+2*(m-1)), data(peakGuess,1), ampHalf); 
        t_halfAmp(m) = interp1(data(peakGuess,1), data(peakGuess,3+2*(m-1)), t_half(m)); 

       end
       
       thalfLive(k,:) = t_half';
       thalfAmpLive(k,:) = t_halfAmp';
       peakAmpLive(k,:) = t_halfPeak';


    figHand = figure(100+k);
    set(figHand, 'Tag', 'LiveCell');
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
    plot([t_half(1) t_half(1)], [t_halfAmp(1), 0], '--', 'linewidth', 1, 'color', DarkenHue(LiveCells));
    plot([data(1,1) t_half(2)], [t_halfAmp(2), t_halfAmp(2)], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(LiveCells, 0.4)));
    plot([t_half(2) t_half(2)], [t_halfAmp(2), 0], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(LiveCells, 0.4)));
    % plot([data(1,1) t_half(3)], [t_halfAmp(3), t_halfAmp(3)], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));
    % plot([t_half(3) t_half(3)], [t_halfAmp(3), .9], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));

    hold off
    set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
    yRange = [max(max(data(28:193, 2:end))) min(min(data(28:193, 2:end)))];
    set(gca, 'ylim', [yRange(2)-0.02*abs(diff(yRange)) yRange(1)+0.02*abs(diff(yRange))])
    xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
%     [legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
%     for lN = 5:2:9
%         hobj(lN).MarkerSize = 10; 
%         hobj(lN).LineWidth = 1.5;
%     end

end

% Assemble into single figure, N x 4 wide
figObs = findobj('tag', 'LiveCell');
nRows = ceil(size(figObs,1)/4);
rowVect = reshape(repmat(1:nRows, 4, 1), [], 1);
togFig = figure(11);
clf(togFig);
post = get(togFig, 'position');
set(togFig, 'position', [50, 50, 340*4+25 200*nRows+50*(nRows-1)], 'color', [1 1 1], 'PaperPositionMode', 'auto', 'units', 'pixels');
LRSwitch = 0;
for m = 1:numel(figObs)
    ax = get(figObs(m), 'Children');
    set(ax, 'Parent', togFig, 'units', 'pixels', 'fontsize', 10);
    set(findobj('parent', ax, 'type', 'line'), 'markersize', 3);
    
    if LRSwitch == 0
        axPost = [70, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    elseif LRSwitch == 1
        axPost = [400, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    elseif LRSwitch == 2
        axPost = [730, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    elseif LRSwitch == 3
        axPost = [1060, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    end

    
    set(ax, 'xtick', logspace(-5, 0, 6));
    
    LRSwitch = LRSwitch + 1;
    if LRSwitch > 3
        LRSwitch = 0;
    end
    
    set(ax, 'Position', axPost);
end

close(figObs);

print(togFig, 'D:\Dropbox\Proposals\FSCS\Figures\SI\FigS1LiveCellAll.png', '-dpng', '-r300');
print(togFig, 'D:\Dropbox\Proposals\FSCS\Figures\SI\FigS1LiveCellAll.eps', '-depsc', '-r300');

%%
% 7KC + Live cells

thalf7KC = zeros(numel(KCFiles), 3);
thalfAmp7KC = zeros(numel(KCFiles), 3);
peakAmp7KC = zeros(numel(KCFiles), 3);

close all

for k = 1:numel(KCFiles);

    dataFile = fullfile(topFolder, files(KCFiles(k,:)));
    data = dlmread(dataFile{1}, '\t', 14, 0);

       t_half = zeros(3, 1);
       t_halfAmp = zeros(3,1);
       t_halfPeak = zeros(3,1);
       for m = 1:2
       ampPeak = mean(data(53:63,3+(m-1))); % 5e-5 to 1e-4
       t_halfPeak(m) = ampPeak;
       ampEnd = mean(data(185:195, 3+(m-1))); % .5 to 1

       ampHalf = ((ampPeak - ampEnd)/2) + ampEnd;

       peakGuess = [find(data(1:end-5,3+2*(m-1)) > ampHalf, 1, 'last'), find(data(1:end-5,3+2*(m-1)) < ampHalf, 1, 'first')]; % Force this to not pick one of the last 5 points 



        t_half(m) = interp1(data(peakGuess,3+2*(m-1)), data(peakGuess,1), ampHalf); 
        t_halfAmp(m) = interp1(data(peakGuess,1), data(peakGuess,3+2*(m-1)), t_half(m)); 

       end
       
       thalf7KC(k,:) = t_half';
       thalfAmp7KC(k,:) = t_halfAmp';
       peakAmp7KC(k,:) = t_halfPeak';


    figHand = figure(100+k);
    set(figHand, 'Tag', '7KC');
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
    plot([t_half(1) t_half(1)], [t_halfAmp(1), 0], '--', 'linewidth', 1, 'color', DarkenHue(Treated_7kc_mbck));
    plot([data(1,1) t_half(2)], [t_halfAmp(2), t_halfAmp(2)], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_7kc_mbck, 0.4)));
    plot([t_half(2) t_half(2)], [t_halfAmp(2), 0], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_7kc_mbck, 0.4)));
    % plot([data(1,1) t_half(3)], [t_halfAmp(3), t_halfAmp(3)], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));
    % plot([t_half(3) t_half(3)], [t_halfAmp(3), .9], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));

    hold off
    set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
    yRange = [max(max(data(28:193, 2:end))) min(min(data(28:193, 2:end)))];
    set(gca, 'ylim', [yRange(2)-0.02*abs(diff(yRange)) yRange(1)+0.02*abs(diff(yRange))])
    xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
%     [legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
%     for lN = 5:2:9
%         hobj(lN).MarkerSize = 10; 
%         hobj(lN).LineWidth = 1.5;
%     end

end

% Assemble into single figure, N x 4 wide
figObs = findobj('tag', '7KC');
nRows = ceil(size(figObs,1)/4);
rowVect = reshape(repmat(1:nRows, 4, 1), [], 1);
togFig = figure(11);
clf(togFig);
post = get(togFig, 'position');
set(togFig, 'position', [50, 50, 340*4+25 200*nRows+50*(nRows-1)], 'color', [1 1 1], 'PaperPositionMode', 'auto', 'units', 'pixels');
LRSwitch = 0;
for m = 1:numel(figObs)
    ax = get(figObs(m), 'Children');
    set(ax, 'Parent', togFig, 'units', 'pixels', 'fontsize', 10);
    set(findobj('parent', ax, 'type', 'line'), 'markersize', 3);
    
    if LRSwitch == 0
        axPost = [70, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    elseif LRSwitch == 1
        axPost = [400, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    elseif LRSwitch == 2
        axPost = [730, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    elseif LRSwitch == 3
        axPost = [1060, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 50, 260, 175];
        
    end

    
    set(ax, 'xtick', logspace(-5, 0, 6));
    
    LRSwitch = LRSwitch + 1;
    if LRSwitch > 3
        LRSwitch = 0;
    end
    
    set(ax, 'Position', axPost);
end

close(figObs);
print(togFig, 'D:\Dropbox\Proposals\FSCS\Figures\SI\FigS27KCAll.png', '-dpng', '-r300');
print(togFig, 'D:\Dropbox\Proposals\FSCS\Figures\SI\FigS27KCAll.eps', '-depsc', '-r300');


%%
% LatB + live cells

thalfLatB = zeros(numel(LatBFiles), 3);
thalfAmpLatB = zeros(numel(LatBFiles), 3);
peakAmpLatB = zeros(numel(LatBFiles), 3);
corrcross = zeros(numel(LatBFiles), 1);
corrmax = zeros(numel(LatBFiles), 1);

close all

for k = 1:numel(LatBFiles);

    dataFile = fullfile(topFolder, files(LatBFiles(k,:)));
    data = dlmread(dataFile{1}, '\t', 14, 0);

       t_half = zeros(3, 1);
       t_halfAmp = zeros(3,1);
       t_halfPeak = zeros(3,1);
       for m = 1:2
       ampPeak = mean(data(53:63,3+(m-1))); % 5e-5 to 1e-4
       t_halfPeak(m) = ampPeak;
       ampEnd = mean(data(185:195, 3+(m-1))); % .5 to 1

       ampHalf = ((ampPeak - ampEnd)/2) + ampEnd;

       peakGuess = [find(data(1:end-5,3+2*(m-1)) > ampHalf, 1, 'last'), find(data(1:end-5,3+2*(m-1)) < ampHalf, 1, 'first')]; % Force this to not pick one of the last 5 points 



        t_half(m) = interp1(data(peakGuess,3+2*(m-1)), data(peakGuess,1), ampHalf); 
        t_halfAmp(m) = interp1(data(peakGuess,1), data(peakGuess,3+2*(m-1)), t_half(m)); 

       end
       
       thalfLatB(k,:) = t_half';
       thalfAmpLatB(k,:) = t_halfAmp';
       peakAmpLatB(k,:) = t_halfPeak';
       cc = data(find(data(28:193,7) > 1, 1, 'first')+28,1);
       if isempty(cc)
           corrcross(k) = 0; % Never crosses 1
       else
           corrcross(k) = cc; % Where it crosses 1
       end
       
       [pkHt, cm] = findpeaks(data(28:193,7), 'minpeakdistance', 100);
       cm(pkHt < 1) = [];
       pkHt(pkHt < 1) = [];
       
       if ~isempty(cm)
           cm = cm(pkHt == max(pkHt(:)));
           corrmax(k) = data(cm+28, 1);
           
       else
           corrmax(k) = 0;
       end


    figHand = figure(100+k);
    set(figHand, 'Tag', 'LatB');
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

    plot([data(1,1) t_half(1)], [t_halfAmp(1), t_halfAmp(1)], '--', 'linewidth', 1, 'color', DarkenHue(Treated_LatB));
    plot([t_half(1) t_half(1)], [t_halfAmp(1), 0], '--', 'linewidth', 1, 'color', DarkenHue(Treated_LatB));
    plot([data(1,1) t_half(2)], [t_halfAmp(2), t_halfAmp(2)], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_LatB, 0.4)));
    plot([t_half(2) t_half(2)], [t_halfAmp(2), 0], '--', 'linewidth', 1, 'color', DarkenHue(ChangeSaturation(Treated_LatB, 0.4)));
    % plot([data(1,1) t_half(3)], [t_halfAmp(3), t_halfAmp(3)], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));
    % plot([t_half(3) t_half(3)], [t_halfAmp(3), .9], '--', 'linewidth', 1, 'color', DarkenHue(DarkenHue(LiveCells, 0.6)));

    hold off
    set(gca, 'fontsize', 12, 'position', [0.1300    0.1402    0.7750    0.8087], 'xlim', [1e-5 1e0]);
    yRange = [max(max(data(28:193, 2:end))) min(min(data(28:193, 2:end)))];
    set(gca, 'ylim', [yRange(2)-0.02*abs(diff(yRange)) yRange(1)+0.02*abs(diff(yRange))])
    xlabel('Lag (\tau, sec)', 'interpreter', 'tex', 'fontsize', 12, 'fontname', 'helvetica'); ylabel('Correlation (G(\tau))', 'fontsize', 12);
%     [legHand, hobj] = legend(s1([1 3 5]), {'L_o Auto', 'L_d Auto', 'Cross'}, 'fontsize', 12);
%     for lN = 5:2:9
%         hobj(lN).MarkerSize = 10; 
%         hobj(lN).LineWidth = 1.5;
%     end

end
figObs = findobj('tag', 'LatB');
close(figObs);
%%

% Assemble into single figure, N x 5 wide
figObs = findobj('tag', 'LatB');
nRows = ceil(size(figObs,1)/5);
rowVect = reshape(repmat(1:nRows, 5, 1), [], 1);
togFig = figure(11);
clf(togFig);
post = get(togFig, 'position');
set(togFig, 'position', [50, 50, 340*5+25 210*nRows+30*(nRows-1)], 'color', [1 1 1], 'PaperPositionMode', 'auto', 'units', 'pixels');
LRSwitch = 0;
for m = 1:numel(figObs)
    ax = get(figObs(m), 'Children');
    set(ax, 'Parent', togFig, 'units', 'pixels', 'fontsize', 10);
    set(findobj('parent', ax, 'type', 'line'), 'markersize', 3);
    
    if LRSwitch == 0
        axPost = [70, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 20, 260, 175];
        
    elseif LRSwitch == 1
        axPost = [400, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 20, 260, 175];
        
    elseif LRSwitch == 2
        axPost = [730, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 20, 260, 175];
        
    elseif LRSwitch == 3
        axPost = [1060, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 20, 260, 175];
        
    elseif LRSwitch == 4
        axPost = [1390, 200*nRows+50*(nRows-1) - (rowVect(m))*230 + 20, 260, 175];
        
    end

    
    set(ax, 'xtick', logspace(-5, 0, 6));
    
    LRSwitch = LRSwitch + 1;
    if LRSwitch > 4
        LRSwitch = 0;
    end
    
    set(ax, 'Position', axPost);
end


% print(togFig, 'D:\Dropbox\Proposals\FSCS\Figures\SI\FigS3LatBAll.png', '-dpng', '-r300');
% print(togFig, 'D:\Dropbox\Proposals\FSCS\Figures\SI\FigS3LatBAll.eps', '-depsc', '-r300');


%% Assembled data
% thalftimes
thL = log10(thalfLive);
th7 = log10(thalf7KC);
thB = log10(thalfLatB);

plotData = {thL(:,1), thL(:,2), th7(:,1), th7(:,2), thB(:,1), thB(:,2)};
meanthalf = [mean(thL(:,1:2)); mean(th7(:,1:2)); mean(thB(:,1:2))];
semthalf = [prctile(thL, [.25 97.5]); prctile(th7, [.25 97.5]); prctile(thB, [.25 97.5])];

errorUpper = reshape((semthalf(1:2:end, 1:2) - meanthalf)', [], 1);
errorLower = reshape((semthalf(2:2:end, 1:2) - meanthalf)', [], 1);
meanPlot = reshape(meanthalf', [], 1);

sumHand = figure(21);
clf(sumHand);
ax = axes('parent', sumHand);
set(sumHand, 'color', [1 1 1], 'paperpositionmode', 'auto');


plotSpread(plotData, ...
    'xNames', {'Live L_o', 'Live L_d', '7KC L_o', '7KC L_d', 'LatB L_o', 'LatB L_d'}, ...
    'distributionMarkers', {LoMarker, LdMarker, LoMarker, LdMarker, LoMarker, LdMarker}, ...
    'distributionColors', {LiveCells, LiveCells, Treated_7kc_mbck, Treated_7kc_mbck, Treated_LatB, Treated_LatB});

hold on
tt = errorbar(1:6, meanPlot, errorLower(:), errorUpper(:), 'ko', 'linestyle', 'none');
hold off

ylabel(sprintf('log_{10}(ACF decay time) \n(log_{10}(t_{^1/_2}))'), 'Fontsize', 12);
xlabel('Condition + Domain', 'Fontsize', 12);
set(ax, 'FontSize', 12, 'xlim', [.6 6.4], 'ylim', [-3.5 -.5]);
set(get(ax, 'children'), 'markersize', 6, 'Linewidth', 1.5);
set(sumHand, 'position', [680   558   400   350]);
ax.XTickLabelRotation = 45;


print(sumHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4GAssTimescales.png', '-dpng', '-r300');
print(sumHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4GAssTimescales.eps', '-depsc', '-r300');

%%
% ratioOfTHalftimes
thL = log10(thalfLive);
th7 = log10(thalf7KC);
thB = log10(thalfLatB);


plotData = {thL(:,2)./thL(:,1), th7(:,2)./th7(:,1), thB(:,2)./thB(:,1)};
meanthalf = [mean(thL(:,2)./thL(:,1)); mean(th7(:,2)./th7(:,1)); mean(thB(:,2)./thB(:,1))];
semthalf = [prctile(thL(:,2)./thL(:,1), [.25 97.5]); prctile(th7(:,2)./th7(:,1), [.25 97.5]); prctile(thB(:,2)./thB(:,1), [.25 97.5])];

errorUpper = reshape((semthalf(:, 1) - meanthalf)', [], 1);
errorLower = reshape((semthalf(:, 2) - meanthalf)', [], 1);
meanPlot = reshape(meanthalf', [], 1);

sumHand = figure(22);
clf(sumHand)
ax = axes('parent', sumHand);
set(sumHand, 'color', [1 1 1], 'paperpositionmode', 'auto');


plotSpread(plotData, ...
    'xNames', {'Live', '7KC', 'LatB'}, ...
    'distributionMarkers', {CrossMarker, CrossMarker, CrossMarker}, ...
    'distributionColors', {LiveCells, Treated_7kc_mbck, Treated_LatB});

hold on
errorbar(1:3, meanPlot(:), errorLower(:), errorUpper(:), 'ko', 'linestyle', 'none');
hold off

ylabel(sprintf('log_{10}(Decay time ratio) \n(log_{10}(t_{^1/_2, L_d}) / log_{10}(t_{^1/_2, L_o}))'), 'Fontsize', 12);
xlabel('Condition', 'Fontsize', 12);
set(ax, 'FontSize', 12,  'xlim', [.6 3.4]);
set(get(ax, 'children'), 'markersize', 6, 'Linewidth', 1.5);

set(sumHand, 'position', [680   558   400   300]);

print(sumHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4FRelTimescales.png', '-dpng', '-r300');
print(sumHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4FRelTimescales.eps', '-depsc', '-r300');

%%
% ratioOfamplitudes
plotData = {peakAmpLive(:,2)./peakAmpLive(:,1), peakAmp7KC(:,2)./peakAmp7KC(:,1), peakAmpLatB(:,2)./peakAmpLatB(:,1)};
meanthalf = [mean(plotData{1}), mean(plotData{2}), mean(plotData{3})];
semthalf = [prctile(plotData{1}, [.25 97.5]);prctile(plotData{2}, [.25 97.5]);prctile(plotData{3}, [.25 97.5]);];

errorUpper = reshape((semthalf(:, 1) - meanthalf')', [], 1);
errorLower = reshape((semthalf(:, 2) - meanthalf')', [], 1);
meanPlot = reshape(meanthalf', [], 1);


sumHand = figure(23);
clf(sumHand)
ax = axes('parent', sumHand);
set(sumHand, 'color', [1 1 1], 'paperpositionmode', 'auto');


plotSpread(plotData, ...
    'xNames', {'Live', '7KC', 'LatB'}, ...
    'distributionMarkers', {CrossMarker, CrossMarker, CrossMarker}, ...
    'distributionColors', {LiveCells, Treated_7kc_mbck, Treated_LatB});

hold on
errorbar(1:3, meanPlot(:), errorLower(:), errorUpper(:), 'ko', 'linestyle', 'none');
hold off

ylabel(sprintf('Amplitude Ratio \n(G_{L_d}(0) / G_{L_o}(0))'), 'Fontsize', 12);
xlabel('Condition', 'Fontsize', 12);
set(ax, 'FontSize', 12,  'xlim', [.6 3.4]);
set(get(ax, 'children'), 'markersize', 6, 'Linewidth', 1.5);

set(sumHand, 'position', [680   558   400   300]);

print(sumHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4HRelTimescales.png', '-dpng', '-r300');
print(sumHand, 'D:\Dropbox\Proposals\FSCS\Figures\Fig4\Fig4HRelTimescales.eps', '-depsc', '-r300');

