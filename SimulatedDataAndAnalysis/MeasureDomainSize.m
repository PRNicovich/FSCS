% Measure generated Ising domain size from saved images using pair
% correlation.
% Domains are approximated as Gaussians and FWHM is reported.

%%% Independent Diffusion
% no domain pattern used
% Figure 2 A+B

%%% Static Domains
% Figure 3B - Lo faster
% Need timescales for bigger, smaller domains
% This is the one that needs to get measured

%%% Motile Domains
% Figure 3D - Med domain fast diffusion

%%% Switching domains +- diffusion
% Figure 4 A+B
% no domains used

%%% Motile domains with diffusion
% Figure 5 - Diffusion slower


%%%%%%%%%%%

filePath = 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCSSimulations\FSCS\Static Domains';
rmax = 500;

fileList = ls(strcat(filePath, '\*.tif'));
fileList(cellfun(@isempty, strfind(cellstr(fileList), '=')), :) = [];

corrLength = zeros(size(fileList, 1), 1);

for k = 1:length(fileList)

    if k == 1
        initGuess = [1 1];
    else
        initGuess = pcFits;
    end
    
    img = imread(fullfile(filePath, fileList(k,:)));
    img = img(1501:3499, 1501:3499); % Crop to middle useful section size

    [G, r, g, dg, mask] = get_autocorr(double(img), ones(size(img, 1), size(img, 2)), rmax);

    % Fit to get correlation length
    % fitOut = lsqcurvefit(@IsingPairCorr, [1 1], r, g);
    [pcFits] = lsqcurvefit(@GaussianPairCorr, [1 1], r(2:end), g(2:end));
    
    corrLength(k) = pcFits(1)*2*sqrt(2*log(2));
    
end

%% Format for processing

assembledCell = cell(size(fileList, 1), 3);

for k = 1:size(fileList, 1)
    
    stringParts = strsplit(fileList(k,:), '_');
    assembledCell{k,1} = stringParts{2}(1:min([9, length(stringParts{2})]));
    
    assembledCell{k,2} = corrLength(k);
    
end

%% Measure peak correlation times for each of these same curves, then add to assembled data

dataFiles = ls(strcat(filePath, '\*.txt'));
timeMaxList = zeros(size(dataFiles, 1)*20, 1);

for k = 1:size(dataFiles, 1)
   
    dataReturn = FSCSSimTXTReader(fullfile(filePath, dataFiles(k,:)));  
    ccorrOnly = dataReturn(:,4:4:end);
    [~, maxCC] = max(ccorrOnly);
    timeMaxCC = dataReturn(maxCC, 1);
    
    timeMaxList(((k-1)*20+1):(k*20)) = timeMaxCC;
    
end

%%
nameCell = cell(size(fileList, 1), 2);

for m = 1:size(dataFiles, 1);
    stringParts = strsplit(dataFiles(m,:), '_');
    for k = 1:20
        nameCell{(m-1)*20 + k,1} = stringParts{2}(1:min([9, length(stringParts{2})]));
        nameCell{(m-1)*20 + k,2} = timeMaxList((m-1)*20 + k, 1);
    end
    
end

[sortNames, Idx] = sort(nameCell(:,1));
[sortAssem, IdxAss] = sort(assembledCell(:,1));

sortCell = {assembledCell{IdxAss,1}; assembledCell{IdxAss, 2}; nameCell{Idx, 2}}';
sortData = [sortCell{:,2}; sortCell{:,3}]';
sortData = [sortData, ones(size(sortData, 1), 1)];
sortData(sortData(:,1) > 10, 3) = 2;
sortData(sortData(:,1) > 35, 3) = 3;

smallMean = mean(sortData(sortData(:,3) == 1, 1:2));
medMean = mean(sortData(sortData(:,3) == 2, 1:2));
bigMean = mean(sortData(sortData(:,3) == 3, 1:2));

smallStd = std(sortData(sortData(:,3) == 1, 1:2));
medStd = std(sortData(sortData(:,3) == 2, 1:2));
bigStd = std(sortData(sortData(:,3) == 3, 1:2));




