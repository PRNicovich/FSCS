function keepTrace = timepointsClosestToMean(clustArray, Nkeep)
% Given 2D cluster array, return the Nkeep times closest to the mean

% Average all time points together
% Quick and dirty straight mean here
clustMean = (mean(clustArray, 2));

% Find deviation of each individual trace from the mean
% Sum of squares absolute deviation.  Never cross-comparing traces so no need to normalize.
traceDev = squeeze(sqrt(sum((clustArray - repmat(clustMean, [1, size(clustArray, 2), 1])).^2, 1)));

% Find the Nkeep traces that are closest to the mean
[~, sortOrd] = sort(traceDev, 2);

keepTrace = (le(sortOrd, Nkeep));
