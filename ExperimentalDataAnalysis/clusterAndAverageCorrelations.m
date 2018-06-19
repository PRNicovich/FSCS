function [meanTrace, stdDevTrace] = clusterAndAverageCorrelations(clustArray, Nkeep, avgRange)

    % Average all time points together
    % Quick and dirty straight mean here
    clustMean = (mean(clustArray, 2));
        
    % Find deviation of each individual trace from the mean
    % Sum of squares absolute deviation.  Never cross-comparing traces so no need to normalize. 
    traceDev = squeeze(sqrt(sum((clustArray - repmat(clustMean, [1, size(clustArray, 2), 1])).^2, 1)));

    % Find the Nkeep traces that are closest to the mean
    [~, sortOrd] = sort(traceDev, 1);
    whichChannel = (var(traceDev, [], 1) == min(var(traceDev, [], 1))); % Only use least-varying channel to do this to avoid basing decision on noise
    keepTrace = (le(sortOrd(:, whichChannel), Nkeep));
     
    % Send the remaining traces for averaging, stddev calc
    meanTrace = AvgTrace(clustArray(:, keepTrace, :));
    stdDevTrace = stdTrace(clustArray(:, keepTrace, :));

    function meanTrace = AvgTrace(cA)

        % Straight mean of remaining traces

        meanTrace = squeeze(mean(cA, 2));

    end


    function stdCorr = stdTrace(cA)

        % Std dev method from Wohland, et al, Biophysical Journal 80(6)
        % 2987–2999, Eqn 20 and 21
        % g(\tau) = \frac{1}{L} \Sigma_{l = 1}^{L} \frac{G_l(\tau) - G_{l,\inf}}{G_l(0) - G_{l,\inf}}
        % \sigma = \sqrt{\frac{1}{L - 1} \Sigma_{l = 1}^{L} ( \frac{G_l(\tau) - G_{l,\inf}}{G_l(0) - G_{l,\inf}} - g(\tau))^2

        G_0 = mean(mean(cA(avgRange, :, :), 1), 2);
        g = (1/size(cA, 2)) * sum( ( cA - 0 )./(repmat(G_0, [size(cA, 1), Nkeep, 1]) - 0), 2);

        stdCorr = squeeze(sqrt( (1/(size(cA, 2) - 1)) * sum(((( cA - 0 )./(repmat(G_0, [size(cA, 1), Nkeep, 1]) - 0)) - repmat(g, [1, Nkeep, 1])).^2, 2)));

    end

end