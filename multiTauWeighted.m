function [autoFin, autotimeFin] = multiTauWeighted(y, w, NcascStart, NcascEnd, Nsub)
%	[autocorr, autotime] = multiTauWeighted(timestamps, weights, nCStart, nCEnd, nSub)
%   Implement cross-correlation, dealing with non-homogeneous weights

    dt = max(y) - min(y);
    mT = max(y);
%     y = round(y(:),0); % Round off to nanoseconds
%     numshape = size(num, 1);
     
    autotime = zeros((NcascEnd+1)*(Nsub+1),1);
    
    if size(w, 2) == 1;
        % Single channel, autocorrelation only 
        auto = zeros((NcascEnd+1)*(Nsub+1), size(y, 2), 1, 'double');
    elseif size(w,2) == 2;
        % Single channel, autocorrelation in each plus cross-corr 
        auto = zeros((NcascEnd+1)*(Nsub+1), size(y, 2), 3, 'double');
    end
    shift = 0;
    delta = 1;
    

    for j = 0:NcascEnd
        
        % Get unique times and bin weighting factors
        [n, bin] = histc(y, unique(y));
        sing = find(n == 1);
        wn = zeros(numel(n), size(w, 2));
        wn(sing, :) = w(ismember(bin, sing), :);
        
        binSums = unique(n); 
        binSums(binSums < 2) = [];
        for k = 1:numel(binSums)
            multi = find(n == binSums(k));
            Mi    = find(ismember(bin, multi));
            switch binSums(k) % Certainly a more general way to do this....
                case 2 % Presumably only this case will ever be called
                    wn(multi, :) = w(Mi(1:2:end), :) + w(Mi(2:2:end), :); 
                case 3
                    wn(multi, :) = w(Mi(1:3:end), :) + w(Mi(2:3:end), :) + w(Mi(3:3:end), :); 
                case 4
                    wn(multi, :) = w(Mi(1:4:end), :) + w(Mi(2:4:end), :) +...
                        w(Mi(3:4:end), :) + w(Mi(4:4:end), :);  
            end
        end
        
        w = wn;
        y = unique(y);

        for k = 1:Nsub
            
            shift = shift + delta; % Calc new lag time in nanoseconds
            lag = round(shift/delta); % Shift to apply in bins
            
            if j >= NcascStart
                
                % As long as j is above start value (skip v short lags)
  
               % Appy lag and compare lists
               % Can hopefully vectorize this to at least cover whole Nsub range at
               % once
               
               % Can this be vectorized?  At least across Nsub?
               
               i1 = ismember(y,y+lag); % Time stamp, lag match
               i2 = ismember(y+lag,y); % Ch2 -> Ch1
               
               auto((k+(j)*Nsub), :, 1) = sum(w(i1, 1).*w(i2, 1))/delta;
               
               if size(w, 2) > 1
               
                auto((k+(j)*Nsub), :, 2) = sum(w(i1,2).*w(i2,2))/delta;
                auto((k+(j)*Nsub), :, 3) = sum([w(i1,1).*w(i2,2); w(i2,1).*w(i1,2)])/(2*delta);
               end
                
            end
            
            % Add shift to list
            autotime(k+(j)*Nsub) = shift; % Shift time in ns
            
        end % End sub-cascade
        
        % Divide y by 2, round off
        y = ceil(y/2); 
        
        % Double value to next shift
        delta = 2*delta;
        
    end % End cascade
    
    % Correcting for decrease in signal length with each lag increase
    for j = 1:size(auto, 1)
         
        auto(j,:,:) = (auto(j,:,:)*dt)/(dt - autotime(j)); % Eqn 3 + 6
        
    end
    
    auto(:,:,1) = auto(:,:,1)*mT./(sum(w(:,1)).^2); 
    
    if size(w,2) == 2;
        auto(:,:,2) = auto(:,:,2)*mT./(sum(w(:,2)).^2); 
        auto(:,:,3) = auto(:,:,3)*mT./(sum(w(:,1))*sum(w(:,2)));
    end
%     autotime = autotime/1000000; % time/1e6 - output in ms


    % Removes trailing zeros.
    autotimeFin = autotime(sum(auto, 3) ~= 0);
    autoFin = auto(sum(auto, 3) ~= 0, :,  :);
    
    
