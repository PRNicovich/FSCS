% Perform FSCS correlation on time-domain data
% Fluor is t x c data, with Fluor(i, j) intensity of channel j at time i
% LoSpectrum and LdSpectrum reference spectra for Lo and Ld phases, of size
% 1 x c
% corrOut is t x 4, with columns [time, LoAutoCorr, LdAutoCorr, CrossCorr]

function corrOut = FSCSintensityCorrelation(spectralFluor, LoSpectrum, LdSpectrum)

    % experimental spectrum
    spec = sum(spectralFluor);
    spec = spec/sum(spec(:));

    % Filter coefficients
    f = ([LoSpectrum, LdSpectrum]'*((diag(spec)^-1))*...
        [LoSpectrum, LdSpectrum])^-1*([LoSpectrum, LdSpectrum]'*(diag(spec)^-1));

    LoTrace = sum(spectralFluor.*repmat(f(1,:), size(spectralFluor, 1), 1), 2);
    LdTrace = sum(spectralFluor.*repmat(f(2,:), size(spectralFluor, 1), 1), 2);

    corrOut = zeros(size(spectralFluor, 1), 4);
    corrOut(:,1) = 0:(size(spectralFluor, 1)-1);

    % Lo autocorr
    corrOut(:,2) = ifft(fft(LoTrace - repmat(mean(LoTrace), size(spectralFluor, 1), 1)).* ...
        conj(fft(LoTrace - repmat(mean(LoTrace), size(spectralFluor, 1), 1))))/...
        (size(spectralFluor, 1).^2) + 1; 

    % Ld autocorr
    corrOut(:,3) = ifft(fft(LdTrace - repmat(mean(LdTrace), size(spectralFluor, 1), 1)).* ...
        conj(fft(LdTrace - repmat(mean(LdTrace), size(spectralFluor, 1), 1))))/...
        (size(spectralFluor, 1).^2) + 1; 

    % Lo + Ld cross-corr
    corrOut(:,4) = ifft(fft(LoTrace - repmat(mean(LoTrace), size(spectralFluor, 1), 1)).* ...
        conj(fft(LdTrace - repmat(mean(LdTrace), size(spectralFluor, 1), 1))))/...
        (size(spectralFluor, 1).^2) + 1; 
    

end