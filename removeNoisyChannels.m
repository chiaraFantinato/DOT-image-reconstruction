function remCh = removeNoisyChannels(nirsData,dRange,SNRrange,K)
    
    % Compute mean and SNR of each channel
    meanValue = mean(nirsData);
    std_nirsData=median(movstd(nirsData,K)); % movstd(A,K)) % if A is a matrix, movstd works along columns 
    SNRValue = mean(nirsData)./std_nirsData;
    
    % Find channels within dRange and > SNR
    remCh = zeros(size(nirsData,2),1); % column vector
    remCh(meanValue > dRange(1) & meanValue < dRange(2) & SNRValue > SNRrange) = 1; % 1 = good channel
                                                                                    % 0 = bad channel
    
    % Channels should be removed even if only 1 wavelength is bad quality
    tmp = [remCh(1:end/2) remCh(end/2 + 1: end)];
    remCh = zeros(size(nirsData,2),1);
    remCh(sum(tmp,2) == 2) = 1; % Keep only channels that have 1 for both wavelengths
    remCh(end/2+1:end) = remCh(1:end/2); % Copy the decision (0 or 1) to the second wavelength

end