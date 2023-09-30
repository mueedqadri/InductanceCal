% Compute the FFT of the current_fea
% f_fft = (0:N-1)*(fs/N); % Frequency vector for FFT
% 
% current_fea_fft = fft(current_fea); % Compute FFT
% current_fea_fft_magnitude = abs(current_fea_fft); % Magnitude of the FFT
% 
% Plot the FFT magnitude spectrum
% figure;
% plot(f_fft, current_fea_fft_magnitude,'red'); 
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('FFT Magnitude Spectrum (Linear Scale)');
% grid on;
% xlim([0 1000]);


% fs = 1/(4e-6);
% 
% % Define window length and overlap
% window_length = length(current_fea)/8; % Example value, can be adjusted
% overlap = window_length/2; % 50% overlap
% nfft = 2^nextpow2(length(current_fea)); % Number of FFT points
% 
% [Pxx, f] = pwelch(current_fea, window_length, overlap, nfft, fs, 'onesided');
% 
% % Find the power at 60Hz
% P_at_60Hz = interp1(f, Pxx, 60, 'linear', 'extrap');
% 
% % Compute the PSD in dB relative to the power at 60Hz
% PSD_dB = 10 * log10(Pxx / P_at_60Hz);
% 
% % Plot the relative PSD in dB
% figure;
% plot(f, PSD_dB);
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz relative to 60Hz)');
% title('PSD with 60Hz normalized to 0 dB');
% % grid on;
% xlim([0 1000]);