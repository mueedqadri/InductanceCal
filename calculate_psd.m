function [f, PSD_dB, P_at_60Hz] = calculate_psd(time_data, current_data, segment_fraction)
    % Calculate the sampling frequency
    incr = time_data(3) - time_data(2);
    fs = 1 / incr;

    % Calculate the segment length
    N = length(current_data);
    segment_length = segment_fraction * N;

    % Calculate the power spectral density using pwelch
    [Pxx, f] = pwelch(current_data, segment_length, [], [], fs, 'onesided');

    % Find the power at 60Hz
    P_at_60Hz = interp1(f, Pxx, 60, 'linear', 'extrap');

    % Normalize the power spectral density and convert to dB
    PSD_dB = 10 * log10(Pxx / P_at_60Hz);
end
