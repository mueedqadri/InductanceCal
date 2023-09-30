% Filename: fitSinusoid.m
function [fitResult, gof] = fitSinusoid(x, y)
    % Perform Fourier analysis to get a better initial guess for 'b'
    fftResult = fft(y);
    fftMagnitude = abs(fftResult(1:length(y)/2));
    [~, peakIndex] = max(fftMagnitude);
    frequencies = linspace(0, 1/(2*(x(2)-x(1))), length(y)/2);
    dominantFrequency = frequencies(peakIndex);

    % Use the dominant frequency as a starting point for the 'b' parameter
    startPoint = [1, 2*pi*dominantFrequency, 0, 0];

    % Create fit options with the 'StartPoint' option
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', startPoint);

    % Create the fit type with the specified options
    fitType = fittype('a*sin(b*x + c) + d', 'independent', 'x', 'dependent', 'y', 'options', fitOptions);

    % Perform the fit
    [fitResult, gof] = fit(x, y, fitType);
end
