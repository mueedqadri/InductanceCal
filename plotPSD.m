function plotPSD(PSD_dB_fea_h, f, PSD_dB_fea_f, f_fea_f, doZoom, zoomRange, legend1, legend2)
    if nargin < 6
        zoomRange = [620, 800];
    end
    blue = [0 0.4470 0.7410];
    red = [0.8500 0.3250 0.0980];


    % Plotting
    figure;
    plot(f, PSD_dB_fea_h, 'Color', blue, 'LineWidth', 1);
    hold on;
    plot(f_fea_f, PSD_dB_fea_f, 'Color', red, 'LineWidth', 1, 'LineStyle', '-.');
    
    xlabel('Frequency (Hz)', 'FontSize', 15);
    ylabel('PSD (dB)', 'FontSize', 15);
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    grid on;
    xlim([0 1000]);
    lgd = legend(legend1, legend2);
    set(lgd, 'FontSize', 14);

    if doZoom
        % Highlight the zoomed region on the main plot
        indexOfInterest1 = (f > zoomRange(1)) & (f < zoomRange(2));
        indexOfInterest2 = (f_fea_f > zoomRange(1)) & (f_fea_f < zoomRange(2));
        rectangle('Position', [zoomRange(1), min(PSD_dB_fea_f(indexOfInterest2)), zoomRange(2)-zoomRange(1), max(PSD_dB_fea_f(indexOfInterest2)) - min(PSD_dB_fea_f(indexOfInterest2))], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

        axes('position', [.60 .55 .25 .18]);
        box on;
        plot(f_fea_f(indexOfInterest2), PSD_dB_fea_f(indexOfInterest2), 'Color', red, 'LineWidth', 2);
        hold on;
        plot(f(indexOfInterest1), PSD_dB_fea_h(indexOfInterest1), 'Color', blue, 'LineWidth', 2);
        axis tight;
        grid on;
        hold off;
        box off;
    end
    
    hold off;
end