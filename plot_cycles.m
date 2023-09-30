function [ time_to_plot, Ir_A_to_plot] = plot_cycles(time_vector, Ir_A, threshold_time, num_zero_crossings)
    % Find indices where time_vector is greater than threshold_time
    indices_gt_threshold = find(time_vector > threshold_time);

    % Extract subset of data based on the identified indices
    Ir_A_sub = Ir_A(indices_gt_threshold);
    time_vector_sub = time_vector(indices_gt_threshold);

    % Find zero crossings by looking where Ir_A changes sign
    zero_crossings = find(diff(sign(Ir_A_sub)));

    % Ensure we have at least num_zero_crossings zero crossings
    if length(zero_crossings) >= num_zero_crossings
        start_idx = zero_crossings(1); % first zero crossing
        end_idx = zero_crossings(num_zero_crossings); % num_zero_crossings-th zero crossing

        time_to_plot = time_vector_sub(start_idx:end_idx);
        Ir_A_to_plot = Ir_A_sub(start_idx:end_idx);

        % % Now plot the data
        % figure;
        % plot(time_to_plot, Ir_A_to_plot);
        % title(['Ir_A over ', num2str(num_zero_crossings / 20), ' cycles for time > ', num2str(threshold_time), 's']);
        % xlabel('Time (s)');
        % ylabel('Ir_A');
        % grid on;
    else
        disp(['Not enough cycles in the data for time > ', num2str(threshold_time), 's']);
    end
end
