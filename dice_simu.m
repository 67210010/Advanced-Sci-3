clc;
clear;

% Initialize the figure
figure;
hold on;

% Theoretical probabilities
theoretical_prob = [1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1] / 36;
values = 2:12;

hand_prob = [1, 2, 8, 11, 17, 13, 12, 13, 12, 8, 3] / 100;

for num_steps = 1:1000
    % Initialize the results array
    results = zeros(num_steps, 1);

    % Simulate rolling two dice
    for i = 1:num_steps
        die1 = randi(6);
        die2 = randi(6);
        results(i) = die1 + die2;
    end

    % Calculate the probabilities
    hist_counts = histcounts(results, 'BinLimits', [1.5, 12.5], 'BinWidth', 1);
    probability = hist_counts / num_steps;

    % Clear the previous plot
    cla;

    % Plot histogram
    histogram(results, 'BinEdges', 1.5:1:12.5, 'Normalization', 'probability', 'FaceAlpha', 0.5);

    % Plot simulated probabilities
    plot(values, probability, 'b', 'LineWidth', 2);

    % Plot theoretical probabilities
    plot(values, theoretical_prob, 'r', 'LineWidth', 2);

    % Plot hand probabilities
    plot(values, hand_prob, 'g', 'LineWidth', 2);

    % Update the title with the current number of steps
    title(['Probability Distribution of Rolling Two Dice - Steps: ', num2str(num_steps)]);
    
    % Update the legend, labels, and grid
    legend('Histogram', 'Simulated', 'Theoretical','Hand rolling for 100 roll');
    xlabel('Sum of Rolling Two Dice');
    ylabel('Probability');
    grid

    % Pause to update the plot
    pause(0.01);
end

hold off;