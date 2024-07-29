clc;
clear;

% Parameters
N = 1000; % number of walkers
n = 100; % number of steps

% Initialize positions
positions = zeros(N, n);

% Random walk simulation
for i = 1:N
    for j = 2:n
        step = randi([0 1]) * 2 - 1; % step can be -1 or 1
        positions(i, j) = positions(i, j-1) + step;
    end
end

% Calculate mean position at final step
mean_position = mean(positions(:, end));

% Plot trajectories
figure;
hold on;
for i = 1:N
    plot(1:n, positions(i, :));
end
hold off;
xlabel('number of steps (n)');
ylabel('position (x)');
title(sprintf('N=%d, n=%d, <x>=%f', N, n, mean_position));

% Plot histogram
figure;
histogram(positions(:, end), 'Normalization', 'count');
xlabel('position');
ylabel('number of particles');
title('100 steps of random walk');
