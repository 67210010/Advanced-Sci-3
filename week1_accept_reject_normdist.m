clc;
clear;

% Parameters
n = 1000; % Number of samples
mu = 0; % Mean of the normal distribution
sigma = 1; % Standard deviation of the normal distribution
x = linspace(-5, 5, n); % x-axis range
target = (1/(sigma*sqrt(2*pi))) * exp(-0.5 * ((x-mu)/sigma).^2); % Normal distribution PDF

% Generate random samples
x_rand = -5 + 10 * rand(1, n);
y_rand = rand(1, n) * max(target); % Scale by max value of target

% Acceptance criteria
accepted = y_rand < (1/(sigma*sqrt(2*pi))) * exp(-0.5 * ((x_rand-mu)/sigma).^2);
rejected = ~accepted;

% Plotting
figure;
hold on;
plot(x, target, 'k', 'LineWidth', 1.5); % Plot target function
plot(x_rand(accepted), y_rand(accepted), 'bo'); % Accepted samples
plot(x_rand(rejected), y_rand(rejected), 'rx'); % Rejected samples
xlabel('x');
ylabel('target(x)');
legend('target', 'accept', 'reject');
title('Monte Carlo Sampling with Normal Distribution');
hold off;

% Adding credit
annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', 'Credit: Poramet', 'FontSize', 10, 'HorizontalAlignment', 'right', 'EdgeColor', 'none');

