clc;
clear;

% Number of simulation steps
n = 10000;

% Initialize the count of points within the quarter circle
inside_circle = 0;

% Generate random points and count how many fall inside the quarter circle
% Plotting the points
figure;
hold on;
theta = linspace(0, pi/2, 100);
plot(cos(theta), sin(theta), 'r', 'LineWidth', 2); % quarter circle
for i = 1:n
    x = rand;
    y = rand;
    if x^2 + y^2 <= 1
        inside_circle = inside_circle + 1;
        plot(x, y, 'g.');
    else
        plot(x, y, 'b.');
    end
end
title('Monte Carlo Method for Estimating \pi');
xlabel('x');
ylabel('y');
legend('Quarter Circle', 'Inside Circle', 'Outside Circle');
axis equal;
hold off;

% Estimate of pi
pi_estimate = (inside_circle / n) * 4;

% Display results
fprintf('Total simulation steps = %d\n', n);
fprintf('Estimated pi value = %.5f\n', pi_estimate);
