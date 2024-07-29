% Initialize parameters
N_values = [100, 100, 1000];
n_values = [1000, 10000, 1000];

for k = 1:length(N_values)
    N = N_values(k);
    n = n_values(k);
    
    % Initialize positions
    positions = zeros(N, n);
    
    % Random walk simulation
    for i = 1:N
        for j = 2:n
            step = randi([0 1]) * 2 - 1; % step can be -1 or 1
            positions(i, j) = positions(i, j-1) + step;
        end
    end
    
    % Calculate <x^2> at each step
    x_squared_mean = mean(positions.^2, 1);
    
    % Plot <x^2> against number of steps
    figure;
    plot(1:n, x_squared_mean, 'LineWidth', 2);
    xlabel('number of steps (n)');
    ylabel('\langle x^2 \rangle');
    title(sprintf('N=%d walkers, n=%d steps, \\langle x^2 \\rangle at final step = %.2f', N, n, x_squared_mean(end)));
end
