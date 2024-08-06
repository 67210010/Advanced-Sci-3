% Parameters for SSA
nrun = 100;   % Number of runs (simulations)
nstep = 10000; % Number of steps (time steps in each simulation)
A0 = 20;     % Initial number of molecules A
dt = 0.005;  % Time step size
k = 0.1;     % Reaction rate constant
p = k * dt;  % Probability of reaction occurring in each time step

% Initialize array to store the number of molecules for each run and step
A = zeros(nrun, nstep + 1); % Note: nstep + 1 to accommodate j+1 indexing
A(:, 1) = A0;  % Set initial number of molecules
t = zeros(1, nstep + 1); % Time vector initialization

% Gillespie algorithm function
function [times, A_values] = gillespie_algorithm(A_0, k, max_time)
    times = [0];
    A_values = [A_0];
    A = A_0;
    t = 0;
    while A > 0
        r1 = rand();
        r2 = rand();
        tau = (1 / (A * k)) * log(1 / r1);
        t = t + tau;
        if t > max_time
            break;
        end
        A = A - 1;
        times = [times t];
        A_values = [A_values A];
    end
end

% Parameters for Gillespie
max_time = 50;  % Maximum simulation time
num_trajectories = 10;  % Number of trajectories

% Generate Gillespie data
gillespie_data = cell(num_trajectories, 1);
for i = 1:num_trajectories
    [times, A_values] = gillespie_algorithm(A0, k, max_time);
    gillespie_data{i} = {times, A_values};
end

% Create a common time vector for mean calculation of Gillespie
common_times = linspace(0, max_time, 100);
mean_A_values_gillespie = zeros(size(common_times));

% Interpolate and calculate the mean for Gillespie
for i = 1:numel(common_times)
    A_values_at_time = zeros(num_trajectories, 1);
    for j = 1:num_trajectories
        times = gillespie_data{j}{1};
        A_values = gillespie_data{j}{2};
        A_values_at_time(j) = interp1(times, A_values, common_times(i), 'previous', 'extrap');
    end
    mean_A_values_gillespie(i) = mean(A_values_at_time);
end

% Run SSA simulations
figure;
hold on;
xlabel('Time (sec)');
ylabel('Number of Molecules A');
title('Comparison of Gillespie and SSA');
hMeanSSA = plot(t, zeros(1, nstep + 1), 'k--', 'LineWidth', 2, 'DisplayName', 'Mean (SSA)');
hMeanGillespie = plot(common_times, mean_A_values_gillespie, 'r--', 'LineWidth', 2, 'DisplayName', 'Mean (Gillespie)');
colors = lines(nrun); % Use MATLAB's built-in colormap function 'lines' for distinct colors

for i = 1:nrun % Loop over each simulation run
    for j = 1:nstep % Loop over each time step
        r = rand; % Generate a random number between 0 and 1
        if A(i, j) > 0 % Check if there are molecules present
            if r < p * A(i, j) % If random number is less than probability p * number of molecules
                A(i, j+1) = A(i, j) - 1; % One molecule reacts, decrease count by 1
            else
                A(i, j+1) = A(i, j); % No reaction, count remains the same
            end
        else
            A(i, j+1) = A(i, j); % If no molecules present, count remains the same
        end
        t(j+1) = t(j) + dt; % Update time vector for the next time step
        
        % Update plot for each run periodically
        if mod(j, 100) == 0 || j == nstep % Update plot every 100 steps or at the end of simulation
            plot(t(1:j+1), A(i, 1:j+1), 'Color', colors(i, :), 'HandleVisibility', 'off'); % Plot individual run in unique color
            A_mean = mean(A(1:i, 1:j+1), 1); % Calculate mean number of molecules for current runs
            set(hMeanSSA, 'XData', t(1:j+1), 'YData', A_mean); % Update mean line plot data
            drawnow; % Force MATLAB to update the figure window
        end
    end
end

% Final mean calculation and plot update for SSA
A_mean = mean(A, 1); % Calculate final mean number of molecules at each time step
set(hMeanSSA, 'XData', t, 'YData', A_mean); % Update mean line plot with final data
legend('Mean (SSA)', 'Mean (Gillespie)', 'Location', 'Northeast'); % Add legend to the plot
hold off; % Release the hold on the current figure
