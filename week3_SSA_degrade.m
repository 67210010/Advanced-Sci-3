clc;
clear;

% Parameters
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

% Prepare the figure for dynamic plotting
figure;
hold on;
xlabel('Time (sec)');
ylabel('Number of Molecules A');
hMean = plot(t, zeros(1, nstep + 1), 'k--', 'LineWidth', 2); % Mean line
title('Fixed time step stochastic simulation');

% Define a colormap for the individual runs
colors = lines(nrun); % Use MATLAB's built-in colormap function 'lines' for distinct colors

% Run simulations
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
            plot(t(1:j+1), A(i, 1:j+1), 'Color', colors(i, :)); % Plot individual run in unique color
            A_mean = mean(A(1:i, 1:j+1), 1); % Calculate mean number of molecules for current runs
            set(hMean, 'XData', t(1:j+1), 'YData', A_mean); % Update mean line plot data
            title(sprintf('Fixed time step stochastic simulation (nrun = %d, nstep = %d)', i, j)); % Update title with current run and step
            drawnow; % Force MATLAB to update the figure window
        end
    end
end

% Final mean calculation and plot update
A_mean = mean(A, 1); % Calculate final mean number of molecules at each time step
set(hMean, 'XData', t, 'YData', A_mean); % Update mean line plot with final data
title(sprintf('Fixed time step stochastic simulation (nrun = %d, nstep = %d)', nrun, nstep)); % Set final title with total runs and steps
legend('Mean of Simulations','Individual Simulations', 'Location', 'Northeast'); % Add legend to the plot
hold off; % Release the hold on the current figure
