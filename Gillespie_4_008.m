clc;
clear;

% Parameters
k1 = 1e-3;
k2 = 1e-2;
k3 = 1.2;
k4 = 1.0;
A0 = 0;
B0 = 0;
num_runs = 10;
num_steps = 500;

% Gillespie Algorithm function
function [t, A, B] = gillespie_algorithm(k1, k2, k3, k4, A0, B0, num_steps)
    t = zeros(1, num_steps + 1);
    A = zeros(1, num_steps + 1);
    B = zeros(1, num_steps + 1);
    A(1) = A0;
    B(1) = B0;
    
    for step = 1:num_steps
        a1 = k1 * A(step) * (A(step) - 1) / 2;
        a2 = k2 * A(step) * B(step);
        a3 = k3;
        a4 = k4;
        a0 = a1 + a2 + a3 + a4;
        
        r1 = rand();
        r2 = rand();
        tau = -log(r1) / a0;
        t(step + 1) = t(step) + tau;
        
        if r2 * a0 < a1
            A(step + 1) = A(step) - 2;
            B(step + 1) = B(step);
        elseif r2 * a0 < a1 + a2
            A(step + 1) = A(step) - 1;
            B(step + 1) = B(step) - 1;
        elseif r2 * a0 < a1 + a2 + a3
            A(step + 1) = A(step) + 1;
            B(step + 1) = B(step);
        else
            A(step + 1) = A(step);
            B(step + 1) = B(step) + 1;
        end
    end
end

% Running simulations
all_times = cell(num_runs, 1);
all_A = cell(num_runs, 1);
all_B = cell(num_runs, 1);

for run = 1:num_runs
    [t, A, B] = gillespie_algorithm(k1, k2, k3, k4, A0, B0, num_steps);
    all_times{run} = t;
    all_A{run} = A;
    all_B{run} = B;
end

% Averaging the results
t_max = max(cellfun(@max, all_times));
t_avg = linspace(0, t_max, num_steps);

A_avg = zeros(1, num_steps);
B_avg = zeros(1, num_steps);
A_count = zeros(1, num_steps);
B_count = zeros(1, num_steps);

for run = 1:num_runs
    A_interp = interp1(all_times{run}, all_A{run}, t_avg, 'linear');
    B_interp = interp1(all_times{run}, all_B{run}, t_avg, 'linear');
    
    valid_A = ~isnan(A_interp);
    valid_B = ~isnan(B_interp);
    
    A_avg(valid_A) = A_avg(valid_A) + A_interp(valid_A);
    B_avg(valid_B) = B_avg(valid_B) + B_interp(valid_B);
    
    A_count(valid_A) = A_count(valid_A) + 1;
    B_count(valid_B) = B_count(valid_B) + 1;
end

A_avg = A_avg ./ A_count;
B_avg = B_avg ./ B_count;

% Plotting results
figure;

subplot(1, 2, 1);
hold on;
for run = 1:num_runs
    h_run = plot(all_times{run}, all_A{run}, 'Color', [0.5, 0.5, 0.5]);
end
h_mean_A = plot(t_avg, A_avg, 'r', 'LineWidth', 2);
xlabel('Time (sec)');
ylabel('Number of Molecules A');
legend([h_run, h_mean_A], {'Run', 'Mean'}, 'Location', 'northeast');
hold off;

subplot(1, 2, 2);
hold on;
for run = 1:num_runs
    h_run = plot(all_times{run}, all_B{run}, 'Color', [0.5, 0.5, 0.5]);
end
h_mean_B = plot(t_avg, B_avg, 'b', 'LineWidth', 2);
xlabel('Time (sec)');
ylabel('Number of Molecules B');
legend([h_run, h_mean_B], {'Run', 'Mean'}, 'Location', 'northeast');
hold off;