% Number of steps
n = 1000;

% Initial position
x = zeros(1, n+1);
y = zeros(1, n+1);

% Create figure
figure;
plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Mark the original position
hold on;
walker = plot(x(1), y(1), 'b', 'LineWidth', 2);
xlabel('x');
ylabel('y');
grid on;

% Random walk with animation
for i = 2:n+1
    direction = randi(4); % 1: up, 2: down, 3: left, 4: right
    switch direction
        case 1
            x(i) = x(i-1);
            y(i) = y(i-1) + 1;
        case 2
            x(i) = x(i-1);
            y(i) = y(i-1) - 1;
        case 3
            x(i) = x(i-1) - 1;
            y(i) = y(i-1);
        case 4
            x(i) = x(i-1) + 1;
            y(i) = y(i-1);
    end
    % Update plot
    set(walker, 'XData', x(1:i), 'YData', y(1:i));
    title(['2D Random Walk; 1 walker and n=', num2str(i-1)]);
    drawnow;
end

hold off;
