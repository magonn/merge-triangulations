data = load('res1.txt');

hold on;
plot(data(:, 1), data(:, 3), 'color', rand(1, 3), 'lineWidth', 1.5);
plot(data(:, 1), data(:, 4), 'color', rand(1, 3), 'lineWidth', 1.5);

xlabel('number of points');
ylabel('time, sec');
title('Speed plot')
legend('bridges', 'fictive edges');
hold off;
