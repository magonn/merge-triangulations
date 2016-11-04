data = load('res.txt');

hold on;
plot(data(:, 1), data(:, 3), 'color', rand(1, 3), 'lineWidth', 2);

xlabel('number of points');
ylabel('time, sec');
title('Speed plot')
legend('fictive edges');
hold off;
