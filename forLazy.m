data = load('good_res.txt');
new = load('res.txt');

hold on;
plot(data(:, 1), data(:, 3), 'color', rand(1, 3), 'lineWidth', 2);
plot(new(:, 1), new(:, 3), 'color', rand(1, 3), 'lineWidth', 2);

xlabel('number of points');
ylabel('time, sec');
title('Speed plot')
legend('fictive edges', 'New');
hold off;