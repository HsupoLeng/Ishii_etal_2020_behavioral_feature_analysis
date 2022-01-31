y = [0.5, 0.6, 0.8, 0.9, 1, ...
    1.1, 1.3, 1.5, 1.5, 1.7, ...
    2.3, 2.3, 2.4, 2.9, 5.0]; 
x = [1, 1, 1, 1, 1, ....
    2, 2, 2, 2, 2, ...
    3, 3, 3, 3, 3];
group_name = {'group 1', 'group 2', 'group 3'}; 

figure()
boxplot(y, x, 'Symbol', '', 'Labels', group_name);
hold on;
scatter(x, y, 12); 
hold off; 
set(gcf, 'Renderer', 'painters'); 