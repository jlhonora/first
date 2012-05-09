load('Optim95phantomnonoise');
col_range = [35:93];
m_ideal = abs(m_ideal);
m_total = abs(m_total);
line = 65;
figure(1), clf
plot(col_range, m_ideal(line,col_range,1), 'r--', col_range, m_total(line,col_range,1), 'b', 'LineWidth', 2.5);
%line([1,1],[1,128]);
xlim([col_range(1) col_range(end)])
ylim([0,6]);
%axis off
title('Water profile');
set(gcf, 'Color', [1 1 1]);
h_legend = legend('IDEAL','FIRST')
set(h_legend, 'FontSize', 20)

figure(2), clf
plot(col_range, m_ideal(line,col_range,2)*0.55, 'r--', col_range, m_total(line,col_range,2), 'b', 'LineWidth', 2.5);
axis off
title('Fat profile');
set(gcf, 'Color', [1 1 1]);
xlim([col_range(1) col_range(end)])
ylim([0, 6]);
axis on
h_legend = legend('IDEAL','FIRST')
set(h_legend, 'FontSize', 20)