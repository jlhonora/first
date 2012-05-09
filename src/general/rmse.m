fo = 0;
rectSignal = [38   63    20    3];
rectNoise = rectSignal;
sim_range = 46:90;
sim_string = 'Optim%dphantomSNR10T2_780';
img_box_rows = rectSignal(2):(rectSignal(2) + rectSignal(4));
img_box_cols = rectSignal(1):(rectSignal(1) + rectSignal(3));
Npoints = length(img_box_rows)*length(img_box_cols);
Nsim = length(sim_range);

for ii = sim_range
    load(sprintf(sim_string, ii));
    standard = m_ideal(:,:,1) + m_ideal(:,:,2);
    fo = fo + mean(mean(standard(img_box_rows, img_box_cols))); 
end
fo = fo/length(sim_range);

rmse_ideal = 0;
rmse_first = 0;

for ii = sim_range
    load(sprintf(sim_string, ii));
    f = m_total(img_box_rows, img_box_cols, 1) + m_total(img_box_rows, img_box_cols, 2);
    rmse_first = rmse_first + sqrt(sum(sum((abs(f-fo)).^2)))/(Npoints*Nsim); 
    f = m_ideal(img_box_rows, img_box_cols, 1) + m_ideal(img_box_rows, img_box_cols, 2);
    rmse_ideal = rmse_ideal + sqrt(sum(sum((abs(f-fo)).^2)))/(Npoints*Nsim);
end

fprintf(sprintf('RMSE First: %.5f\n', rmse_first));
fprintf(sprintf('RMSE Ideal: %.5f\n', rmse_ideal));

fprintf(sprintf('RMSE perc First: %.5f\n', rmse_first/abs(fo)*100));
fprintf(sprintf('RMSE perc Ideal: %.5f\n', rmse_ideal/abs(fo)*100));