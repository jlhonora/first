line = 65;
point = 38;
acq_range = [46:90];
Nacq = length(acq_range);

values_fat_ideal = zeros(1,Nacq);
values_fat_first = zeros(1,Nacq);
values_water_ideal = zeros(1,Nacq);
values_water_first = zeros(1,Nacq);

for ii = acq_range
   str = sprintf('Optim%dphantomSNR10T2_780', ii);
   load(str);
   m_ideal = abs(m_ideal);
   m_total = abs(m_total);
   values_fat_ideal(ii-acq_range(1)+1) = m_ideal(line,point,2);
   values_water_ideal(ii-acq_range(1)+1) = m_ideal(line,point,1);
   values_fat_first(ii-acq_range(1)+1) = m_total(line,point,2);
   values_water_first(ii-acq_range(1)+1) = m_total(line,point,1);
end

mean_fat_ideal = mean(values_fat_ideal);
mean_water_ideal = mean(values_water_ideal);
mean_fat_first = mean(values_fat_first);
mean_water_first = mean(values_water_first);

std_fat_ideal = std(values_fat_ideal);
std_water_ideal = std(values_water_ideal);
std_fat_first = std(values_fat_first);
std_water_first = std(values_water_first);

snr_fat_ideal = mean_fat_ideal/std_fat_ideal
snr_water_ideal = mean_water_ideal/std_water_ideal
snr_fat_first = mean_fat_first/std_fat_first
snr_water_first = mean_water_first/std_water_first

load('Optim95phantomnonoise');
m_ideal = abs(m_ideal);
m_total = abs(m_total);
orig_fat_ideal = m_ideal(line, point, 2);
orig_water_ideal = m_ideal(line, point, 1);
orig_fat_first = m_total(line, point, 2);
orig_water_first = m_total(line, point, 1);

disp('Fat');
disp('        mean     std     SNR');
fprintf(sprintf('First: %3.4f  %3.4f  %3.4f\n', mean_fat_first, std_fat_first, snr_fat_first));
fprintf(sprintf('Ideal: %3.4f  %3.4f  %3.4f\n', mean_fat_ideal, std_fat_ideal, snr_fat_ideal));

disp('Water');
disp('        mean     std     SNR');
fprintf(sprintf('First: %3.4f  %3.4f  %3.4f\n', mean_water_first, std_water_first, snr_water_first));
fprintf(sprintf('Ideal: %3.4f  %3.4f  %3.4f\n', mean_water_ideal, std_water_ideal, snr_water_ideal));

figure(1)
plot(1:Nacq, values_fat_ideal, 'ro', ...
     1:Nacq, values_fat_first, 'bo', ...
     1:Nacq, repmat(orig_fat_ideal,1,Nacq)', 'r', ...
     1:Nacq, repmat(orig_fat_first,1,Nacq)', 'b');
title('Fat');
xlabel('Acq')
ylabel('Signal');
legend('Ideal', 'First', 'Orig ideal', 'Orig first');

figure(2)
plot(1:Nacq, values_water_ideal, 'ro', ...
     1:Nacq, values_water_first, 'bo', ...
     1:Nacq, repmat(orig_water_ideal,1,Nacq)', 'r', ...
     1:Nacq, repmat(orig_water_first,1,Nacq)', 'b');
title('Water');
xlabel('Acq')
ylabel('Signal');
legend('Ideal', 'First', 'Orig ideal', 'Orig first');





