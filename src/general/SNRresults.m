clc
load('Optim48brainnoiseSNR10'); load('BRAINWATERRECT');
%load('Optim25'); load('BRAINWATERRECT');
[SNR_W_I, mu, stdev, rectSignal, rectNoise] = SNRcalculator(m_ideal(:,:,1), rectSignal, rectNoise);
[SNR_F_I, mu, stdev, rectSignal, rectNoise] = SNRcalculator(m_ideal(:,:,2), rectSignal, rectNoise);
[SNR_WF_I, mu, stdev, rectSignal, rectNoise] = SNRcalculator(m_ideal(:,:,1) + m_ideal(:,:,2), rectSignal, rectNoise);

disp('IDEAL:');
fprintf(sprintf('W: %5.4g, F: %5.4g, W+F: %5.4g\n', SNR_W_I,SNR_F_I,SNR_WF_I)); 

[SNR_W_F, mu, stdev, rectSignal, rectNoise] = SNRcalculator(m_total(:,:,1), rectSignal, rectNoise);
[SNR_F_F, mu, stdev, rectSignal, rectNoise] = SNRcalculator(m_total(:,:,2), rectSignal, rectNoise);
[SNR_WF_F, mu, stdev, rectSignal, rectNoise] = SNRcalculator(m_total(:,:,1) + m_total(:,:,2), rectSignal, rectNoise);

disp('FIRST:');
fprintf(sprintf('W: %5.4g, F: %5.4g, W+F: %5.4g\n', SNR_W_F,SNR_F_F,SNR_WF_F));

results_mat(6,:) = [SNR_W_I,SNR_F_I,SNR_WF_I, SNR_W_F,SNR_F_F,SNR_WF_F];


