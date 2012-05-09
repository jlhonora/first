function [SNR, mu, stdev, rectSignal, rectNoise] = SNRcalculator(Im, varargin)
% El siguiente comando permite elegir rectangulos de senal y ruido
% a traves de getrect
% [SNR, mu, stdev, rectSignal, rectNoise] = SNRcalculator(Im)
%
% Para pasar como parametro los rectangulos de Senal y de Ruido se puede usar
% [SNR, mu, stdev, rectSignal, rectNoise] = SNRcalculator(Im, rectSignal, rectNoise)

Im = abs(Im);

clf
imshow(Im,[]);

if(length(varargin)==2)
    rectSignal = varargin{1};
    rectNoise = varargin{2};
    hSignal = rectangle('Position', rectSignal, 'EdgeColor', 'b');
    hNoise = rectangle('Position', rectNoise, 'EdgeColor', 'r');    
elseif(isempty(varargin))
    disp('Select signal rectangle...');
    rectSignal = getrect;
    hSignal = rectangle('Position', rectSignal, 'EdgeColor', 'b');
    %set(hSignal, 'Color', [1 0 0]);

    disp('Select noise rectangle...');
    rectNoise = getrect;
    hNoise = rectangle('Position', rectNoise, 'EdgeColor', 'r');
    %set(hNoise, 'Color', [0 0 1]);  
else
    error('Unrecognized number of inputs');
end

SignalValues = Im(ceil(rectSignal(2)):floor(rectSignal(2) + rectSignal(4)), ceil(rectSignal(1)):floor(rectSignal(1) + rectSignal(3)));
mu = mean(SignalValues(:));

NoiseValues = Im(ceil(rectNoise(2)):floor(rectNoise(2) + rectNoise(4)), ceil(rectNoise(1)):floor(rectNoise(1) + rectNoise(3)));
stdev = std(NoiseValues(:));

SNR = mu/stdev;

fprintf(sprintf('The SNR is %f, mean = %f, st. dev. = %f\n', SNR, mu, stdev));
