function [] = PlotGradientFunction(x_positions, optim_vector)
% Plots the objective function gradient.
% Needs to be checked, probably won't work.

global Mspecies

Nsamples = length(optim_vector)/(2*Mspecies + 1);
N = Nsamples;

[xx gg] = GradientFunction(optim_vector); %!!

if(length(gg)~=((2*Mspecies + 1)*Nsamples))
    error('Incorrect length of gradient vector');
end

g_obj = zeros(Mspecies, Nsamples);
g_field = zeros(1,Nsamples);

for ii = 1:Mspecies
    initial_index = ((ii - 1)*2*N);
    g_obj(ii,:) = gg((initial_index + 1):(initial_index + N)) + ...
                1i*gg((initial_index + N + 1):(initial_index + 2*N));
end
g_field = gg(2*Mspecies*N+1:end);


figure(3);
for ii = 1:Mspecies
    subplot(Mspecies, ii, 1)
    plot(x_positions,abs(g_obj(ii,:)));
    title(sprintf('Absolute Value of Object Gradient for species #%d',ii));
end
drawnow

figure(4);
plot(x_positions,g_field);title('Absolute Value of \psi Gradient');
drawnow

set(gcf, 'Color', [1 1 1]);