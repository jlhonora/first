function [m_actual, field_inhomogeneity] = get_object(x_positions,Mspecies)

N = length(x_positions);
m_actual = zeros(Mspecies, N);
x = x_positions;

FOV = x(end) - x(1);


sigma = FOV/4;
center = 0;
amplitude = 20;
field_inhomogeneity = [diff(exp(-(pi/2)*((x-center)/sigma).^2)) 0];
field_inhomogeneity = field_inhomogeneity./(max(abs(field_inhomogeneity))/amplitude);
%field_inhomogeneity = zeros(1,N);
%field_inhomogeneity = ones(1,N).*amplitude;

offsets = [-1.5 2];
offsets = [0 0];

%%%%%%%%%%%% first element %%%%%%%%%%%%%%
% % Part One: rect
%cj = -0.35*FOV; wj = 0.2*FOV; aj = 0.5; ph = 4/FOV; offset = offsets(1);%cj = center, wj = width, aj = amplitude, ph = phase (to generate a complex object)
cj = -0.35*FOV; wj = 0.2*FOV; aj = 2; ph = 4/FOV; offset = offsets(1);%cj = center, wj = width, aj = amplitude, ph = phase (to generate a complex object)
% Mj = aj*wj*sinc(wj*(kx-ph)).*exp(-1i*2*pi*cj*(kx-ph));
% M = M + Mj.*exp(-1i*2*pi*(delta_f(1) + phi_x)*kf);
m_actual(1,:) = m_actual(1,:) + aj*rect((x-cj)/wj).*exp(1i*2*pi*x*ph) + repmat(offset,1,N).*exp(1i*2*pi*x*ph);
% % Part Two: triangle
cj = 0.1*FOV; wj = 0.2*FOV; aj = 2; ph = 0/FOV; offset = offsets(1);
%m_actual(1,:) = m_actual(1,:) + aj*triang((x-cj)/wj).*exp(1i*2*pi*x*ph) + repmat(offset,1,N).*exp(1i*2*pi*x*ph);

%%%%%%%%%%%% second element %%%%%%%%%%%%%%
% % Part One: rect
cj = 0.15*FOV; wj = 0.2*FOV; aj = 2; ph = 6/FOV; offset = offsets(2);%cj = center, wj = width, aj = amplitude, ph = phase
% Mj = aj*wj*sinc(wj*(kx-ph)).*exp(-1i*2*pi*cj*(kx-ph));
% M = M + Mj.*exp(-1i*2*pi*(delta_f(2) + phi_x)*kf);
m_actual(2,:) = m_actual(2,:) + aj*rect((x-cj)/wj).*exp(1i*2*pi*x*ph) + repmat(offset,1,N).*exp(1i*2*pi*x*ph);
% % Part Two: triangle
cj = 0.1*FOV; wj = 0.2*FOV; aj = 2.5; ph = 0/FOV; offset = offsets(2);
%m_actual(2,:) = m_actual(2,:) + aj*triang((x-cj)/wj).*exp(1i*2*pi*x*ph) + repmat(offset,1,N).*exp(1i*2*pi*x*ph);
