function [m,m_opt,field_inhomogeneity,phi_x_opt] = SpeciesSeparationN(params)
% Chemical Species Separation using k-space formulation
% Used for compiling with mcc
%
% WARNING : This script is outdated and only includes a simulation.
%           To work properly, pleasy copy the script SpeciesSeparation2DMAIN
%           
% As a suggestion, make a detail of all parameters in a txt file and then
% read them in this script. Pass the name of the parameters' text file as
% the script's argument.
%
% Use this line: 
%
% mcc -d ./binaries_directory -R 'ErrorLog.txt' -m -o SpeciesSeparationN
%
% Then go to the binaries directory from command line and call
%
% binaries_directory> SpeciesSeparationN parameters.txt
%
%
% JLH
% Nov. 2010

%clear all
close all
clc

disp('Chemical Species Separation using k-space formulation');
disp('JLH');
disp('Dec. 2010');
error('Scrpit outdated, pleasy modify');

global kx kf Mspecies Nsamples M frequency_offset x_positions gammabar

%%%% General Parameters %%%%
%Nsamples = 128; % # points
if(ischar(N))
    Nsamples = str2double(N);
else
    Nsamples = N;
end

%%%% General Parameters %%%%
%Nsamples = 128; % # points
N = Nsamples;
Mspecies = 2;  % # species
frequency_offset = [0 -3.4]; % ppm, for frequency deviation 

fprintf(sprintf('Starting Species Separation for %d Species and %d Points\n', Mspecies, Nsamples));

%%%% Gradient Parameters %%%%
b = 0.3*1e-3; % s
time_step = 2*b/Nsamples; % s
G = 0.5*1e-3; % T/m
field_strength = 1.5; % T
gammabar = 42.577*1e6; % hz/T %CAMBIAR!!!!!
dkx = gammabar*G*time_step; % 1/m %REVISAR

frequency_offset = frequency_offset_ppm*gammabar*1e-6*field_strength; % Hz

% Automatic calculation for two-point Dixon
% Considers first two species
Delta_f = abs(frequency_offset(1) - frequency_offset(2));
TE = zeros(1,3);
TE(1) = 1000*pi/(2*pi*Delta_f); %ms
TE(2) = 2*TE(1); %ms
TE(3) = 2.5*TE(2); %ms

%TE = [2.3 4.6 5]; % Two - Point Dixon
%TE = [1 2 4]; % ms % 
TE = [0.9 1.9 2.9]; % ms % As described in Reeder 2004, pg. 38
%TE = [3.96 4.76 5.55]; % ms % Aligned combination
%TE = [3.37 4.17 4.96]; % ms % Quadrature Combination

TE = TE.*1e-3; % s
TE

InitialDelay = TE - 2*b;

% 1. Traj. Type (only == 1)
% 2. Time Step [time_step] (s)
% 3. Initial Delay (s)
% 4. Negative Lobe duration (s) [b]
% 5. Gradient Amplitude (T/m)
% 6. Number of samples
gradient_parameters = [ 1 time_step InitialDelay(1) b G Nsamples;
                        1 time_step InitialDelay(2) b G Nsamples;
                        1 time_step InitialDelay(3) b G Nsamples];           


%%%% Space Parameters
FOV = 1/dkx; % m
dx = FOV/Nsamples; % m
x_positions = ((-FOV/2):dx:(FOV/2 - dx)); %m

[m, field_inhomogeneity] = get_object(x_positions,Mspecies); % Supports Mspecies = 2 

% Number of acquisitions
Nacq = size(gradient_parameters,1);
kx = zeros(Nacq, Nsamples); % Hz/m
kf = zeros(Nacq, Nsamples); 
M = zeros(size(kx));       
                    
m_fft_estimate = zeros(Nacq, Nsamples);      

kSpaceColors = rand(Nacq, 3); % 3 is for RGB

% Acquire and save k-space
for ii = 1:Nacq
    [kx_aux, kf_aux] = get_trajectory(gradient_parameters(ii,:));        
    % kf_aux = mean(kf_aux)*ones(size(kf_aux));
    Maux = get_acquisition(m, x_positions, frequency_offset, field_inhomogeneity, kx_aux, kf_aux);
    M(ii,:) = Maux;
    kx(ii,:) = kx_aux;
    kf(ii,:) = kf_aux;
    m_fft_estimate(ii, :) = fftshift(ifft(fftshift(M(ii,:))));
    %m_fft_estimate(ii, :) = fftshift(ifft(M(ii,:))); % Probably not
    %m_fft_estimate(ii, :) = ifft(M(ii,:)); % DISCARDED
    figure(1)
    hold on
        plot(kx_aux, kf_aux, '.-', 'Color', kSpaceColors(ii,:))
    hold off
    set(gcf, 'Color', [1 1 1]);
end

clear kx_aux kf_aux Maux

figure(1)
legend(num2str((1:Nacq)'));
title('k-space trajectories');
ylabel('k_f');
xlabel('k_x');
drawnow

figure(2)
clf;
subplot(3,1,1);
title('Phantom and Initial estimate of object');
plot(x_positions,real(m.'),'b-',...
    x_positions,real(m_fft_estimate(1, :)),'ro',...
    x_positions,real(m_fft_estimate(2, :)),'k+');
legend('re(m)', 're(m)', 're(m1)', 're(m2)');
subplot(3,1,2);
plot(x_positions,imag(m.'),'b-',...
    x_positions,imag(m_fft_estimate(1, :)),'ro',...
    x_positions,imag(m_fft_estimate(2, :)),'k+');
legend('im(m)', 'im(m)', 'im(m1)', 'im(m2)');
subplot(3,1,3);
plot(x_positions,abs(m.'),'b-',...
    x_positions,abs(m_fft_estimate(1, :)),'ro',...
    x_positions,abs(m_fft_estimate(2, :)),'k+');
legend('|m1|','|m2|', '|m1est|', '|m2est|');
set(gcf, 'Color', [1 1 1]);
xlabel('[cm]');
drawnow

%pause

%%%%%%%%%%%%%%%% 2 Point Dixon %%%%%%%%%%%%%%%
if(TE(2)/TE(1)==2)
    disp('Dixon separation detected');
    m_dixon = [(m_fft_estimate(1,:) + m_fft_estimate(2,:)); 
               (-m_fft_estimate(1,:) + m_fft_estimate(2,:));].*0.5;
    m_dixon = ExtTwoPointDixon(m_fft_estimate);       
    %figure(8)
    %clf
    %PlotObject(m, m_dixon, field_inhomogeneity, field_inhomogeneity, 'Two-Point Dixon' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% IDEAL %%%%%%%%%%%%%%%%%%%
TEp = 0.5.*kf(:,length(kf)/2) + 0.5.*kf(:,length(kf)/2 + 1);
TE = TEp';
%TEp = kf(:,length(kf)/2)./gammabar; 

[m_ideal, phi_ideal] = IDEALseparation(m_fft_estimate, TE', frequency_offset, field_inhomogeneity);
%figure(6)
%clf
%PlotObject(m, m_ideal, field_inhomogeneity, phi_ideal, 'IDEAL with known Field Map');
%figure(5)
%pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% We need to estimate real + imag components
% for each species (2*M*N), plus the field
% inhomogenity (N).
% The order is [spe1_real(1xN) spe1_imag(1xN) ... speM_real(1xN) speM_imag(1xN) phi_x(1xN)]
optimization_vector = zeros(1,(2*Mspecies + 1)*Nsamples);
initial_optimization_vector = zeros(1,(2*Mspecies + 1)*Nsamples);

% Fill positions with real and imag values of initial estimate
%initial_optimization_vector(1:(2*Mspecies*Nsamples)) = ...
%                repmat([real(m_fft_estimate(1,:)) imag(m_fft_estimate(1,:))],1,Mspecies);
% The field inhomogeneity initial estimate is the one from IDEAL
initial_optimization_vector((2*Mspecies*Nsamples + 1):end) = phi_ideal; 

for ii = 1:Mspecies            
    initial_optimization_vector((2*(ii-1)*Nsamples + 1):(2*ii*Nsamples)) = ...
                    [real(m_ideal(ii,:)) imag(m_ideal(ii,:))];
end

% for ii = 1:Mspecies            
%     initial_optimization_vector((2*(ii-1)*Nsamples + 1):(2*ii*Nsamples)) = ...
%                     [real(m(ii,:)) imag(m(ii,:))];
% end
%initial_optimization_vector((2*Mspecies*Nsamples + 1):end) = field_inhomogeneity;
%initial_optimization_vector((2*Mspecies*Nsamples + 1):end) = zeros(1,Nsamples);

% TypicalX = zeros(size(optimization_vector));
% for ii = 1:Mspecies
%    TypicalX(((ii - 1)*2*Nsamples + 1):(ii*2*Nsamples)) =  [real(m(ii,:)) imag(m(ii,:))];
% end
% 
% if(length(TypicalX)~=((2*Mspecies + 1)*Nsamples))
%     warning('Initial optimization vector seems to be incomplete');
% end    

ubound = repmat(max(abs(m_fft_estimate(1,:))), 1, (2*Mspecies + 1)*Nsamples).*1.5;
max_field_inhomogeneity = max(abs(frequency_offset))/3;
ubound((2*Mspecies*Nsamples + 1):(2*Mspecies*Nsamples + Nsamples)) = repmat(max_field_inhomogeneity,1,Nsamples);
lbound = -ubound;

Niter = 100;
GradOffiter = 0; %!!
GradOniter = 10;

clear m_ideal phi_ideal

close all;

for ii = 1:Niter    
%     % First 'medium-scale: SQP, Quasi-Newton, line-search'
%     options = optimset('TypicalX',TypicalX,...
%         'TolFun',1e-5,'MaxIter',GradOffiter,'TolX',1e-5,'MaxFunEvals',5000,...
%         'GradObj','off',...
%         'Display','iter',...
%         'Algorithm', 'interior-point');
%     disp('Optimization Gradients OFF');
%     tic;
%     [optimization_vector,fval,exiflag,output] = fmincon(@(x)GradientFunction2(x, kx, kf, M, frequency_offset, x_positions),...
%         initial_optimization_vector,[],[],[],[],lbound,ubound,[],options);
%     toc; output
%     
%     initial_optimization_vector = optimization_vector; % Ready for next optimization
%   
%     m_opt = zeros(Mspecies, N);
% 
%     for jj = 1:Mspecies
%         initial_index = ((jj - 1)*2*N);
%         m_opt(jj,:) = optimization_vector((initial_index + 1):(initial_index + N)) + ...
%                 1i*optimization_vector((initial_index + N + 1):(initial_index + 2*N));
%     end
%     phi_x_opt = optimization_vector(2*Mspecies*N+1:end).';
%     
     ShouldIPlotGradients = 0;    
%     
%     % Plotting
%     figure(5)
%     PlotObject(m,m_opt,field_inhomogeneity,phi_x_opt.');        
%     if ShouldIPlotGradients==1,
%         PlotGradientFunction(x_positions, optimization_vector);
%     end            
%     pause(1);
    
    % Second: 'large-scale: trust-region reflective Newton'
    options = optimset('TolFun',1e-5,'MaxIter',GradOniter,'TolX',1e-5,'MaxFunEvals',5000,...
        'GradObj','on',...
        'Display','iter',...
        'Algorithm', 'interior-point');
    disp('Optimization Gradients ON');
    tic;
    [optimization_vector,fval,exiflag,output] = fmincon(@(x)GradientFunction2(x, kx, kf, M, frequency_offset, x_positions),...
        initial_optimization_vector,[],[],[],[],lbound,ubound,[],options);
    toc; output
    
    if ShouldIPlotGradients==1,
        PlotGradientFunction(x_positions, optimization_vector)
    end
    
    initial_optimization_vector = optimization_vector; % Ready for next optimization
  
    m_opt = zeros(Mspecies, N);

    for jj = 1:Mspecies
        initial_index = ((jj - 1)*2*N);
        m_opt(jj,:) = optimization_vector((initial_index + 1):(initial_index + N)) + ...
                1i*optimization_vector((initial_index + N + 1):(initial_index + 2*N));
    end
    phi_x_opt = optimization_vector(2*Mspecies*N+1:end).';

    % Plotting
    %figure(5)
    %PlotObject(m,m_opt,field_inhomogeneity,phi_x_opt.');    
    if ShouldIPlotGradients==1,
        PlotGradientFunction(x_positions, optimization_vector);
    end            
    pause(1);
    %save('CurrentOptimizationState','m','m_opt','field_inhomogeneity','phi_x_opt');
end;

disp('Max. number of iterations reached. Exiting');