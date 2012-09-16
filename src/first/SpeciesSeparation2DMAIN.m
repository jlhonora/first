% % % Chemical Species Separation using k-space formulation % % %
%
% JLH - jlhonora@ing.puc.cl
% Biomedical Imaging Center
% Pontificia Universidad Catolica - Chile
% Last Updated: Jun. 2011
% 2D, R2 capability. Multipeak version is available but not documented.
%
% DISCLAIMER: The following code is provided as is and cannot be used 
% for patient diagnosis. No support is guaranteed, although you can contact 
% the authors at the given e-mail adresses or at the Biomedical Imaging 
% Center:
%     Centro de Imagenes Biomedicas
%     Pontificia Universidad Catolica de Chile
%     Vicuna Mackenna 4860, Macul.
%     Santiago, Chile.
%     Phone: +56 2 3548468.
%
% % % Description
%     The following script is used as the main function for separating
%     chemical species with the FIRST algorithm. Implementations for
%     two-point Dixon and IDEAL are also included. A routine for generating
%     a R2 decay map is computed. This script makes extensive use of
%     the Parallel Computing, Optimization and Image Processing toolboxes
%     provided with MATLAB.
%     
% % % FIRST Algorithm
%     FIRST uses a complete MR signal model considering different chemical
%     species, with addition of R2* correction and multipeak modelling of 
%     fat (not in this script). The estimation of the unknown variables is
%     done by fitting the raw data of the acquisitions to the signal model 
%     using an interior point algorithm available in the fmincon function. 
%     The derivates of the objective function are calculated internally by
%     MATLAB or provided with an explicit calculation. We will refer to the
%     first case as GradOff and to the second case as GradOn.
%         Two modes of solving the minimization problem can be used. The
%     first one is separating the problem in two steps. The first step
%     consists in estimating the species and field map. For the second
%     step, the field map is considered as correct and a single R2* value
%     is considered for all species in a single pixel, so the estimation of
%     R2* and the correction of the species intensity by this decay is done
%     simultaneously. The second way of solving the problem is to estimate
%     species, field inhomogeneity and R2* simultaneously.
%         Because it fits the model in the k-space domain, the acquisition
%     trajectories and time map are needed for solving the problem. This
%     script works for cartesian acquisitions with no speed improvements,
%     so the time map needs to be the same for every line in the phase
%     encoding direction. Nevertheless, it can use different trajectories
%     for every acquisition. This script assumes a horizontal readout. For
%     vertical readout a transposition of the k-space data is recommended.
%     Please note that transposition for matrixes with complex values is
%     done by Atransposed = A.'; <- The point is strictly necessary.
%         The script takes about 5 hours to solve a data set of 4 images of
%     size 256 x 256 pixels, so be patient!
%
% % % List of third-party functions (not provided by MATLAB)
%
% - Get2DAcquisition 
%     Returns a structure containing the dicom info, acquisition images,
%     k-space trajectories, raw data, TE and all other data needed
%     throughout the script.
% - Plot2DObject
%     Useful function for displaying a stack of images, with every stack
%     thought as a different chemical species.
% - Plot1DObject
%     Same as Plot2DObject but for 1D objects. Useful for printing only one
%     line of the object.
% - PlotGradientFunction
%     Plots the calculated gradients of the objective function. Its use is
%     not recommended.
% - GetR2fit
%     Calculates the R2 signal decay of a stack of images.
% - ExtTwoPointDixon
%     Two-point Dixon separation with field inhomogeneity calculation.
% - IDEALT2sep
%     IDEAL algorithm with T2* estimation capabilities. The original
%     version of IDEAL can be calculated with this script.
% - GradientFunction5
%     Auxiliar function for fmincon to estimate species and field
%     inhomogeneity.
% - GradientFunctionR2
%     Auxiliar function for fmincon to estimate species and R2 decay.
%     The v2 version considers a variable time map for the R2* decay and it
%     is slightly faster.
% - GradientFunction7
%     Auxiliar function for fmincon to estimate species, field inhomogeneity
%     and R2 decay. The v2 version considers a variable time map for the R2* 
%     decay and it is slightly faster.
%
%     Further considerations
%     The main script is highly documented, but we suggest to browse
%     through the third party functions written for specific functioning
%     and customization.


disp('Chemical Species Separation using k-space formulation');
disp('2D Version - JLH, Jun 2011');

path(path, '../general');
path(path, '../ideal');
path(path, '../data');

global kx ky kf Mspecies Nx Ny M frequency_offset x_positions y_positions gammabar n_doneiters

%%%% General Parameters %%%%
% These parameters modify the behaviour of the script.

% Number of iterations for the whole process.
Niter = 1;

% Number of iterations for the first step without explicitly calculating 
% the gradients of the objective function.
GradOffiter = 0; 
% Number of iterations for the first step explicitly calculating 
% the gradients of the objective function.
GradOniter = 280;
% Number of iterations for the second step without explicitly calculating 
% the gradients of the objective function.
GradOffiterR2 = 0;
% Number of iterations for the second step explicitly calculating 
% the gradients of the objective function.
GradOniterR2 = 0;
% Number of iterations for jointly estimating all variables.
GradOniterJ = 780;

% Range of lines for which the algorithm will be applied. For example, 
% if range = [1 256], the algorithm will be applied to the 
% lines from 1 to 256
%%%%
range = [1 128];
row_range = [1:128];
col_range = [1:128];
%%%%

% Parallelization of the solution. Refer to the documentation for
% customization. 
UseParallel = 0; Npools = 4;

% If set to 1, the output from the iterations is printed, along with
% warnings and others. If set to 0, warnings are deactivated.
HighConsoleOutput = 1;
% If parallel computation is used, HighConsoleOutput is set to false.
HighConsoleOutput = logical(HighConsoleOutput*(~UseParallel));

% Plots images and 1D objects for each iteration
WithGraphics = 1;

% If set to 1, pauses before entering the main optimization loop.
PauseBeforeOptim = 0;

% For plotting gradients of the objective function.
ShouldIPlotGradients = false; 

% Windowing of the raw data
ApplyHammingWindow = 0; 

% Multiply each line of k-space in order to warp the relevance of each
% k-space point. Set to false - use is not recommmended.
ApplyCoefficients = false;

% If set to 1, computes IDEAL with T2* estimation.
T2IDEAL = 1;

% Multipeak modelling of fat
MultiPeak = 0;

% Generate field map from acquisitions. If AcqForB0 = [3 2] considers
% images 3 and 2 for the calculation.
CalculateB0fromAcq = true; AcqForB0 = [3 2];

% Get the acquisition data. If set to 0 uses the acquisition struct
% available in workspace (if any).
IncludeSimulation = 1; %%% !!

% Adds AWGN to the input data;
WithNoise = 0;
NoiseSNR = 10;

% If no high console output, warnings are turned off.
warning on;
if(~HighConsoleOutput)
    warning off;
end

% Number of species in the acquisition.
Mspecies = 2;  % # species
% Resonating frequencies for each species are needed (in ppm). Fat is
% consider to have a resonating frequency of -3.5ppm from water (0 ppm).
frequency_offset_ppm = [0 -3.5];% -5.1]; % ppm, for frequency deviation 
% Field strength of the main coil.
field_strength = 1.5; % T
% Gamma-bar MR parameter
gammabar = 42.577*1e6; % hz/T 
% Frequency offset in Hertz
frequency_offset = frequency_offset_ppm*gammabar*1e-6*field_strength; % Hz

% Requested TEs in miliseconds. Note that the nominal TE for each 
% acquisition is different from the exact TE sampling value. 
%TE = [2.3 4.6 5]; % Two - Point Dixon
%TE = [1 2 4]; % ms % 
%TE = [0.9 1.9 2.9];% 4.5]; % ms % As described in Reeder 2004, pg. 38
%TE = [3.96 4.76 5.55]; % ms % Aligned combination
%TE = [3.37 4.17 4.96]; % ms % Quadrature Combination
%TE = [4.6 6.9 4.8]; % ms  % 11-01-17 % YY-MM-DD
%TE = [4.9 5.2 5.5 6.2]; %[6.2 4.6 5.2]; % ms  % 11-03-01
TE = [4.6 5.2 5.5 5.9 6.3 7.0]; % 4.75 4.9 5.05 5.35 5.5]; % ms  % 11-03-10
% TE = [4.4 4.6 5.4 6.5 8]; % 4.8 5.2 5.4 5.6 5.9 6.2 6.9 7.4 8 15]; % ms % 11-04-04
%TE = [4.6 4.8 6.2 7.5]; % 4.6 4.8 5 5.4 5.8 6.2 6.7 6.9 7.5 8 11 14]; % ms  % 11-04-11
TE = TE.*1e-3; % s

% If set to true, asummes the first two echo times as in-phase and
% out-of-phase images. 
IsDixon = false;
% Number of acquisitions.
Nacq = length(TE);

if(IncludeSimulation) % If IncludeSimulation == 0, then the input data is 
    % taken from the workspace. If not, the data is obtained from dicom and
    % raw data files.
    %
    % The following variables are used throughout the whole script: 
    %   AcqStructs : Explained above.
    %   M : A stack of complex matrixes where each layer is the fourier
    %       space for a certain acquisition.
    %   Mkxy : A stack where each layer is the corresponding layer of M but
    %          after a Fourier transform of the columns.
    %   m_fft_estimate : Stack of the 2D FFT of each layer of M. 
    %   kx : A matrix of Nacquisitions lines and Nx columns containing the
    %        kx trajectory of each acquisition (asuming cartesian acquisitions).
    %        For 2D optimization a stack of gridded trajectories will be
    %        needed.
    %   ky : Similar to kx.
    %   kf : Time map of the acquisition. Similar to kx. 
    %   x_positions : Vector of Ncols samples that range from -FOVx/2 to
    %                 FOVx/2
    %   y_positions : Similar to x_positions
    AcqStructs = [];
    M = [];
    Mkxy = [];
    m_fft_estimate = [];
    field_inhomogeneity = [];
    kx = [];
    ky = [];
    kf = [];
    x_positions = [];
    y_positions = [];
    
    % The acquisition directory contains the dicom images and the raw data
    % for each acquisitions. Along with the acquisition files, an "info" 
    % plain text file (with no extension) is needed with the following
    % format:
    %
    % TE MAG RE IM RAW
    % 4.60 IM_0249 IM_0250 IM_0251 raw_000
    % 4.80 IM_0255 IM_0256 IM_0257 raw_001
    % 5.00 IM_0261 IM_0262 IM_0263 raw_002
    %
    % In this example, the first line indicates an acquisition of nominal
    % TE of 4.6 ms, which magnitude, real and imaginary image files are
    % IM_0249, 250 and 251, respectively. If real and imaginary files are
    % not available, any string will do. Finally, raw_000 is the raw data
    % for this acquisition. 
    %
    %AcqDirec = '../../data/set1'; disp('InVivo Brain');
    AcqDirec = '../../data/set2'; disp('Bottle Phantom');
    %AcqDirec = 'D:\Dropbox\MR\Pruebas Resonador\InVivo 11-04-19-1'; disp('Thigh Acquisition');
    %AcqDirec = 'D:\Dropbox\MR\Pruebas Resonador\InVivo 11-04-19-2'; disp('Brain Acquisition');
    RequestedTEs = TE;
    [AcqStructs] = Get2DAcquisition(AcqDirec, RequestedTEs, 'NoDicom'); %'DicomOnly', 'NoDicom'
    Nrows = size(AcqStructs(1).Real,1); Ny = Nrows;
    Ncols = size(AcqStructs(1).Real,2); Nx = Ncols;
    M = zeros(Nrows, Ncols, Nacq);
    Mkxy = zeros(Nrows, Ncols, Nacq);
    m_fft_estimate = zeros(Nrows, Ncols, Nacq);
    kx = zeros(Nacq, Ncols);
    ky = zeros(Nacq, Nrows);
    kf = zeros(Nacq, Ncols);
    
    if(ApplyHammingWindow)
        HammingWindow = repmat(hamming(Nx).'.^0.3, Ny,1);
        HammingWindow = HammingWindow.*HammingWindow'; 
    end
    
    for ii = 1:Nacq
        M(:,:,ii) = AcqStructs(ii).KSpace;
            % Windowing of the k-space data
            if(ApplyHammingWindow)
                M(:,:,ii) = M(:,:,ii).*HammingWindow;
            end
        % Add AWGN to the input data
        if(WithNoise)
            for jj = 1:Nrows
                M(jj,:,ii) = awgn(M(jj,:,ii), NoiseSNR);
            end
        end
        % Fourier transform for columns only    
        Mkxy(:,:,ii) = fftshift(ifft(fftshift(M(:,:,ii))));
            
        % Fill k-space trajectories. I suggest changing for a stack of
        % matrixes to allow 2D optimization
        kx(ii,:) = AcqStructs(ii).KTraj.kx;
        ky(ii,:) = AcqStructs(ii).KTraj.ky;
        kf(ii,:) = AcqStructs(ii).KTraj.kf;
        % The echo time is the middle of kf
        TE(ii) = AcqStructs(ii).KTraj.kf(Nx/2+1);
        % Fourier transform of the raw data
        m_fft_estimate(:,:,ii) = fftshift(ifft2(fftshift(M(:,:,ii))));
    end
    x_positions = AcqStructs(1).Pos.x_pos;
    y_positions = AcqStructs(1).Pos.y_pos;
    % Initial estimate of field inhomogeneity (arbitrary)
    field_inhomogeneity = repmat(max(abs(frequency_offset))/3, Nrows, Ncols);
end

% format long
% TE
% format short

clear kx_aux ky_aux kf_aux Maux HammingWindow

% Calculate the field map estimate from the acquisitions.
if(CalculateB0fromAcq)
    idx1 = AcqForB0(1);
    Im1 = AcqStructs(idx1).Real + 1i*AcqStructs(idx1).Imag; 
    TE1 = AcqStructs(idx1).TE;
    idx2 = AcqForB0(2);
    Im2 = AcqStructs(idx2).Real + 1i*AcqStructs(idx2).Imag; 
    TE2 = AcqStructs(idx2).TE;    
    field_inhomogeneity = angle(Im1.*conj(Im2))./(2*pi*(TE1 - TE2)*1e-3); % Hz
end

% Mask the estimation and assume field_map = 0 if |field_map| > 110
fmask = AcqStructs(1).Mag>5e-2.*mean(mean(AcqStructs(1).Mag));
fmask = bwareaopen(fmask,10,4);
fmask = imfill(fmask,'holes');
field_inhomogeneity(~fmask) = 0;
field_inhomogeneity(field_inhomogeneity>110) = 0;
field_inhomogeneity(field_inhomogeneity<-260) = 0;
% Remove phase wrap
field_inhomogeneity(field_inhomogeneity<-110) = field_inhomogeneity(field_inhomogeneity<-110) - frequency_offset(2);
% Filter
field_inhomogeneity = imfilter(field_inhomogeneity, fspecial('gaussian', [3 3], 2));

if(WithGraphics)
    Plot2DObject(m_fft_estimate, field_inhomogeneity, 'i2DFT-estim', 'Acquisition');
    m_fft2 = zeros(size(m_fft_estimate));
    for ii = 1:Nacq
        m_fft2(:,:,ii) = AcqStructs(ii).Real + 1i*AcqStructs(ii).Imag;
    end
    Plot2DObject(m_fft2, field_inhomogeneity, 'i2DFT', 'Acquisition');
end

%%%%%%%%%%%%% R2* Fitting from Acq %%%%%%%%%%%
try
    R2_fitted; % check if there's one in workspace
catch
    fprintf('Calculating R2* from exponential fit\n');
    tic
    R2_fitted = GetR2fit(AcqStructs, TE);
    R2_fitted(R2_fitted<0) = 0;
    R2_fitted = R2_fitted.*(AcqStructs(1).Mag>(0.01*(mean(mean(AcqStructs(1).Mag)))));
    toc
end

%%%%%%%%%%%%%%%% 2 Point Dixon %%%%%%%%%%%%%%%
if(TE(2)/TE(1)==2 || IsDixon)
    disp('Dixon separation detected');
    [m_dixon, phi_dixon] = ExtTwoPointDixon(m_fft_estimate);   
    if(WithGraphics)
        Plot2DObject(m_dixon, phi_dixon, 'Two-Point Dixon' );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% IDEAL %%%%%%%%%%%%%%%%%%%
% TEp = 0.5.*kf(:,length(kf)/2) + 0.5.*kf(:,length(kf)/2 + 1);
% TE = TEp';

%[m_ideal, phi_ideal] = IDEALsep(m_fft_estimate, TE, frequency_offset);
%[m_ideals, phi_ideals, t2_map] = IDEALT2sep(m_fft_estimate, TE, frequency_offset, phi_ideal);

%clear m_ideal;

% If IDEAL already available in workspace then it is not calculated again.
% Due to low rank matrix inversion, IDEAL separation throws a lot of
% warnings, which are disabled.
warning off
try
    m_ideals = m_ideal; phi_ideals = phi_ideal; t2_map = zeros(size(phi_ideal));
catch
    tic
    if(T2IDEAL)    
        [m_ideal, phi_ideal, R2_map] = IDEALT2sep(m_fft_estimate, TE, frequency_offset, field_inhomogeneity, 'SmoothFieldMap', R2_fitted);
    else
        [m_ideal, phi_ideal] = IDEALsep(m_fft_estimate, TE, frequency_offset, field_inhomogeneity, 'SmoothFieldMap');%, field_inhomogeneity);
        R2_map = zeros(size(phi_ideal));
    end
    % Total processing time of IDEAL
    t_ideal = toc;
end

% If console output is needed, the warning is turned on again.
if(HighConsoleOutput)
    warning on;
end

% Maximum field inhomogeneity in terms of the resonance frequency of the
% species.
max_field_inhomogeneity = max(max(abs(frequency_offset)))/1.8;
% Field values saturation
phi_ideal((phi_ideal)>max_field_inhomogeneity) = max_field_inhomogeneity;
phi_ideal((phi_ideal)<-max_field_inhomogeneity) = -max_field_inhomogeneity;

m_ideal(abs(m_ideal)>10) = 0; %%% !!!!!! %%%

if(WithGraphics)
    if(T2IDEAL)
        Plot2DObject(m_ideal, phi_ideal, 'T2-IDEAL with unknown Field Map');
        figure(gcf + 1)
        set(gcf, 'Color', [1 1 1]);
        imshow(R2_map,[0 300]);
        colorbar;
        title('T2-IDEAL - R2* Map');
    else
        Plot2DObject(m_ideal, phi_ideal, 'IDEAL with unknown Field Map');
    end
end

%m_ideal = m_ideals;
%phi_ideal = phi_ideals;
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Main Optimization Loop %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsamples = Nx;
N = Nsamples;
% Upper and lower bounds for species and field inhomogeneity.
ubound = repmat(max(max(abs(m_fft_estimate(:,:,1)))), 1, 2*(Mspecies + 1)*Nsamples).*1.5;
ubound(2*Mspecies*Nsamples + (1:Nsamples)) = repmat(max_field_inhomogeneity*0.9,1,Nsamples);
lbound = -ubound;

% Upper and lower bounds for R2*
max_R2 = 350;
min_R2 = 0;
ubound(2*(Mspecies)*N + N + (1:N)) = repmat(max_R2,1,N);
lbound(2*(Mspecies)*N + N + (1:N)) = repmat(min_R2,1,N);

if(PauseBeforeOptim)
    %return
    if(WithGraphics)
        figure(5)
    end
    pause
end

% Stack of complex matrixes for each species
m_total = zeros(size(m_ideal));
% Separate stacks for water and fat are used, needed for the Parallel
% Computing Toolbox to work
m_total_water = zeros(size(phi_ideal));
m_total_fat = zeros(size(phi_ideal));
% Matrixes for field and R2
phi_total = zeros(size(field_inhomogeneity));
R2_total = zeros(size(field_inhomogeneity));
% Stores final objective value for each column
vobj = zeros(1,Ny);

%m = zeros(Ny+1,Nx,Mspecies);
%m(2:end,:,:) = m_ideal;
%phi_aux = zeros(Ny+1,Nx);
%field_inhomogeneity = phi_ideal;

% Not used - See the ApplyCoeficients variable.
coefs = repmat((kx(1,:)).^2./(kx(end)/1.8) + 1, Nx, 1);
if(~ApplyCoefficients)
    coefs = ones(Nx, Ny);    
end

disp('Entering main 1D optimization loop');

% TypicalX for Object. Used by fmincon to get an idea of the magnitude of
% the variables.
TypicalX = zeros(1, 2*(Mspecies + 1)*Nsamples);
for ii = 1:Mspecies            
    TypicalX((2*(ii-1)*Nsamples + 1):(2*ii*Nsamples)) = ...
                    [real(m_ideal(Ny/2,:,ii)) imag(m_ideal(Ny/2,:,ii))];
end
% TypicalX for Field
TypicalX(2*Mspecies*Nsamples + (1:Nsamples)) = field_inhomogeneity(Ny/2,:);
% TypicalX for R2
R2fat = 120; % 
TypicalX(2*(Mspecies)*Nsamples + Nsamples + (1:Nsamples)) = repmat(R2fat,1,Nx);

if(UseParallel)
    if(matlabpool('size')>0)
        matlabpool close;
    end
    matlabpool('open', Npools);
    n_doneiters = 0;
    fprintf('Using %d parallel workers\n',matlabpool('size'));
else
end
tic

%for yy = range
for yy = range(1):range(2) % Ranges of the processed lines.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    Myy = zeros(Nacq, Nsamples);
    % Select the corresponding line of k-space
    for ii = 1:Nacq
        Myy(ii,:) = Mkxy(yy,:,ii);
    end
    
    % We need to estimate real + imag components
    % for each species (2*M*N), plus the field
    % inhomogenity (N).
    % The order is:
    % [spe1_real(1xN) spe1_imag(1xN) ... speM_real(1xN) speM_imag(1xN) phi_x(1xN) R2_x(1xN)]
    optimization_vector = zeros(1,2*(Mspecies + 1)*Nsamples);
    initial_optimization_vector = zeros(1,2*(Mspecies + 1)*Nsamples);
    R2_opt = zeros(1,Nsamples);
    phi_opt = zeros(1,Nsamples);
    m_opt = zeros(Mspecies,Nsamples);
    % Fill positions with real and imag values of initial estimate
    %initial_optimization_vector(1:(2*Mspecies*Nsamples)) = ...
    %               repmat([real(m_fft_estimate(yy,:,1)) imag(m_fft_estimate(yy,:,1))],1,Mspecies);            
    % The field inhomogeneity initial estimate:
    %initial_optimization_vector((2*Mspecies*Nsamples)+ (1:N)) = field_inhomogeneity(yy,:); % zeros(1,Nsamples); % phi_ideal(yy,:); % field_inhomogeneity(yy,:);
    %initial_optimization_vector((2*Mspecies*Nsamples + 1):end) = initial_optimization_vector((2*Mspecies*Nsamples + 1):end) ...
    %                + (rand(1,Nsamples) - 0.5)*10;
    %initial_optimization_vector((2*Mspecies*Nsamples) + (1:N)) = phi_ideal(yy,:);
    for ii = 1:Mspecies            
        initial_optimization_vector((2*(ii-1)*Nsamples + 1):(2*ii*Nsamples)) = ...
                        [real(m_ideal(yy,:,ii)) imag(m_ideal(yy,:,ii))];
    end
%     for ii = 1:Mspecies            
%         initial_optimization_vector((2*(ii-1)*Nsamples + 1):(2*ii*Nsamples)) = ...
%                         [real(m_total(yy,:,ii)) imag(m_total(yy,:,ii))];
%     end

%    initial_optimization_vector((2*Mspecies*Nsamples) + (1:N)) = field_inhomogeneity(yy,:);
    %phi_total(yy,:); % zeros(1,Nsamples); % phi_ideal(yy,:); % field_inhomogeneity(yy,:);

    % R2
    %initial_optimization_vector(2*(Mspecies)*Nx + Nx + (1:Nx)) = repmat(R2fat,1,Nx);
    initial_optimization_vector(2*(Mspecies)*Nx + Nx + (1:Nx)) = R2_fitted(yy,:);
    
    m_init = zeros(Mspecies, N);
    for jj = 1:Mspecies
        initial_index = ((jj - 1)*2*N);
        m_init(jj,:) = initial_optimization_vector((initial_index + 1):(initial_index + N)) + ...
                1i*initial_optimization_vector((initial_index + N + 1):(initial_index + 2*N));
    end
    phi_init = initial_optimization_vector(2*Mspecies*N+(1:N)).';

    R2_init = initial_optimization_vector(2*(Mspecies)*Nx + Nx + (1:Nx));
    
    for ii = 1:Niter        
        if(GradOffiter>0) % Calculate species and field without calculating gradients       
            options = optimset('TolFun',1e-6,'MaxIter',GradOffiter,'TolX',1e-8,'MaxFunEvals',30e4,...
                'GradObj','off',...
                'TypicalX', TypicalX,...
                ...'Display','iter',...
                'Hessian','lbfgs',...
                'Algorithm', 'interior-point');
            if(HighConsoleOutput)
                options.Display = 'iter'; % Prints each iteration
            end
            disp('Optimization Gradients OFF - no R2');
            %tic;
            [optimization_vector,fval,exiflag,output] = fmincon(@(x)GradientFunction5(x, kx, kf, Myy, frequency_offset, x_positions),...
                initial_optimization_vector,[],[],[],[],lbound,ubound,[],options);
            %toc;            
            if(HighConsoleOutput)
                output % Prints the optimization output
            end            
            initial_optimization_vector = optimization_vector; % Ready for next optimization
            m_opt = zeros(Mspecies, N);
            for jj = 1:Mspecies
                initial_index = ((jj - 1)*2*N);
                m_opt(jj,:) = optimization_vector((initial_index + 1):(initial_index + N)) + ...
                        1i*optimization_vector((initial_index + N + 1):(initial_index + 2*N));
            end
            phi_x_opt = optimization_vector(2*Mspecies*N+(1:N)).';
            %R2_opt = optimization_vector(2*Mspecies*Nsamples + N + (1:N));         
            if(WithGraphics)
                figure(5)
                %Plot1DObject(m(yy,:,:),m_opt,field_inhomogeneity(yy,:),phi_x_opt.');
                Plot1DObject(m_init,m_opt,phi_init.',phi_x_opt.');
                %Plot1DObject(m_ideal(yy,:,:),m_opt,field_inhomogeneity(yy,:),phi_x_opt.');
                pause(1);
            end
            if(output.iterations < GradOffiter)
                disp('Warning: less gradient-off iterations than requested');
                %break;
            end        
        end
        
        if(GradOniter>0) % Calculate species and field calculating gradients
            % Second: 'large-scale: trust-region reflective Newton'
            options = optimset('TolFun',1e-6,'MaxIter',GradOniter,'TolX',1e-8,'MaxFunEvals',5000,...
                'GradObj','on',...
                'TypicalX', TypicalX,...
                ...'Display','iter',...
                'Hessian','lbfgs',...
                ...'SubproblemAlgorithm', 'cg',...
                ...'DerivativeCheck', 'on',...
                ...'FinDiffType', 'central',...
                ...'LargeScale', 'on',...
                ...'ScaleProblem', 'none',...
                'Algorithm', 'interior-point');
            if(HighConsoleOutput)
                options.Display = 'iter';
            end
            disp('Optimization Gradients ON - No R2');
            [optimization_vector,fval,exiflag,output] = fmincon(@(x)GradientFunction5(x, kx, kf, Myy, frequency_offset, x_positions),...
                    initial_optimization_vector,[],[],[],[],lbound,ubound,[],options);
            if(HighConsoleOutput)
                output
            end
            if(ShouldIPlotGradients)
                PlotGradientFunction(x_positions, optimization_vector)
            end
            initial_optimization_vector = optimization_vector; % Ready for next optimization
            m_opt = zeros(Mspecies, Nsamples);
            for jj = 1:Mspecies
                initial_index = ((jj - 1)*2*Nsamples);
                m_opt(jj,:) = optimization_vector((initial_index + 1):(initial_index + Nsamples)) + ...
                        1i*optimization_vector((initial_index + Nsamples + 1):(initial_index + 2*Nsamples));
            end
            phi_x_opt = optimization_vector(2*Mspecies*N+(1:N)).';
            %R2_opt = optimization_vector(2*(Mspecies)*N + N + (1:N));
            % Plotting
            if(WithGraphics)
                figure(5)
                %Plot1DObject(m(yy,:,:),m_opt,field_inhomogeneity(yy,:),phi_x_opt.');
                Plot1DObject(m_init,m_opt,field_inhomogeneity(yy,:),phi_x_opt.'); %phi_init.'
                pause(1);
            end
            if(ShouldIPlotGradients)
                PlotGradientFunction(x_positions, optimization_vector);
            end
            if(output.iterations < GradOniter)
                disp('Warning: less gradient-on iterations than requested');
                break;
            end
%             if(SaveIntermediate == 1 && ~UseParallel)
%                 save('LastOptimS2', 'm_total', 'phi_total', 'R2_total', 'vobj'); 
%             end
        end

        %%%%%%%%%% R2 Optim %%%%%%%%%
        if(GradOffiterR2>0) % Calculate species and R2* without calculating gradients
            % First 'medium-scale: SQP, Quasi-Newton, line-search'
            options = optimset('TolFun',1e-6,'MaxIter',GradOffiterR2,'TolX',1e-8,'MaxFunEvals',30e4,...
                'GradObj','off',...
                'TypicalX', TypicalX,...
                ...'Display','iter',...
                'Hessian','lbfgs',...
                'Algorithm', 'interior-point');
            if(HighConsoleOutput)
                options.Display = 'iter';
            end
            disp('Optimization Gradients OFF - with R2');
            %initial_optimization_vector(2*(Mspecies)*Nsamples + N + (1:N)) = R2_map(yy,:);
            [optimization_vector,fval,exiflag,output] = fmincon(@(x)GradientFunctionR2v2(x, kx, kf, Myy, frequency_offset, x_positions, phi_x_opt),...
                initial_optimization_vector,[],[],[],[],lbound,ubound,[],options);
            %toc;             
            if(HighConsoleOutput)
                output
            end            
            initial_optimization_vector = optimization_vector; % Ready for next optimization
            m_opt = zeros(Mspecies, N);
            for jj = 1:Mspecies
                initial_index = ((jj - 1)*2*N);
                m_opt(jj,:) = optimization_vector((initial_index + 1):(initial_index + N)) + ...
                        1i*optimization_vector((initial_index + N + 1):(initial_index + 2*N));
            end
            %phi_x_opt = optimization_vector(2*Mspecies*N+(1:N)).';
            R2_opt = optimization_vector(2*(Mspecies)*Nsamples + N + (1:N));
            % Plotting       
            if(WithGraphics)
                figure(5)
                %Plot1DObject(m(yy,:,:),m_opt,field_inhomogeneity(yy,:),phi_x_opt.');
                Plot1DObject(m_init,m_opt,phi_init.',phi_x_opt.');
                %Plot1DObject(m_ideal(yy,:,:),m_opt,field_inhomogeneity(yy,:),phi_x_opt.');
                pause(1);
            end
            if(output.iterations < GradOffiter)
                disp('Warning: less gradient-off iterations than requested');
                %break;
            end        
        end    
        
        %%%%%%%%%% R2 Optim  - GRAD ON %%%%%%%%%
        if(GradOniterR2>0) % Calculate species and R2* calculating gradients
            options = optimset('TolFun',1e-6,'MaxIter',GradOniterR2,'TolX',1e-12,'MaxFunEvals',30e4,...
                'GradObj','on',...
                'TypicalX', TypicalX,...
                ...'Display','iter',...
                ...'DerivativeCheck', 'on',...
                ...'LargeScale', 'on',...
                'Hessian','lbfgs',...
                'Algorithm', 'interior-point');
            if(HighConsoleOutput)
                options.Display = 'iter';
            end
            disp('Optimization Gradients ON - with R2');            
            if(GradOffiterR2==0)
                %initial_optimization_vector(2*(Mspecies)*N + N + (1:N)) = zeros(1,N);
                %optimization_vector(2*(Mspecies)*N + N + (1:N)) = zeros(1,N);
            end
            [optimization_vector,fval,exiflag,output] = fmincon(@(x)GradientFunctionR2v2(x, kx, kf, Myy, frequency_offset, x_positions, phi_x_opt),...
                initial_optimization_vector,[],[],[],[],lbound,ubound,[],options);        
            if(HighConsoleOutput)
                output
            end            
            initial_optimization_vector = optimization_vector; % Ready for next optimization
            m_opt = zeros(Mspecies, N);
            for jj = 1:Mspecies
                initial_index = ((jj - 1)*2*N);
                m_opt(jj,:) = optimization_vector((initial_index + 1):(initial_index + N)) + ...
                        1i*optimization_vector((initial_index + N + 1):(initial_index + 2*N));
            end
            %phi_x_opt = optimization_vector(2*Mspecies*N+(1:N)).';
            R2_opt = optimization_vector(2*(Mspecies)*Nsamples + N + (1:N));
            if(WithGraphics)
                figure(5)
                Plot1DObject(m_init,m_opt,phi_init.',phi_x_opt.');
                pause(1);
            end
            if(output.iterations < GradOniterR2)
                %disp('Warning: less gradient-on iterations than requested');
                %break;
            end
%             if(SaveIntermediate == 1 && ~UseParallel)
%                 save('LastOptimS4', 'm_total', 'phi_total', 'R2_total', 'vobj'); 
%             end
        end
        %%%%%%%%%% R2 and PHI Optim  - GRAD ON %%%%%%%%%
        if(GradOniterJ>0) % Calculate species, field and R2* calculating gradients
            options = optimset('TolFun',1e-6,'MaxIter',GradOniterJ,'TolX',1e-12,'MaxFunEvals',30e4,...
                'GradObj','on',...
                'TypicalX', TypicalX,...
                ...'Display','iter',...
                ...'DerivativeCheck', 'on',...
                ...'LargeScale', 'on',...
                'Hessian','lbfgs',...
                'Algorithm', 'interior-point');
            if(HighConsoleOutput)
                options.Display = 'iter';
            end
            disp('Optimization Gradients ON - with R2 and Phi');            
            if(GradOffiterR2==0)
                %initial_optimization_vector(2*(Mspecies)*N + N + (1:N)) = zeros(1,N);
                %optimization_vector(2*(Mspecies)*N + N + (1:N)) = zeros(1,N);
            end
            [optimization_vector,fval,exiflag,output] = fmincon(@(x)GradientFunction7v2(x, kx, kf, Myy, frequency_offset, x_positions),...
                initial_optimization_vector,[],[],[],[],lbound,ubound,[],options);        
            if(HighConsoleOutput)
                output
            end            
            initial_optimization_vector = optimization_vector; % Ready for next optimization
            m_opt = zeros(Mspecies, N);
            for jj = 1:Mspecies
                initial_index = ((jj - 1)*2*N);
                m_opt(jj,:) = optimization_vector((initial_index + 1):(initial_index + N)) + ...
                        1i*optimization_vector((initial_index + N + 1):(initial_index + 2*N));
            end
            phi_x_opt = optimization_vector(2*Mspecies*N+(1:N)).';
            R2_opt = optimization_vector(2*(Mspecies)*Nsamples + N + (1:N));
            if(WithGraphics)
                figure(5)
                Plot1DObject(m_init,m_opt,phi_init.',phi_x_opt.');
                pause(1);
            end
            if(output.iterations < GradOniterJ)
                %disp('Warning: less gradient-on iterations than requested');
                %break;
            end
%             if(SaveIntermediate == 1 && ~UseParallel)
%                 save('LastOptimS5', 'm_total', 'phi_total', 'R2_total', 'vobj'); 
%             end
        end
            %save('CurrentOptimizationState','m_opt','phi_x_opt');
    end;    
    
    % WARNING: ONLY WORKS WITH M == 2!!!
    m_total_water(yy,:) = m_opt(1,:);
    m_total_fat(yy,:) = m_opt(2,:);
    
%     for ii = 1:Mspecies
%         m_total(yy,:,ii) = m_opt(ii,:);
%     end
    phi_total(yy,:) = reshape(phi_x_opt,[1, Nsamples]);
    R2_total(yy,:) = reshape(R2_opt,[1, Nsamples]);
    
    %m(yy+1,:,:) = m_total(yy,:,:);
    %phi_aux(yy+1,:) = phi_total(yy,:); 
    %pause
    
    vobj(yy) = fval; % Store the objective value
    fprintf('Done for position y(%d)\n', yy);
end

t_total = toc; % Total time of FIRST

m_total(:,:,1) = m_total_water;
m_total(:,:,2) = m_total_fat;

% Printing
fprintf('Done for positions %i to %i\n', range(1), range(end));
disp('Saving optimization results');
save('LastOptim', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');

fprintf('\n');
disp('Processing time performance:');
fprintf('IDEAL: %5.2f\n', t_ideal);
fprintf('Optim: %5.2f\n', t_total);
fprintf('Ratio: %5.2f\n', t_total/t_ideal);
if(UseParallel)
    matlabpool close
end

if(WithGraphics)
    Plot2DObject(m_total,phi_total,'Optimization');
end
