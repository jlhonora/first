 function [m_ideal phi_ideal, R2_map] = IDEALT2sep(m_fft, TE, frequency_offset, varargin)
% Calculates IDEAL separation with unknown field map and T2* estimation. 
% Parameters:
%   
%       m_fft : Nrows x Ncols x Nacquisition complex matrix with intensity
%               values of each acquisition
%       TE    : The echo times for each acquisition
%       frequency_offset : The resonating frequency (Hz) of each species 
%                          (-220 for fat)    
%       
% Full options:
%
%       [m, p] = IDEALT2sep(m_fft, TE, frequency_offset,...
%                           MyKnownFieldMap,...
%                           'SmoothFieldMap',...
%                           sigma,...
%                           R2_starting_point);
%
%       If 'sigma' is not specified, default value 1.5 is used. 
%       R2_starting_point is used as a starting point for T2*-IDEAL

Mspecies = length(frequency_offset);
Is2D = 1;
CalculateR2Map = 1;

if(length(size(m_fft))>2)
    Is2D = true;
    Nacq = length(TE);
    Nrows = size(m_fft,1);
    Ncols = size(m_fft,2);
    Npoints = Nrows*Ncols;
else
    Npoints = size(m_fft,2);    
end

DimType = '1D';
if(Is2D)
    DimType = '2D';
end
fprintf(sprintf(...
    ['\tComputing ' DimType ' IDEAL species separation with unknown field map\n']));

m_fft_aux = zeros(Nacq, Nrows*Ncols);
for ii = 1:Nacq
    m_fft_aux_mat = m_fft(:,:,ii);
    m_fft_aux(ii,:) = m_fft_aux_mat(:);
end
m_fft = m_fft_aux;

[SmoothFieldMap, SigmaSmoothFieldMap, KnownFieldMap, phi_0, R2_map] = GetOptions(varargin, Npoints);
KnownFieldMap = false;

phi_0 = reshape(phi_0, 1, Npoints);
R2_map = reshape(R2_map, 1, Npoints);
R2_map = zeros(size(phi_0));

if(CalculateR2Map)
    if(Mspecies~=2)
        fprintf('\tUnable to calculate R2* map, Mspecies != 2\n');
        CalculateR2map = false;
    else
        fprintf('\tIncluding R2* estimation\n'); 
    end
end

m_max = max(max(abs(m_fft)));

if(size(TE,1)<size(TE,2))
    TE = TE'; % Must be a column vector
end

cjn = cos((2*pi).*(TE*frequency_offset));
djn = sin((2*pi).*(TE*frequency_offset));

A = zeros(2*Nacq, 2*Mspecies);
A(1:Nacq, 1:2:(2*Mspecies)) = cjn;
A(Nacq + (1:Nacq), 2:2:(2*Mspecies)) = cjn;
A(1:Nacq, 2:2:(2*Mspecies)) = -djn;
A(Nacq + (1:Nacq), 1:2:(2*Mspecies)) = djn;

B = zeros(2*Nacq, 2*Mspecies + 1);
B(1:Nacq, 2:2:(2*Mspecies + 1)) = cjn;
B(Nacq + (1:Nacq), 3:2:(2*Mspecies + 1)) = cjn;
B(1:Nacq, 3:2:(2*Mspecies + 1)) = -djn;
B(Nacq + (1:Nacq), 2:2:(2*Mspecies + 1)) = djn;

m_ideal = zeros(Mspecies, Npoints);
phi_ideal = zeros(size(phi_0));

max_iter = 100;
min_delta_phi = 0.01;

signum = -1;

row_range = [63:66];
col_range = [1:Ncols];

point_range = zeros(1, length(row_range)*length(col_range));
count = 1;
for ii = row_range
    for jj = col_range
        point_range(count) = sub2ind([Nrows, Ncols], ii, jj);
        count = count + 1;
    end    
end


for ii = point_range
    nValues = m_fft(:,ii);
    current_phi = phi_0(ii);
    n_iter = 0;
    delta_phi = 1000;
    if(~KnownFieldMap)
        while(abs(delta_phi)>min_delta_phi && n_iter < max_iter)

    %         % If the object is too small, field inhom. = 0
    %         if(max(abs(nValues))<(m_max*1e-2))
    %             delta_phi = 0;
    %             break;
    %         end

            S_es = nValues.*exp(signum*(2i*pi*current_phi).*TE); %!!! Field. Inhom. sign changed !!!%
            S_est = [real(S_es); imag(S_es);];
            rho_est = A\S_est; %[Re(rho_1) Im(rho_1) ... Re(rho_M) Im(rho_M)]'

            gjn =  (2*pi).*TE.*(-djn*rho_est(1:2:end) -cjn*rho_est(2:2:end)) +...
                   (2i*pi).*TE.*(cjn*rho_est(1:2:end) -djn*rho_est(2:2:end));
            B(1:Nacq) = real(gjn);
            B(Nacq + (1:Nacq)) = imag(gjn);

            S_est2 = zeros(2*Nacq,1);
            S_est2(1:Nacq) = S_est(1:Nacq) - (cjn*rho_est(1:2:end) -djn*rho_est(2:2:end));
            S_est2(Nacq + (1:Nacq)) = S_est(Nacq + (1:Nacq)) - (djn*rho_est(1:2:end) +cjn*rho_est(2:2:end));

            %y_est = inv(B'*B)*B'*S_est2; %[dphi dRe(rho_1) dIm(rho_1) ... dRe(rho_M) dIm(rho_M)]'
            y_est = B\S_est2;
            delta_phi = y_est(1);
            %delta_phi = 0;
            %current_phi = field_inhomogeneity(ii);
            current_phi = current_phi - signum*delta_phi;
            n_iter = n_iter + 1;
            %pause
            current_phi(current_phi>100) = 100;
            current_phi(current_phi<-100) = -100;
        end
    end
    
    if(n_iter == max_iter)
        %warning(sprintf('\tMax. iterations reached\n'))
    end
    
    % Simple phase unwrapping
%     if(abs(current_phi)>max(abs(frequency_offset))) 
%          index = abs(frequency_offset)==max(abs(frequency_offset));  
%          current_phi = mod(current_phi, frequency_offset(index));
%     end
    
    S_es = nValues.*exp(signum*(2i*pi*current_phi).*TE); %!!! Field. Inhom. sign changed !!!%
    S_est = [real(S_es); imag(S_es);];
    %rho_est = inv(A'*A)*A'*S_est; %[Re(rho_1) Im(rho_1) ... Re(rho_M) Im(rho_M)]'
    rho_est = A\S_est;
    m_ideal(:,ii) = rho_est(1:2:end) + 1i*rho_est(2:2:end);
    phi_ideal(ii) = current_phi;

    %%%%%%%%%%%%%% T2* Estimation %%%%%%%%%%%%%%
    if(CalculateR2Map)
        [water, fat, field_map, R2_map_i] = GetT2map(m_fft(:,ii), phi_ideal(ii), TE, frequency_offset, R2_map(ii));
        phi_ideal(ii) = field_map;
        R2_map(ii) = R2_map_i;
        if((abs(water) ~= 0) || (abs(fat) ~= 0))
            m_ideal(:, ii) = [water; fat]; 
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

phi_ideal = phi_ideal.';

SmoothR2 = 1;
sigmaR2 = 0.0001;

if(SmoothFieldMap && Is2D) % Re-estimate values considering a field map estimation
    sigma = SigmaSmoothFieldMap; % For smoothing field map
    sigma = 0.0001;
    fprintf(sprintf('\tField Map smoothed with sigma = %2.4f\n', sigma));
    h_filter = fspecial('gaussian', [3 3], sigma);
    phi_ideal_filt = reshape(phi_ideal, Nrows, Ncols);
    phi_ideal_filt = imfilter(phi_ideal_filt, h_filter);
    phi_ideal_filt = reshape(phi_ideal_filt,1,Npoints);
    if(SmoothR2)
        fprintf(sprintf('\tR2* Map smoothed with sigma = %2.4f\n', sigmaR2));
        phi_ideal = phi_ideal_filt;
        h_filter = fspecial('gaussian', [5 5], sigmaR2);
        R2_map_filt = reshape(R2_map, Nrows, Ncols);
        R2_map_filt = imfilter(R2_map_filt, h_filter);
        R2_map_filt = reshape(R2_map_filt,1,Npoints);
        R2_map = R2_map_filt;    
        Delta_f = abs(frequency_offset(1) - frequency_offset(2)); % Frequency deviation in Hz.
        Delta_f = -Delta_f;
        AR2 = ones(Nacq, 2);
        AR2(:,2) = exp(2i*pi*Delta_f*TE);
    end
    
    for ii = point_range
        nValues = m_fft(:,ii);
        current_phi = phi_ideal(ii);
        S_es = nValues.*exp(signum*(2i*pi*current_phi).*TE); %!!! Field. Inhom. sign changed !!!%
        S_est = [real(S_es); imag(S_es);];
        %rho_est = inv(A'*A)*A'*S_est; %[Re(rho_1) Im(rho_1) ... Re(rho_M) Im(rho_M)]'
        rho_est = A\S_est;
        m_ideal(:,ii) = rho_est(1:2:end) + 1i*rho_est(2:2:end);
        
        if(SmoothR2)
            phi_est = (current_phi + 1i/(2*pi)*R2_map(ii));
            %phi_est = -phi_est;
            P = diag(exp(-2i*pi*phi_est*TE));
            rho_est = inv(AR2'*AR2)*AR2'*P*m_fft(:,ii); % [water; fat]
            m_ideal(:,ii) = rho_est;
        end
    end
end  

if(Is2D)
    m = zeros(Nrows,Ncols,Mspecies);
    for ii = 1:Mspecies
        m(:,:,ii) = reshape(m_ideal(ii,:),Nrows,Ncols);
    end
    m_ideal = m;
    phi_ideal = reshape(phi_ideal, Nrows, Ncols);
    R2_map = reshape(R2_map, Nrows, Ncols);    
end

fprintf('\tDone\n');
 end

function [water, fat, phi, R2] = GetT2map(m_fft, phi_ideal, TE, frequency_offset, R2_start)

    curve_fitting = 0;
    
    if(curve_fitting)
        curve = fit(TE, abs(m_fft), 'exp1', 'StartPoint', [abs(m_fft(1)),-1000]);
        t2 = -1/curve.b;
        water = 0;
        fat = 0;
        phi = phi_ideal;
        if(t2<0) 
            t2 = 2;
        elseif(t2(t2>2))
            t2 = 2;
        end
    else
        % Estimates T2 map from the formulation in Yu, 2007.
        % Only for fat and water separation
        max_iter = 100;
        min_delta_phi = 0.01;
        min_delta_R2 = 0.0001;

        n_iter = 0;
        delta_phi = 100;
        delta_R2 = 100;
        max_R2 = 350; % min_t2 = 0.00285; 2.85 miliseconds
        min_R2 = 0; %max_t2 = inf;

        Nacq = size(m_fft,1);

        Delta_f = abs(frequency_offset(1) - frequency_offset(2)); % Frequency deviation in Hz.
        Delta_f = -Delta_f;  

        phi = phi_ideal;
        R2 = R2_start;

        A = ones(Nacq, 2);
        A(:,2) = exp(2i*pi*Delta_f*TE);

        while(((abs(delta_phi) > min_delta_phi) || (abs(delta_R2) > min_delta_R2)) && (n_iter < max_iter))
            phi_est = (phi + 1i/(2*pi)*R2);
            %phi_est = -phi_est;
            P = diag(exp(-2i*pi*phi_est*TE));
            rho_est = inv(A'*A)*A'*P*m_fft;
            B = [(rho_est(1) + rho_est(2)*exp(2i*pi*Delta_f*TE)).*(2i*pi*TE) ones(Nacq,1) exp(2i*pi*Delta_f*TE)];
            delta = inv(B'*B)*B'*(P*m_fft - A*rho_est);
            phi_est = phi_est + delta(1);
            delta_phi = real(delta(1));
            delta_R2 = imag(delta(1))*2*pi;
            phi = real(phi_est);
            R2 = (imag(phi_est)*2*pi);
            n_iter = n_iter + 1;
            R2(R2>max_R2) = max_R2;
            R2(R2<min_R2) = min_R2;
        end

    %         % Simple phase unwrapping
        if(abs(phi)>max(abs(Delta_f))) 
    %           max_f = abs(Delta_f);
    %           phi(phi>=max_f) = mod(phi(phi>=max_f), max_f);
    %           phi(phi<=-max_f) = mod(phi(phi<=-max_f), -max_f);
       end

        %t2 = 1/R2;
        %t2 = R2;
        water = rho_est(1);
        fat = rho_est(2);
    end
end


function [SmoothFieldMap, SigmaSmoothFieldMap, KnownFieldMap, phi_0, R2_map] = GetOptions(vars, Npoints)

varargin = vars;

n_argin = length(varargin);

phi_0 = zeros(1,Npoints);
SmoothFieldMap = false;
KnownFieldMap = false;
SigmaSmoothFieldMap = 1.5; % Default sigma for field smoothing.

switch(n_argin) % Deals with variable input arguments. Possibilities:
                %    1. Known Field map (Nrows x Ncols double matrix)
                %    2. String 'SmoothFieldMap'
                %    3. Extra parameter, 'sigma' for smoothing field map.        
    case 1,
        cl = class(varargin{1});
        switch(cl)
            case 'char',
                if(strcmp(lower(varargin{1}),'smoothfieldmap'))
                    SmoothFieldMap = true;
                end                    
            case 'double',
                phi_0 = varargin{1};
                KnownFieldMap = true;
        end
    case {2,3},
        cl = class(varargin{1});
        switch(cl)
            case 'char',
                if(strcmp(lower(varargin{1}),'smoothfieldmap'))
                    SmoothFieldMap = true;
                    if(length(varargin{2})==1)
                        SigmaSmoothFieldMap = varargin{2};
                    else
                        phi_0 = varargin{2};
                        KnownFieldMap = true;
                    end
                end  
                try
                    phi_0 = varargin{3};
                    KnownFieldMap = true;
                catch
                end
            case {'double', 'single'},
                phi_0 = double(varargin{1});
                KnownFieldMap = true;
                if(strcmp(lower(varargin{2}),'smoothfieldmap'))
                    SmoothFieldMap = true;
                    try
                        %SigmaSmoothFieldMap = varargin{3};
                        R2_map = varargin{3};
                        fprintf('\tUsing R2* fitting from acquisitions as starting point\n');
                    catch
                    end
                end    
        end
        
    otherwise,
        warning('Unrecognized number of parameters in IDEAL separation, please re-run.')
end
end
