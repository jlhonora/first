 function [m_ideal phi_ideal n_iter_map] = IDEALsep(m_fft, TE, frequency_offset, varargin)
% Calculates IDEAL separation with unknown field map. 
%
% WARNING: The script IDEALT2sep should be used instead of this outdated
% version.
%
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
%       [m, p] = IDEALsep(m_fft, TE, frequency_offset,...
%                           MyKnownFieldMap,...
%                           'SmoothFieldMap',...
%                           sigma);
%
%       If 'sigma' is not specified, default value 1.5 is used. 
% 
Mspecies = length(frequency_offset);
Nacq = length(TE);

Nrows = size(m_fft,1);
Ncols = size(m_fft,2);
Npoints = Nrows*Ncols;

m_fft_aux = zeros(Nacq, Npoints);
for ii = 1:Nacq
    m_fft_aux_mat = m_fft(:,:,ii);
    m_fft_aux(ii,:) = m_fft_aux_mat(:);
end
m_fft = m_fft_aux;

[SmoothFieldMap, SigmaSmoothFieldMap, KnownFieldMap, phi_0] = GetOptions(varargin, Npoints);
KnownFieldMap = false;

if(KnownFieldMap)
    fprintf('\tField map is known, estimating species directly\n');
end

phi_0 = reshape(phi_0, 1, Npoints);

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

max_iter = 200;
min_delta_phi = 0.1;
max_phi = max(abs(frequency_offset))/2;

m_ideal = zeros(Mspecies, Npoints);
phi_ideal = zeros(size(phi_0));

signum = -1;
signum_delta = 1;

n_iter_map = zeros(size(phi_0));

fprintf('\tComputing IDEAL species separation with unknown field map\n');
for ii = 1:Npoints
    nValues = m_fft(:,ii);
    current_phi = phi_0(ii);
%    current_phi = 0;
%     try
%         if(abs(phi_ideal(ii-1))~=max_phi)
%             current_phi = phi_ideal(ii-1);
%         else
%             current_phi = 0;
%         end
%     catch
%         current_phi = phi_0(ii);
%     end
    n_iter = 0;
    delta_phi = 10;

%     if(ii == Npoints/2)
%         a = 0;
%     end
    
    if(~KnownFieldMap)
        while(abs(delta_phi)>min_delta_phi && n_iter < max_iter)

            % If the object is too small, field inhom. = 0
            if(max(abs(nValues))<(m_max*1e-2))
                delta_phi = 0;
                break;
            end
            
            S_es = nValues.*exp(signum*(2i*pi*current_phi).*TE); %!!! Field. Inhom. sign changed !!!%
            S_est = [real(S_es); imag(S_es);];
            %rho_est = inv(A'*A)*A'*S_est;
            rho_est = A\S_est; %[Re(rho_1) Im(rho_1) ... Re(rho_M) Im(rho_M)]'

            gjn =  (2*pi).*TE.*(-djn*rho_est(1:2:end) -cjn*rho_est(2:2:end)) +...
                   (2i*pi).*TE.*(cjn*rho_est(1:2:end) -djn*rho_est(2:2:end));
            B(1:Nacq,1) = real(gjn);
            B(Nacq + (1:Nacq),1) = imag(gjn);

            S_est2 = zeros(2*Nacq,1);
            S_est2(1:Nacq) = S_est(1:Nacq) - (cjn*rho_est(1:2:end) -djn*rho_est(2:2:end));
            S_est2(Nacq + (1:Nacq)) = S_est(Nacq + (1:Nacq)) - (djn*rho_est(1:2:end) +cjn*rho_est(2:2:end));

            %y_est = inv(B'*B)*B'*S_est2; %[dphi dRe(rho_1) dIm(rho_1) ... dRe(rho_M) dIm(rho_M)]'
            y_est = B\S_est2;
            delta_phi = y_est(1);
            %delta_phi = 0;
            %current_phi = field_inhomogeneity(ii);
            current_phi = current_phi + signum_delta*delta_phi;
            
%             num = exp(signum*(2i*pi*current_phi).*TE(3));
%             current_phi = atan2(imag(num), real(num))/(signum*(2*pi).*TE(3));
            n_iter = n_iter + 1;
            %pause
            
            if(abs(current_phi)>max_phi)
                current_phi = max_phi*sign(current_phi);
            end
            
        end
        n_iter_map(ii) = n_iter; 
    end
    
%     if(n_iter == max_iter)
%         warning(sprintf('\tMax. iterations reached\n'))
%     end
    
    % Simple phase unwrapping
%    if(abs(current_phi)>max(abs(frequency_offset))) 
%           index = abs(frequency_offset)==max(abs(frequency_offset));  
%           max_f = abs(frequency_offset(index));
%           current_phi(current_phi>=max_f) = mod(current_phi(current_phi>=max_f), max_f);
%           current_phi(current_phi<=-max_f) = mod(current_phi(current_phi<=-max_f), -max_f);
%    end
    
    S_es = nValues.*exp(signum*(2i*pi*current_phi).*TE); %!!! Field. Inhom. sign changed !!!%
    S_est = [real(S_es); imag(S_es);];
    %rho_est = inv(A'*A)*A'*S_est; %[Re(rho_1) Im(rho_1) ... Re(rho_M) Im(rho_M)]'
    rho_est = A\S_est;
    m_ideal(:,ii) = rho_est(1:2:end) + 1i*rho_est(2:2:end);
    phi_ideal(ii) = current_phi;
    %m_ideal(:,ii)
    %m(:,ii)
end

if(SmoothFieldMap) % Re-estimate values considering a field map estimation
    sigma = SigmaSmoothFieldMap; % For smoothing field map
    fprintf('\tField Map smoothed with sigma = %2.4f\n', sigma);
    h_filter = fspecial('gaussian', [3 3], sigma);
    phi_ideal_filt = reshape(phi_ideal, Nrows, Ncols);
    phi_ideal_filt = imfilter(phi_ideal_filt, h_filter);
    phi_ideal_filt = reshape(phi_ideal_filt,1,Npoints);
    phi_ideal = phi_ideal_filt;
    for ii = 1:Npoints
        nValues = m_fft(:,ii);
        current_phi = phi_ideal(ii);
        S_es = nValues.*exp(signum*(2i*pi*current_phi).*TE); %!!! Field. Inhom. sign changed !!!%
        S_est = [real(S_es); imag(S_es);];
        %rho_est = inv(A'*A)*A'*S_est; %[Re(rho_1) Im(rho_1) ... Re(rho_M) Im(rho_M)]'
        rho_est = A\S_est;
        m_ideal(:,ii) = rho_est(1:2:end) + 1i*rho_est(2:2:end);
    end
end
    
phi_ideal = phi_ideal';

m = zeros(Nrows,Ncols,Mspecies);

for ii = 1:Mspecies
    m(:,:,ii) = reshape(m_ideal(ii,:),Nrows,Ncols);
end

m_ideal = m;
phi_ideal = reshape(phi_ideal, Nrows, Ncols);
n_iter_map = reshape(n_iter_map, Nrows, Ncols);

fprintf('\tDone\n');
 end

 
function [SmoothFieldMap, SigmaSmoothFieldMap, KnownFieldMap, phi_0] = GetOptions(vars, Npoints)

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
    case 0,
    case 1,
        cl = class(varargin{1});
        switch(cl)
            case 'char',
                if(strcmpi(varargin{1},'smoothfieldmap'))
                    SmoothFieldMap = true;
                end                    
            case {'single','double'},
                phi_0 = double(varargin{1});
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
            case {'single','double'},
                phi_0 = varargin{1};
                KnownFieldMap = true;
                if(strcmp(lower(varargin{2}),'smoothfieldmap'))
                    SmoothFieldMap = true;
                    try
                        SigmaSmoothFieldMap = varargin{3};
                    catch
                    end
                end    
        end
    otherwise,
        warning('Unrecognized number of parameters in IDEAL separation, please re-run.')
end

end
