function [AcqStructs] = Get2DAcquisition(direc, RequestedTE, varargin)
% Returns file-stored acquisitions at specified TEs
%
%
%

ImagesFromDicom = true;
DicomOnly = false;

if(~isempty(varargin))
    if(strcmpi(varargin{1},'NoDicom'))
        ImagesFromDicom = false;
    elseif(strcmpi(varargin{1},'DicomOnly'))
        DicomOnly = true;
    end
end

if(ImagesFromDicom)
    disp('Images based on Dicom files');
else
    disp('Images based on k-space data');
end

if(DicomOnly)
    disp('K-space data based on images from Dicom files');
end

fprintf('\n');
fprintf('Acq#  R.TE[ms]   TE[ms]      TR[ms]   FA[deg]  BW[Hz/pix]  NSA   Size[RxC] RawSize[RxC]\n');
fprintf('__________________________________________________________________________________________\n');
formatstring = ' %02d   %2.2f   %2.11f   %3.1f     %2.1f     %3.4f    %i    %i x %i  %i x %i\n';

n_req = length(RequestedTE);

InfoStruct = struct('TE', [], 'Mag', [], 'Real', [], 'Imag', [], 'KSpace', []);

AcqStructs = struct('TE', [], 'Mag', [], 'Real', [], 'Imag', [], 'KSpace', [], 'Pos', [], 'KTraj', [], 'DicomInfo', []);
AcqStructs = repmat(AcqStructs,n_req,1);

[TE, mag_name, re_name, im_name, raw_name] = ...
    textread([direc '/info'],'%f %s %s %s %s','headerlines',1);

InfoStruct.TE = TE.*1e-3;
InfoStruct.Mag = mag_name;
InfoStruct.Real = re_name;
InfoStruct.Imag = im_name;
InfoStruct.KSpace = raw_name;

for ii = 1:n_req
    index = find(InfoStruct.TE==RequestedTE(ii));
    if(isempty(index))
        error(sprintf('Requested TE = %2.4f not found',RequestedTE(ii)));
        %continue;
    end
    % Mag
    [I,dicom_info] = GetDicom([direc '/' char(InfoStruct.Mag(index))]);
    RealTE = dicom_info.EchoTime;
    AcqStructs(ii).TE = RealTE;
    AcqStructs(ii).Mag = I;
    AcqStructs(ii).DicomInfo = dicom_info;  
    Nx = size(I,2);
    Ny = size(I,1);
    % Water-Fat shift
    try
        global frequency_offset
        idx = max(abs(frequency_offset))==abs(frequency_offset);
        OffResonanceFreq = abs(frequency_offset(idx));
    catch
        OffResonanceFreq = 217.1427;
    end
    WaterFatShiftPixels = dicom_info.Private_2001_1022; % pixels
    HertzPerPixel = OffResonanceFreq/WaterFatShiftPixels; % Hz/pixels
    MaxFrequency = OffResonanceFreq*Nx/2/WaterFatShiftPixels;
    T2 = 1/(2*MaxFrequency);
    %HertzPerPixel = dicom_info.Private_2005_1033;
    %HertzPerPixel = dicom_info.Private_2005_100c;
    %HertzPerPixel = 100;
    % kx
    Tx = dicom_info.PixelSpacing(2);
    [x_positions, kx] = GetKSpace(Nx,Tx*1e-3);
    % ky    
    Ty = dicom_info.PixelSpacing(1);
    [y_positions, ky] = GetKSpace(Ny,Ty*1e-3);    
    AcqStructs(ii).Pos = struct('x_pos', x_positions, 'y_pos', y_positions);
    % kf
    T = 1/(2*HertzPerPixel*Nx/2); % Sampling time in readout
    kf = RealTE*1e-3 + ((-Nx/2):(Nx/2-1))*T;
    AcqStructs(ii).KTraj = struct('kx', kx, 'ky', ky, 'kf', kf);    
    %Raw
    if(~DicomOnly)
        filename = [direc '/' char(InfoStruct.KSpace(index))];
        [kspace, petable, frx, noi] = raw_READ(filename,'NoConsoleOutput');
        %[raw_info] = raw_read_info(filename, petable);
        kspacep = mean(kspace,3);
        kspacep = kspacep.';
        NrowsIm = size(I,1);
        NcolsIm = size(I,2);
        NrowsK = size(kspacep,1);
        NcolsK = size(kspacep,2);    
        oversampling_x = NcolsK/NcolsIm;
        oversampling_y = NrowsK/NrowsIm;
        
%         frxi = interp1(1:length(frx), frx,...
%             1:1/oversampling_x:length(frx)+1/oversampling_x, 'spline', 'extrap');
%         frxi = repmat(frxi,1, NcolsK);
%        kspacep = kspacep.*repmat(frx,1, NcolsK);

        switch(oversampling_x)
            case 1, 
                Im = fftshift(ifft2(fftshift(kspacep)));
            case 2,
                if(oversampling_y==2)
                    Im = fftshift(ifft2(fftshift((kspacep))));
                    Im = [Im((3*end/4+1):end, (end/4+1):3*end/4); Im(1:end/4, (end/4+1):3*end/4)];
                    Im = -1i*Im;
                    %Im = -real(Im) + 1i*imag(Im);
                else
                    Im = fftshift(ifft2(fftshift(kspacep)),2);                    
                    Im = Im(:,end/4+1:3*end/4);
                    
                    %Im = conj(1i*Im); % No funciona para acq in vivo
                    % 128x256, como en 11-04-11
                    Im = -Im;
                end                
            otherwise
                error('Oversampling factor unknown');
        end
        %Im = Im.';

        kspace = fftshift(fft2(fftshift(Im)));
        AcqStructs(ii).KSpace = kspace;   

        Im = fftshift(ifft2(fftshift(kspace)));
    end
    
    if(ImagesFromDicom)
        %Real
        [I] = GetDicom([direc '\' char(InfoStruct.Real(index))]);
        AcqStructs(ii).Real = I;
        %Imag
        [I] = GetDicom([direc '\' char(InfoStruct.Imag(index))]);
        AcqStructs(ii).Imag = I;    
    else
        AcqStructs(ii).Mag  = abs(Im);
        AcqStructs(ii).Real = real(Im);
        AcqStructs(ii).Imag = imag(Im);
    end
    
    if(DicomOnly) % Reconstruct k-space based on images
        Mdicom = AcqStructs(ii).Real - 1i*AcqStructs(ii).Imag;
        AcqStructs(ii).KSpace = fftshift(fft2(fftshift(Mdicom)));
        mdicom = fftshift(ifft2(fftshift(AcqStructs(ii).KSpace))); 
        AcqStructs(ii).Real = real(mdicom);
        AcqStructs(ii).Imag = imag(mdicom);
        AcqStructs(ii).Mag = abs(mdicom);
    end
    FlipAngle = dicom_info.FlipAngle;
    TR = dicom_info.RepetitionTime;
    SliceThickness = dicom_info.SliceThickness;
    Nx = dicom_info.Width;
    Ny = dicom_info.Height;
    NSA = dicom_info.NumberOfAverages;
    fprintf(sprintf(formatstring, ii, RequestedTE(ii)*1000, RealTE , TR, FlipAngle, HertzPerPixel, NSA, Ny, Nx, length(y_positions), length(x_positions)));
end
fprintf('\n');

function [pos, k] = GetKSpace(N,T)
pos = ((-N/2):(N/2-1))*T;
U = 1/(N*T);
k = ((-N/2):(N/2-1))*U;
