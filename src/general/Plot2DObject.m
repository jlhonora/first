function [ObjectsCell] = Plot2DObject(m, phi, varargin)
% Plots the magnitude value of a species-separated object, 
% in tiles of 1 x Mspecies (third dimension of 'm'). Arguments: 
%
%       m   : Stack of images where each pile is a different species.
%       phi : Field Map.
%
%       Additional options:
%           -  Method (str) : Shows string in the title as the used method.
%           - 'NoFieldMap'  : Does not plot the field map.
%           - 'Colorbar'    : Draws an intensity colorbar for each subimage.
%           - 'ShowAll'     : Shows Real and Imaginary parts too
%       

PlotReal = false;
PlotImag = false;
PlotMag = true;

Mspecies = size(m, 3);

Nrows = size(m,1);
Ncols = size(m,2);

MethodLabel = '';
SpeciesLabel = 'Species';

if(~isempty(phi))
    ShowFieldMap = true;
else
    ShowFieldMap = false;
end
ShowColorbar = false;

ObjectsCell = cell(1,7);

Ranges = []; % Image Ranges
cell_notstring = cellfun('isclass', varargin, 'double');
if(sum(cell_notstring)~=0)
   Ranges = varargin{cell_notstring};
   varargin(cell_notstring) = [];
end

% Possibilities for variable arguments
ArginPosib = {'NoFieldMap'; 'Colorbar'; 'ShowAll'; 'Ranges'}; 
n_argin_posib = length(ArginPosib);
n_argin = length(varargin);
vec_argin = 1:n_argin;

ArginPosib = repmat(lower(ArginPosib),1,n_argin);
Rvarargin = repmat(lower(varargin),n_argin_posib,1);

[ArginMatchesP, ArginMatchesR] = find(strcmpi(ArginPosib, Rvarargin));

n_matches = length(ArginMatchesP);

if(n_argin>0)
    for ii = 1:n_matches
        switch(ArginMatchesP(ii))
            case 1,
                ShowFieldMap = false;
            case 2,
                ShowColorbar = true;
            case 3,
                PlotReal = true;
                PlotImag = true;
            case 4,
                fprintf('\tUsing user-supplied ranges for images\n');
        end
    end
    vec_argin(ArginMatchesR) = [];    
    n_argin = length(vec_argin);    
    switch(n_argin)
        case 0,
        case 1, 
            MethodLabel = [varargin{vec_argin(1)} ' - '];
        case 2,
            MethodLabel = [varargin{vec_argin(1)} ' - '];
            SpeciesLabel = varargin{vec_argin(2)};
        otherwise,
            MethodLabel = [varargin{vec_argin(1)} ' - '];
            SpeciesLabel = varargin{vec_argin(2)};
            warning('Unrecognized species label at PlotObject, using second string');
    end
end

TotalPlots = PlotReal + PlotImag + PlotMag;

try
    global x_positions;
catch
    warning('Object positions variable x_positions not found.');
    x_positions = 1:N;
end
    
try
    global y_positions;
catch
    warning('Object positions variable y_positions not found.');
    y_positions = 1:N;
end

RangePhi = [];

if(isempty(Ranges))
    abs_intensity = [min(min(min(abs(m)))) max(max(max(abs(m))))];
    real_intensity = [min(min(min(real(m)))) max(max(max(real(m))))];
    imag_intensity = [min(min(min(imag(m)))) max(max(max(imag(m))))];
else
    abs_intensity = [Ranges(3,1) Ranges(3,2)];
    real_intensity = [Ranges(1,1) Ranges(1,2)];
    imag_intensity = [Ranges(2,1) Ranges(2,2)];
    if(size(Ranges,1)>3)
        RangePhi = Ranges(4,:);
    end
end

x_pos = x_positions;
y_pos = y_positions;

x_min = min(x_positions);
x_max = max(x_positions);

y_min = min(y_positions);
y_max = max(y_positions);

current_fig = gcf;
if(current_fig==1)
    figure(1)
else
    figure(current_fig + 1);
end

clf

for ii = 1:Mspecies
    ObjectsCell{ii} = real(m(:,:,ii));
    ObjectsCell{ii + Mspecies} = imag(m(:,:,ii));
    ObjectsCell{ii + 2*Mspecies} = abs(m(:,:,ii));
    if(PlotReal)        
        subplot(TotalPlots,Mspecies,ii);
        imshow(ObjectsCell{ii}, real_intensity);
        if(ShowColorbar)
            colorbar;
        end    
        if(ii==1)
            ylabel('Real');
        end
    end
         
%    axis([x_min x_max y_min y_max]);       
    if(PlotImag)
        subplot(TotalPlots,Mspecies,ii + (PlotImag)*Mspecies);    
        imshow(ObjectsCell{ii + Mspecies}, imag_intensity);
        if(ShowColorbar)
            colorbar;
        end
        if(ii==1)
            ylabel('Imag');
        end
    end
%    axis([x_min x_max y_min y_max]);    
    if(PlotMag)
        subplot(TotalPlots,Mspecies,ii + (PlotImag + PlotReal)*Mspecies);    
        imshow(ObjectsCell{ii + 2*Mspecies}, abs_intensity);
        if(ShowColorbar)
            colorbar;
        end
        if(ii==1)
            ylabel('Abs');
        end
    end
    subplot(TotalPlots,Mspecies,ii);
    title(sprintf('%s%s #%d',MethodLabel, SpeciesLabel, ii));
end

set(gcf, 'Color', [1 1 1]);
truesize(gcf, [Nrows, Ncols]);

ObjectsCell{3*Mspecies + 1} = phi;
if ShowFieldMap
    figure(gcf + 1)
    set(gcf, 'Color', [1 1 1]);
    imshow(phi, RangePhi); truesize(gcf, [Nrows, Ncols]);
    colorbar;
    title(sprintf('%sField Inhomogeneity',MethodLabel));
end

drawnow

