function [] = Plot1DObject(m, m_opt, phi, phi_opt, varargin)
% Plots the magnitude value of a 1D species-separated object, 
% in tiles of 1 x Mspecies (third dimension of 'm'). Arguments: 
%
%       m       : Matrix with each line a different species. Used as the
%               reference value for the species.
%       m_opt   : Same as m, but it is the result of the optimization.
%       phi     : Field map. Used as reference.
%       phi_opt : Field map result of the optimization.
%
%       Additional options:
%           -  Method (str) : Shows string in the title as the used method.
%  

PlotReal = 0;
PlotImag = 0;
PlotMag = true;
TotalPlots = PlotReal + PlotImag + PlotMag;

Mspecies = size(m_opt, 1);

if(length(size(m_opt))>2)
    m_aux = zeros(size(phi_opt));
    Mspecies = size(m_opt,3);    
    for ii = 1:Mspecies;
        m_aux(ii,:) = m_opt(:,:,ii);
    end
    m_opt = m_aux;
end

if(length(size(m))>2)
    m_aux = zeros(size(phi_opt));
    for ii = 1:Mspecies;
        m_aux(ii,:) = m(:,:,ii);
    end
    m = m_aux;
end



MethodLabel = '';

if(length(varargin)>0)
    MethodLabel = [varargin{1} ' - '];           
end

try
    global x_positions;
catch
    warning('Object positions variable x_positions not found.');
    x_positions = 1:N;
end

x_pos = x_positions;


phi_scale = max(max(abs(m_opt)))/max(max([abs(phi_opt) abs(phi)]));
phi_scale = min([phi_scale 20]);

if(phi_scale<1e-4)
    phi_scale = 1e-4;
end

x_min = min(x_positions);
x_max = max(x_positions);

max_range = max(max(abs([m; m_opt;])));

y_min = -max_range*1.1;
y_max = max_range*1.1;


for ii = 1:Mspecies    
    if(PlotReal)       
        subplot(TotalPlots,Mspecies,ii);
        plot(x_pos,real(m(ii,:).'),'b-',...
            x_pos,real(m_opt(ii,:).'),'ro',...
            x_pos,phi'*phi_scale,'k-',...
            x_pos,phi_opt'*phi_scale,'k*'); grid on;
            legend('re(m)', 're(m_{opt})', '\phi', '\phi_{opt}');
        axis([x_min x_max y_min y_max]);                
    end
    if(PlotImag)
        subplot(TotalPlots,Mspecies,ii + 2*PlotReal);
        plot(x_pos,imag(m(ii,:).'),'b-',...
            x_pos,imag(m_opt(ii,:).'),'ro',...
            x_pos,phi'*phi_scale,'k-',...
            x_pos,phi_opt'*phi_scale,'k*'); grid on;
        legend('im(m)', 'im(m_{opt})', '\phi', '\phi_{opt}');
        axis([x_min x_max y_min y_max]);    
    end
    if(PlotMag)
        subplot(TotalPlots,Mspecies,ii + 2*PlotReal + 2*PlotImag);
        plot(x_pos,abs(m(ii,:).'),'b-',...
            x_pos,abs(m_opt(ii,:).'),'ro',...
            x_pos,phi'*phi_scale,'k-',...
            x_pos,phi_opt'*phi_scale,'k*'); grid on;
        legend('|m|', '|m_{opt}|', '\phi', '\phi_{opt}');
        axis([x_min x_max y_min y_max]);
    end
    xlabel(sprintf('x [m] \nField inhom. divided by %3.4g',1/phi_scale));
    subplot(TotalPlots,Mspecies,ii);
    title(sprintf('%sSpecies #%d',MethodLabel, ii));
end

set(gcf, 'Color', [1 1 1]);
drawnow

