function [m,phi] = Get2DObject(Nrows,Ncols, Mspecies)
% Used for 2D simulation purposes. Reads artificial phantom images. 

%Mspecies = 3;

Iwater = rgb2gray(imread('./Images/water.png'));
Iwater = double(imresize(Iwater, [Nrows Ncols], 'nearest'));

Ifat = rgb2gray(imread('./Images/fat.png'));
Ifat = double(imresize(Ifat, [Nrows Ncols], 'nearest'));

[xx,yy] = meshgrid(linspace(-3,3,Ncols), linspace(-3,3,Nrows));
phase_water = sin(2*(xx+yy));
phase_fat = sin(2*(xx-yy));

%phase_water = zeros(size(xx));
%phase_fat = zeros(size(xx));

m = zeros(Nrows,Ncols,Mspecies);
m(:,:,1) = Iwater.*exp(1i.*phase_water);
m(:,:,2) = Ifat.*exp(1i.*phase_fat);

if(Mspecies==3)
    Isilicone = rgb2gray(imread('./Images/silicone.png'));
    m(:,:,3) = double(imresize(Isilicone, [Nrows Ncols], 'nearest'));
end

phi = zeros(Nrows,Ncols);

phi = peaks(xx,yy);
phi_amplitude = 10;
phi = phi./(max(max(abs(phi)))/phi_amplitude);
%phi = zeros(Nrows,Ncols);