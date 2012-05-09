function R2fit = GetR2fit(AcqStructs, TE)
% Calculates the R2* decay for a stack of images.
% Can be easily modified to work with a stack of matrixes rather than the
% acquisition structs used in the SpeciesSeparation script.

Nacq = length(TE);
[Nrows, Ncols] = size(AcqStructs(1).Mag);

R2fitv = zeros(1,Nrows*Ncols);
R2fit = zeros(Nrows, Ncols);
MagSet = zeros(Nrows, Ncols, Nacq);

for jj = 1:Nacq
   MagSet(:,:,jj) = log(AcqStructs(jj).Mag); 
end

% for nn = 1:(Nrows*Ncols)
%     [ii,jj] = ind2sub([Nrows,Ncols],nn);
%     output = polyfit(TE, reshape(MagSet(ii,jj,:), 1, Nacq), 1);
%     R2fitv(nn) = -output(1);
% end
% R2fit = reshape(R2fitv,Nrows, Ncols);

row_range = [62:67];
col_range = [1:Ncols];

for ii = row_range
    for jj = col_range
        %[ii,jj] = ind2sub([Nrows,Ncols],nn);
        output = polyfit(TE, reshape(MagSet(ii,jj,:), 1, Nacq), 1);
        R2fit(ii,jj) = -output(1);
    end
end

R2fit(R2fit>250) = 0;
R2fit(R2fit<0) = 0;

se = fspecial('gaussian', [5 5] , 1);
R2fit = imfilter(R2fit, se);
