function Mask = GetMask(I, fraction, minPix)
% Applies a threshold to the input image I defined as a fraction of the
% range of the image. Then all isles of less than minPix are eliminated.

% fraction = 0.015;
% minPix = 20;
% I = AcqStructs(1).Mag;
if(size(I,3)>1)
    I = rgb2gray(I);
end

range = max(max(I)) - min(min(I));
minValue = range*fraction;

Mask = I>minValue;
Mask = imfill(Mask,'holes');

Mask = bwareaopen(Mask, minPix, 8);

se = strel('disk',5);
Mask = imdilate(Mask,se);

