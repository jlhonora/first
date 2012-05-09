% Prints data of an acquisition.
% Shows it ready to copy and paste to the "info" file needed en each
% acquisition folder. This is essential for the Get2DAcquisition function

% Acquired TEs
TE = [4.6 4.8 5 5.4 5.8 6.2 6.9 7.5 8 11 14];

% Indexes of dicom images (IM_0001 ...)
ImgIndex = [249:6:309];

% Raw index of each acquisition
RawIndex = [0:numel(TE)];
Nimages = numel(TE);

for jj = 1:Nimages
    ii = ImgIndex(jj); disp([sprintf('%2.2f ',TE(jj)) 'IM_' sprintf('%04d',ii) ' ' 'IM_' sprintf('%04d',ii + 1) ' ' 'IM_' sprintf('%04d',ii + 2) ' ' 'raw_' sprintf('%03d', RawIndex(jj))])
end