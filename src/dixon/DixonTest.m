% Dixon test

InPhaseRe = 'C:\Users\JLH\Desktop\ImagenesDixon\InPhase\00010102';
InPhaseIm = 'C:\Users\JLH\Desktop\ImagenesDixon\InPhase\00010170';
OutPhaseRe = 'C:\Users\JLH\Desktop\ImagenesDixon\OutPhase\00010101';
OutPhaseIm = 'C:\Users\JLH\Desktop\ImagenesDixon\OutPhase\00010169';

InPhaseRe = double(dicomread(InPhaseRe));
InPhaseIm = double(dicomread(InPhaseIm));
OutPhaseRe = double(dicomread(OutPhaseRe));
OutPhaseIm = double(dicomread(OutPhaseIm));

InPhase = InPhaseRe + 1i*InPhaseIm;
OutPhase = OutPhaseRe + 1i*OutPhaseIm;

[water, fat, phi] = ExtTwoPointDixon2(InPhase, OutPhase);

water = water - min(min(water));
fat = fat - min(min(fat));

figure(1), imshow(water,[0, 800]); title('Water'), colorbar
figure(2), imshow(fat,[0, 800]); title('Fat'), colorbar
figure(3), imshow(phi,[]); title('Phi'), colorbar

figure(4), imshow(fat./(fat+water),[]); title('Fraction'), colorbar

