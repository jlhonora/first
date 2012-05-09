%%%% Testeo field maps u otros

resta1 = field_inhomogeneity - phi_total;

resta1 = resta1.*(abs(m_total(:,:,1)) + abs(m_total(:,:,2))>0.5);

resta1b = field_inhomogeneity - phi_ideal;

resta1b = resta1b.*(abs(m_ideal(:,:,1)) + abs(m_ideal(:,:,2))>15);

I = resta1;
Ir = imresize(I,4, 'nearest');

figure(9)
imshow(Ir,[]); colorbar
SetWhiteBackground

Ib = resta1b;
Ibr = imresize(Ib,4, 'nearest');

figure(10)
imshow(Ibr,[]); colorbar
SetWhiteBackground

resta2 = m - m_total;

I2a = resta2(:,:,1);
I2b = resta2(:,:,2);

I2a = imresize(I2a,4,'nearest');
I2b = imresize(I2b,4,'nearest');

resta3 = m - m_ideal;

I3a = resta3(:,:,1);
I3b = resta3(:,:,2);

vals = [resta2 resta3];
ext = max(abs(vals(:)))

I3a = imresize(I3a,4,'nearest');
I3b = imresize(I3b,4,'nearest');

figure(11)
imshow(I2a,[-ext ext]), colorbar
SetWhiteBackground

figure(12)
imshow(I2b,[-ext ext]), colorbar
SetWhiteBackground

figure(13)
imshow(I3a,[-ext ext]), colorbar
SetWhiteBackground

figure(14)
imshow(I3b,[-ext ext]), colorbar
SetWhiteBackground