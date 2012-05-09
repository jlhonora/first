function [water, fat, phi] = ExtTwoPointDixon2(InPhase, OutPhase)

Nrows = size(InPhase,1);
Ncols = size(OutPhase,2);

Nsamples = size(InPhase(:));

m = zeros(2, Nsamples);

M0 = InPhase(:);
M1 = OutPhase(:);

M1p = M1.*conj(M0)./abs(M0);
M1p2 = M1p.^2;

phi = 0.5.*angle(M1p2);

pM1 = M1p.*exp(-1i.*phi);

m(1,:) = 0.5.*(abs(M0) + pM1);
m(2,:) = 0.5.*(abs(M0) - pM1);

Mspecies = 2;
m_dixon = zeros(Nrows,Ncols,Mspecies);

for ii = 1:Mspecies
    m_dixon(:,:,ii) = reshape(m(ii,:),Nrows,Ncols);
end

phi = abs(reshape(phi, Nrows, Ncols));
fat = abs(m_dixon(:,:,1));
water = abs(m_dixon(:,:,2));


