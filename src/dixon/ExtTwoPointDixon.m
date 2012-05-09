function [m_dixon, phi_dixon] = ExtTwoPointDixon(M_acq)

Nacq = size(M_acq,3);
Nrows = size(M_acq,1);
Ncols = size(M_acq,2);

M_acq_aux = zeros(Nacq,Nrows*Ncols);

for ii = 1:Nacq
    Ma = M_acq(:,:,ii);
    M_acq_aux(ii,:) = Ma(:);
end

M_acq = M_acq_aux;

Nsamples = size(M_acq,2);

m = zeros(2, Nsamples);

M0 = M_acq(1,:);
M1 = M_acq(2,:);

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

phi_dixon = reshape(phi, Nrows, Ncols);


