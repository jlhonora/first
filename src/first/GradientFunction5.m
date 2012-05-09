function [f g] = GradientFunction5(x, kx, kf, M, frequency_offset, x_positions)
% Calculates objective function as
%
% abs( M - sum_j sum_p m_j  exp(-i2pi( x*kx + (delta_f + phi_x)*kf ))
%
% No R2 estimation in this function.

if nargout == 1
    f = gf5(x,kx,kf,M,frequency_offset, x_positions);
else
    [f,g] = gf5(x,kx,kf,M,frequency_offset, x_positions);
end
return

Mspecies = numel(frequency_offset);
Nsamples = numel(x_positions);
N = Nsamples;

Nacq = size(kx,1);
m = zeros(Mspecies, Nsamples);

phi_x = x(2*Mspecies*Nsamples+(1:N));
%sign_phi = 1;

%Mestimate = zeros(Nacq, Nsamples);
Mestimate = zeros(1,Nacq*Nsamples);

expon2 = zeros(N,N*Nacq,Mspecies);

midospi = -2i*pi;

for jj = 1:Mspecies
    m(jj,:) = x((jj-1)*2*N + (1:N)) + ...
                1i*x((jj-1)*2*N + N + (1:N));
    expon2(:,:,jj) = exp(midospi.*(-kf(:)*(frequency_offset(jj) + phi_x)...
                                  + kx(:)*(x_positions))).';        
    Mestimate =  Mestimate + (m(jj,:)*expon2(:,:,jj));   
end

% f = sum(sum(abs(fft(M - Mestimate)).^2));
f = sum((abs((M(:) - Mestimate(:))).^2));

if nargout > 1 % gradient required
    g = zeros(size(x'));
    diff_m = (M(:) - Mestimate(:));
    for jj = 1:Mspecies   
        deriv_x   = expon2(:,:,jj);
        deriv_y   = 1i*deriv_x;
        deriv_phi = (-midospi*(kf(:)*m(jj,:)).').*deriv_x;
%         g((1:2*N)+2*N*(jj-1)) = g((1:2*N)+2*N*(jj-1)) - [
%             (conj(deriv_x)*diff_m+deriv_x*conj(diff_m));
%             (conj(deriv_y)*diff_m+deriv_y*conj(diff_m))];
        g((1:2*N)+2*N*(jj-1)) = g((1:2*N)+2*N*(jj-1)) - [
            (diff_m.'*deriv_x'+diff_m'*deriv_x.'),(diff_m.'*deriv_y'+diff_m'*deriv_y.')].';
        g((1:N)+2*N*(Mspecies)) = g((1:N)+2*N*(Mspecies)) - ...
            (diff_m.'*deriv_phi'+diff_m'*deriv_phi.').';
    end
    %g = 2*g;
    %g((1:N)+2*N*(Mspecies)) = real(g((1:N)+2*N*(Mspecies))) + imag(g((1:N)+2*N*(Mspecies)));
end
     

%a=1;


%     gx = zeros(size(p));
%     gy = gx;
%     gv = gx;
%     
%     if use(1),
%         
%         cs0 = cos(2*pi*(v*k0 + p*q/N));
%         sn0 = sin(2*pi*(v*k0 + p*q/N));
%         rem0 = repmat(real(M0 - Me0),[N 1]);
%         imm0 = repmat(imag(M0 - Me0),[N 1]);
%         
%         gx = gx - 2*sum(rem0.*cs0,2) + 2*sum(imm0.*sn0,2);
%         gy = gy - 2*sum(rem0.*sn0,2) - 2*sum(imm0.*cs0,2);
%         gv = gv - 4*pi*sum(rem0.*((imag(m)*k0).*cs0 - (real(m)*k0).*sn0),2)...
%                 + 4*pi*sum(imm0.*((imag(m)*k0).*sn0 + (real(m)*k0).*cs0),2);
%             
%     end;
% 
%     if use(2),
%         
%         cs1 = cos(2*pi*(v*k1 + p*q/N));
%         sn1 = sin(2*pi*(v*k1 + p*q/N));      
%         rem1 = repmat(real(M1 - Me1),[N 1]);
%         imm1 = repmat(imag(M1 - Me1),[N 1]);
%         
%         gx = gx - 2*sum(rem1.*cs1,2) + 2*sum(imm1.*sn1,2);
%         gy = gy - 2*sum(rem1.*sn1,2) - 2*sum(imm1.*cs1,2);
%         gv = gv - 4*pi*sum(rem1.*((imag(m)*k1).*cs1 - (real(m)*k1).*sn1),2)...
%                 + 4*pi*sum(imm1.*((imag(m)*k1).*sn1 + (real(m)*k1).*cs1),2);
%             
%     end;
%     
%     if use(3),
%         
%         cs2 = cos(2*pi*(v*k2 + p*q/N));
%         sn2 = sin(2*pi*(v*k2 + p*q/N));      
%         rem2 = repmat(real(M2 - Me2),[N 1]);
%         imm2 = repmat(imag(M2 - Me2),[N 1]);
%         
%         gx = gx - 2*sum(rem2.*cs2,2) + 2*sum(imm2.*sn2,2);
%         gy = gy - 2*sum(rem2.*sn2,2) - 2*sum(imm2.*cs2,2);
%         gv = gv - 4*pi*sum(rem2.*((imag(m)*k2).*cs2 - (real(m)*k2).*sn2),2)...
%                 + 4*pi*sum(imm2.*((imag(m)*k2).*sn2 + (real(m)*k2).*cs2),2);
%             
%     end;
%     
%     g = [gx;gy;gv];


