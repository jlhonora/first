function [M] = get_acquisition(m, x_pos, y_pos, delta_f, phi_xy, kx, ky, kf)

Mspecies = size(m,3);
Nx = length(kx);
Ny = length(ky);

M = zeros(Ny,Nx);
Ma = M;
Maa = M;

delta_f_error = rand(1,Mspecies);
delta_f_error = zeros(1,Mspecies);
delta_f = delta_f + delta_f_error;

kf = kf.';
%kf = zeros(size(kf));
kx = kx.';
ky = ky.';
%phi_x = phi_x.';
if(size(x_pos,1)<size(x_pos,2))
   x_pos = x_pos';
end
if(size(y_pos,1)<size(y_pos,2))
   y_pos = y_pos';
end

FOV = x_pos(end) - x_pos(1);

% for yy = 1:Ny 
% for xx = 1:Nx 
%     for jj = 1:Mspecies
%         current_species = m(yy,xx,jj);
%         Mjj = current_species*exp(-2i*pi*(-(delta_f(jj) + phi_x(yy,xx))*kf(yy) + x_pos(xx)*kx(xx) + y_pos(yy)*ky(yy)));    
%         M(yy,xx) = M(yy,xx) + Mjj;     
%     end    
% end
% end

% for jj = 1:Mspecies
%     current_species = m(:,:,jj);
%     Mjj = exp(-2i*pi*y_pos*ky)*current_species*exp(-2i*pi*(-(delta_f(jj) + phi_xy)*repmat(kf,Ny,1) + x_pos*kx));      
%     M = M + Mjj;     
% end    

% Ma = M;

% M = zeros(Ny,Nx); 
% for jj = 1:Mspecies
%     current_species = m(:,:,jj);
%     for kyy = 1:Ny
%         for kxx = 1:Nx
%             Mjj = exp(-2i*pi*y_pos'*ky(kyy))*current_species*...
%                 exp(-2i*pi*(-(delta_f(jj) + phi_xy)*repmat(kf(kxx),Ny,1) + x_pos*kx(kxx)));      
%             M(kyy,kxx) = M(kyy,kxx) + Mjj;
%         end
%     end
% end    

% 
% load('Mexacto');
% %n = norm(M-Mexacto(:,:,1))
% n2 = norm(Ma - M)
% error('quit');

% M = zeros(Ny,Nx);
% for jj = 1:Mspecies
%     current_species = m(:,:,jj);
%     Mjj = exp(-2i*pi*y_pos*ky)*current_species*exp(-2i*pi*(-(delta_f(jj) + phi_xy)*repmat(kf,Ny,1) + x_pos*kx));      
%     M = M + Mjj;
% end    

% M = zeros(Ny,Nx); %%% Standard %%%
% for jj = 1:Mspecies    
%     for yy = 1:Ny
%         y_val = y_pos(yy);        
%         for xx = 1:Nx
%             current_species = m(yy,xx,jj);
%             x_val = x_pos(xx);
%             phi_val = phi_xy(yy,xx);
%             for kyy = 1:Ny
%                 for kxx = 1:Nx                    
%                     Mjj = current_species*exp(-2i*pi*(-(delta_f(jj) + ...
%                         phi_val)*kf(kxx)' + x_val*kx(kxx) + y_val*ky(kyy)));  
%                     M(kyy,kxx) = M(kyy,kxx) + Mjj;
%                 end
%             end
%         end
%     end
% end

% M = zeros(Ny,Nx);  %%% OK %%%
% for jj = 1:Mspecies
%     current_species = m(:,:,jj);
%     for yy = 1:Ny
%         for xx = 1:Nx
%             Mjj = current_species(yy,xx)*exp(-2i*pi*(-(delta_f(jj) + phi_xy(yy,xx))*...
%                 repmat(kf,Ny,1) + repmat(x_pos(xx)*kx,Ny,1) + repmat(y_pos(yy)*ky',1,Nx)));      
%             M = M + Mjj;
%         end
%     end
% end  

%Ma = M;

% M = zeros(Ny,Nx); 
% for jj = 1:Mspecies
%     current_species = m(:,:,jj);
%     for kyy = 1:Ny
%         for kxx = 1:Nx
%             Mjj = current_species.*exp(-2i*pi*(...
%                 -(delta_f(jj) + phi_xy).*repmat(kf(kxx),Ny,Nx)...
%                 + repmat(x_pos'*kx(kxx),Ny,1)...
%                 + repmat(y_pos*ky(kyy),1,Nx)));
%             M(kyy,kxx) = M(kyy,kxx) + sum(sum(Mjj));
%         end
%     end
% end
% 
% Ma = M;

M = zeros(Ny,Nx); %%% OK %%%  - 25% Faster than above
for jj = 1:Mspecies
    for yy = 1:Ny    
    current_species = m(yy,:,jj);
    for kyy = 1:Ny        
        Mjj = current_species*exp(-2i*pi*(-(delta_f(jj) + ...
            phi_xy(yy,:)')*kf + x_pos*kx + y_pos(yy)*ky(kyy)));  
        M(kyy,:) = M(kyy,:) + Mjj;
    end
    end
end

% load('Mexacto');
%n = norm(Ma-M)
% n = norm(M-Mexacto(:,:,2))
% n = norm(M-Mexacto(:,:,3))
%error('quitting');
