function [kx ky kf t g] = get_trajectory(gradient_parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Trajectory


traj_type = gradient_parameters(1);
dt = gradient_parameters(2); % time step
A = gradient_parameters(3); % Initial delay
b = gradient_parameters(4); % negative lobe duration
G = gradient_parameters(5); % field strength (T/m)
Nx = gradient_parameters(6); % number of samples in x direction
Ny = gradient_parameters(7); % number of samples in y direction

global gammabar;

	switch traj_type
		
		case 1, 			
            % Delay with zeroes
			t = (dt:dt:A) - dt;
			g = zeros(size(t));
			kx = g; kf = g;
			% Half negative
			taux = A + (1:Nx/2)*dt;
			[gaux kxaux kfaux] = create_traj(Nx/2,taux,dt,-G,0,0,gammabar);
			t = [t taux]; g = [g gaux]; kx = [kx kxaux]; kf = [kf kfaux];
			% Full positive
			taux = A + b + (1:Nx)*dt;
			[gaux kxaux kfaux] = create_traj(Nx,taux,dt,G,kx(end),kf(end),gammabar);
			t = [t taux]; g = [g gaux]; kx = [kx kxaux]; kf = [kf kfaux];
            kf = t;
	end;
    t_c = t;
    % Sample the trajectory
    sampix = zeros([0 0]); % indeces to whatever the sampling is

    %l=0*b; L=1*b; sampix = [sampix find((t_c>(A+l))&(t_c<=(A+L)))]; 

    l=1*b; L=2*b; sampix = [sampix find((t_c>(A+l))&(t_c<=(A+L)))]; 
    l=2*b; L=3*b; sampix = [sampix find((t_c>(A+l))&(t_c<=(A+L)))]; 

    %l=3*b; L=4*b; sampix = [sampix find((t_c>(A+l))&(t_c<=(A+L)))]; 
    %l=4*b; L=5*b; sampix = [sampix find((t_c>(A+l))&(t_c<=(A+L)))]; 

    %l=5*b; L=6*b; sampix = [sampix find((t_c>(A+l))&(t_c<=(A+L)))]; 
    %l=6*b; L=7*b; sampix = [sampix find((t_c>(A+l))&(t_c<=(A+L)))];

    %GLOBAL
    kx = kx(sampix).'; % k-space will be a column
    kf = kf(sampix).';

    if(Ny<=Nx)
        ky = kx(round(1:(Nx/Ny):Nx));
    else
        ky = interp1(1:Nx, kx, linspace(1,Nx,Ny), 'linear', 'extrap');
    end    
    
    kx = kx - kx(round(Nx/2) + 1); % Center in N/2 + 1;
    ky = ky - ky(round(Ny/2) + 1); % Center in N/2 + 1;
    ky = kx;
    kf = kf - (kf(2) - kf(1)); % Center in N/2 + 1;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sub function for trajectories
function [g kx kf] = create_traj(N,t,dt,Gamp,kxstart,kfstart,gammabar)
%
g = Gamp*ones(size(t));
kx = kxstart + (gammabar).*cumsum(g)*dt;
kf = kfstart + (1:N)*dt.*(gammabar);
%kx = kxstart + gamma*cumsum(g)*dt;
%kf = kfstart + (1:N)*dt;