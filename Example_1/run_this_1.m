
N = 4;               % Number of parallel chains in MCMC
T = 2000;            % Number of iterations
Nx = 5;              % Dimension of model parameters
Ny = 731;            % Dimension of model responses
t1 = 80;             % After which iteration the Kalman proposal is used
t2 = ceil(0.3*T);    % After which iteration the Kalman proposal is not used
Ne = 100;            % Number of archive samples for the Kalman proposal 

xmin = [1.0 0.10 0.10 0.00 0.10];      % Lower bounds of the parameters
xmax = [500 2.00 0.99 0.10 0.99];      % Upper bounds of the parameters
range = [xmin' xmax'];                 % Range of the parameters

xreal = prior(1,Nx,range);             % The reference parameters drawn from the prior distribution
yreal = forwardmodel(xreal);           % The true measurements
sd = yreal*0.1;                        % Standard deviation of measurement errors
Obs = yreal + sd.*randn(size(yreal));  % The measurements perturbed with white noise

p_k = 0.3;                             % The probability of using the Kalman proposal distribution
Z1  = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne);
fid_x = fopen('x.bin'); x_kzs = fread(fid_x,[Nx inf],'double'); fclose(fid_x); delete x.bin y.bin p.bin

p_k = 0;                               % Using the original dream_zs
Z2  = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne);
fid_x = fopen('x.bin'); x_zs = fread(fid_x,[Nx inf],'double'); fclose(fid_x); delete x.bin y.bin p.bin

save results
delete DREAM_KZS.mat                   % The intermediate results saved when running dream_kzs
