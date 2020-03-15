
N = 20;            % Number of parallel chains in MCMC
T = 4000;          % Number of iterations
Nx = 120;          % Dimension of model parameters
Ny = 243;          % Dimension of model responses
t1 = 100;          % After which iteration the Kalman proposal is used
t2 = ceil(0.3*T);  % After which iteration the Kalman proposal is not used
Ne = 200;          % Number of archive samples for the Kalman proposal

currentdir = pwd;
cd([currentdir,'\example']);
copyexample(N);    % Copy files for parallel computation
cd(currentdir);

range = repmat([-5 5],Nx,1);          % Range of the model parameters
xreal = importdata('xreal.dat');      % The reference parameters
yreal = forwardmodel(xreal);          % The true measurements
sd = 0.01*ones(Ny,1);                 % Standard deviation of measurement errors
Obs = yreal + sd.*randn(size(yreal)); % The measurements perturbed with white noise

p_k = 0.3;                            % The probability of using the Kalman proposal distribution
Z1  = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne);
fid_x = fopen('x.bin'); x_kzs = fread(fid_x,[Nx inf],'double'); fclose(fid_x); delete x.bin
fid_y = fopen('y.bin'); y_kzs = fread(fid_y,[Ny inf],'double'); fclose(fid_y); delete y.bin p.bin

p_k = 0;                               % Using the original dream_zs
Z2  = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne);
fid_x = fopen('x.bin'); x_zs = fread(fid_x,[Nx inf],'double'); fclose(fid_x); delete x.bin
fid_y = fopen('y.bin'); y_zs = fread(fid_y,[Ny inf],'double'); fclose(fid_y); delete y.bin p.bin

cd([currentdir,'\example']);
copyexample(N,-1);                     % Delete the files for parallel computation
cd(currentdir);

save results
delete DREAM_KZS.mat                   % The intermediate results saved when running dream_kzs
