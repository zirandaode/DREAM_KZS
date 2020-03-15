
N = 20;                % Number of parallel chains in MCMC
T = 5000;              % Number of iterations
Nx = 110;              % Dimension of model parameters
Ny = 150;              % Dimension of model responses
t1 = 100;              % After which iteration the Kalman proposal is used
t2 = ceil(0.3*T);      % After which iteration the Kalman proposal is not used
Ne = 300;              % Number of archive samples for the Kalman proposal

currentdir = pwd;
cd([currentdir,'\example']);
copyexample(N);        % Copy files for parallel computation
cd(currentdir);

load obsdata
sd = [];               % The standard deviation of measurement errors is unknown
range = [range;[0 0.01];[0 0.01]]; % Range of the model and error model parameters
Obs = yreal + err;     % The measurements perturbed with white noise

p_k = 0.3;             % The probability of using the Kalman proposal distribution
Z1  = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne);
fid_x = fopen('x.bin'); x_kzs = fread(fid_x,[Nx inf],'double'); fclose(fid_x); 
movefile('x.bin','x_kzs.bin'); movefile('y.bin','y_kzs.bin'); 
movefile('p.bin','p_kzs.bin'); movefile('DREAM_KZS.mat','TEMP_KZS.mat');
save results_dream_kzs.mat

p_k = 0;  T = 50000;   % Using the original dream_zs with longer chains
Z2  = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne);
fid_x = fopen('x.bin'); x_zs = fread(fid_x,[Nx inf],'double'); fclose(fid_x); 
movefile('x.bin','x_zs.bin'); movefile('y.bin','y_zs.bin'); 
movefile('p.bin','p_zs.bin'); movefile('DREAM_KZS.mat','TEMP_ZS.mat');

cd([currentdir,'\example']);
copyexample(N,-1);     % Delete the files for parallel computation
cd(currentdir);

save results
delete DREAM_KZS.mat   % The intermediate results saved when running dream_kzs
