%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
% DiffeRential Evolution Adaptive Metropolis with sampling from past archive
% and snooker update, originally written by Jasper A. Vrugt(jasper@uci.edu)
% Modified by Jiangjiang Zhang(zhangjiangjiang.dz@163.com) through introducing
% the Kalman update
%
% SYNOPSIS:                                                                                                
% Z = dream_kzs(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne)
% Z = dream_kzs_restart(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne)
%
% INPUT ARGUMENTS:
%         N: number of parallel chains in MCMC(N should be larger than or equal to 3)
%         T: number of iterations
%        Nx: dimension of model parameters
%        Ny: dimension of model responses
%       Obs: measurement data
%        sd: standard deviation of measurement error
%            if not available, sd = []
%     range: range of model parameters
%       p_k: probability of using the Kalman update, suggested value: 0.3
%            if p_k = 0, dream_kzs reduces to dream_zs
%        t1: after which iteration the Kalman update is used
%        t2: after which iteration the Kalman update is not used
%        Ne: number of archive samples for the Kalman update
%
% SAVED ARGUMENTS:
%         x: two-dimensional array with N chain trajectories (model parameters)
%            fid_x = fopen('x.bin'); x = fread(fid_x,[Nx inf],'double'); fclose(fid_x);
%            chain = permute(reshape(x',[Nx,N,T]),[3,1,2]); --> change 2-d to 3-d          
%         y: two-dimensional array with N chain trajectories (model outputs)
%            fid_y = fopen('y.bin'); y = fread(fid_y,[Ny inf],'double'); fclose(fid_y);
%         p: two-dimensional array with N chain trajectories (log-densities)
%            fid_p = fopen('p.bin'); p = fread(fid_p,[1 inf],'double'); fclose(fid_p);
%
% OUTPUT ARGUMENTS:  
%         Z: thinned chain history
%
% EXAMPLES:
% example 1: estimating model parameters of HYMOD
% example 2: estimating model parameters of hmodel
% example 3: estimating model parameters(Npar = 108) for a groundwater 
%            contaminant transport problem
% example 4: estimating model parameters(Npar = 120) for a three-dimensional 
%            groundwater flow problem
%
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
% Open parellel pool with parpool or matlabpool here

% Make sure that each time the random number is different
rand(floor(sum(100*clock)),1);clear;clc;

% Which example to run
Example = 1;

% Set paths
curr_dir = pwd; dream_dir = [curr_dir,'/DREAM_Dir']; 
exa_dir = [curr_dir,'/Example_',num2str(Example)];
addpath(dream_dir); cd(exa_dir)

% Run the example
run_this_1; 
% 'run_this_2.m' when sd is unknown(set as sd = []) and the error model 
% parameters(defined in 'error_model.m') are estimated together with the
% unknown model parameters

% Analyze the results
results_analysis

% Set paths
rmpath(dream_dir); cd(curr_dir)