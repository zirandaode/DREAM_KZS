function [nn,R_stat,MR_stat] = Convergence(x)

[T,din,N] = size(x);
tt = ceil(T/50):ceil(T/50):T; nn = tt*N;
R_stat = nan(length(tt),din);
MR_stat = nan(length(tt),1);
for t = 1:length(tt)
    t_end = tt(t);
    t_start = max(1,floor(0.5*t_end));
    [R_stat(t,:),MR_stat(t)] = Gelman(x(t_start:t_end,:,:),t_end,'DREAM_ZS');
end

end

function [R_stat,MR_stat] = Gelman(chain,t,ID)
% Calculates the \hat{R}-convergence diagnostic
% ----------------------------------------------------
% For more information please refer to:
% Gelman, A. and D.R. Rubin, (1992) Inference from Iterative Simulation
%      Using Multiple chain, Statistical Science, Volume 7, Issue 4,
%      457-472.
% Brooks, S.P. and A. Gelman, (1998) General Methods for Monitoring
%      Convergence of Iterative Simulations, Journal of Computational and
%      Graphical Statistics. Volume 7, 434-455. Note that this function
%      returns square-root definiton of R (see Gelman et al., (2003),
%      Bayesian Data Analsyis, p. 297).

% Written by Jasper A. Vrugt
% Los Alamos, August 2007
% ----------------------------------------------------

% Compute the dimensions of chain
[n,d,N] = size(chain); warning off; %#ok<*WNOFF>

if (n < 10)
    % Set the R-statistic to a large value
    [R_stat,MR_stat] = deal(NaN(1,d),NaN);
else
    
    % ------------ Univariate statistics ------------
    
    % STEP 0: Determine the chain means
    mean_chains = mean(chain); mean_chains = reshape(mean_chains(:),d,N)';
    
    % STEP 1: Determine the variance between the chain means
    B_uni = n * var(mean_chains);
    
    % STEP 2: Compute the variance of the various chains
    for ii = 1 : N, var_chains(ii,:) = var( chain(:,:,ii) ); end %#ok<*AGROW>
    
    % STEP 3: Calculate the average of the within chain variances
    W_uni = mean(var_chains);
    
    % STEP 4: Estimate the target variance
    sigma2 = ( (n-1)/n) * W_uni + (1/n) * B_uni;
    
    % STEP 5: Compute the R-statistic
    R_stat = sqrt( (N+1)/N * ( sigma2 ./ W_uni ) - (n-1)/(N*n) );
    
    % ---------- End univariate statistics ----------
    
    % ----------- Multivariate statistic ------------
    
    % STEP 1: Calculate the mean covariance W_mult of the m covariances of the chains
    W_mult = 0; for ii = 1 : N, W_mult = W_mult + cov( chain(1:n,1:d,ii) ) + eps * eye(d); end; W_mult = W_mult/N;
    
    % STEP 2: Calculate the covariance B of the m different means of the chains
    B_mult = cov(mean_chains) + eps * eye(d); % eps avoids problems with eig if var = 0
    
    % STEP 3: Calculate multivariate scale reduction factor, \hat{R}^{d}
    R = max ( abs ( eig ( W_mult \ B_mult ) ) );
    
    % STEP 4: Calculate the multivariate scale reduction factor, \hat{R}^d
    MR_stat = sqrt( (N+1)/N * R + (n-1)/n  );
    
    % --------- End multivariate statistic ----------
    
    % Now calculate the multivariate variant of this statistic
    % MR_stat = mpsrf_brooks(chain);
    
    % --------- Now check whether to write warning to warning_file ------
    
    [msgstr,msgid] = lastwarn;
    % if msgstr not empty
    if ~isempty(msgstr)
        % open file
        fid = fopen('warning_file.txt','a+');
        % now write to warning_file.txt
        evalstr = char(strcat(ID,{' '},'WARNING:',{' '},msgid,{' '},'in multivariate R_statistic of Brooks and Gelman at',{' '},num2str(t),{' '},'generations\n'));
        % Now print warning to file
        fprintf(fid,evalstr);
        % Close file
        fclose(fid);
        % Now remove warning
        lastwarn('');
    end
    
    % --------- End check whether to write warning to warning_file ------
    
end

end