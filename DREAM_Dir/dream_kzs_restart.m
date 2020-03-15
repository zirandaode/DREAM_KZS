function Z = dream_kzs_restart(N,T,Nx,Ny,Obs,sd,range,p_k,t1,t2,Ne)
% Please modify prior.m, pdf.m, error_model.m for your own problem

load DREAM_KZS
t_last = t; %#ok<*NODEF>

for t = t_last:T                                                       % Dynamic part: Evolution of N chains
    if mod(t,k) == 0, disp([num2str(t/T*100,'%4.2f'),'% finished']); end   % Show the process
    Xp = nan(N,Nx); p_acc = nan(N,1);                                      % Algorithmic variables
    dX = zeros(N,Nx);                                                      % Set N jump vectors to zero
    lambda = unifrnd(-c,c,N,1);                                            % Draw N lambda values
    std_X = std(X);                                                        % Compute std of each dimension
    R = randsample(1:m,N*n_d,'false'); R = reshape(R,n_d,N);               % Sample N*n_d integer values from [1,...,m]
    if t < t1 || t > t2                                                    % Mix of parallel and snooker updates
        method = randsample({'parallel' 'snooker'},1,'true',[1-p_s p_s]);
    else                                                                   % Mix of parallel, snooker and Kalman updates
        method = randsample({'parallel' 'snooker','Kalman'},1,'true',[1-p_s-p_k p_s p_k]);
    end
    for i = 1:N
        D = randsample(1:delta,1,'true');                                  % Select delta (equal selection probability)
        a = R(1:D,i); b = R(D+1:2*D,i); c_sn = R(2*D+1:3*D,i);             % Define a and b (parallel) + c (snooker)
        if strcmp(method,'parallel')                                       % Parallel direction update
            id(i) = randsample(1:n_CR,1,'true',p_CR);                      % Select index of crossover value
            z = rand(1,Nx);                                                % Draw Nx values from U[0,1]
            A = find(z < CR(id(i)));                                       % Derive subset A selected dimensions
            d_star = numel(A);                                             % How many dimensions sampled?
            if d_star == 0, [~,A] = min(z); d_star = 1; end                % A must contain at least one value
            gamma_d = 2.38/sqrt(2*D*d_star);                               % Calculate jump rate
            g = randsample([gamma_d 1],1,'true',[1-p_g p_g]);              % Select gamma: 80/20 mix [default 1]
            dX(i,A) = c_star*randn(1,d_star) + ...
                (1+lambda(i))*g*sum(Z(a,A)-Z(b,A),1);                      % Compute the ith jump differential evolution
        elseif strcmp(method,'snooker')                                    % Snooker update
            id(i) = n_CR;                                                  % Full crossover
            F = X(i,1:Nx)-Z(a,1:Nx); D_e = max(F*F',1e-300);               % Define projection X(i,1:d) - Z(a,1:d)
            zP = F*( sum((Z(b,1:Nx)-Z(c_sn,1:Nx)).*F )/D_e );              % Orthogonally project zR1 and zR2 onto F
            g = 1.2 + rand;                                                % Determine jump rate
            dX(i,1:Nx) = c_star*randn(1,Nx) + (1+lambda(i))*g*zP;          % Calculate the ith jump snooker update
        elseif strcmp(method,'Kalman')                                     % Kalman update
            f = fopen('x.bin'); x = fread(f,[Nx inf],'double'); fclose(f); % Read saved chains (parameters)
            M = size(x,2); if M < Ne, Ne = M; end                          % If Ne is larger than the number of saved samples
            a1 = randsample(ceil(M/2)-1:M,ceil(Ne/2));
            a2 = randsample(setdiff(1:M,a1),Ne-ceil(Ne/2));
            a3 = [a1,a2];                                                  % Generate archive samples for the Kalman proposal
            xf = x(:,a3); clear x;                                         % Ensemble of the model parameters
            f = fopen('y.bin'); y = fread(f,[Ny inf],'double'); fclose(f); % Read saved chains (responses)
            yf = y(:,a3); clear y;                                         % Ensemble of the model responses
            if isempty(sd)                                                 % Standard deviation of measurement errors is not available
                sd_kf = error_model(X(i,:),Obs);                           % Estimate the standard deviation, you may modify error_model.m
            else
                sd_kf = sd;                                                % Standard deviation of measurement errors is available
            end
            Cd = diag(sd_kf.^2);                                           % Covariance of the measurement errors
            meanxf = repmat(mean(xf,2),1,Ne);                              % Mean of xf
            meanyf = repmat(mean(yf,2),1,Ne);                              % Mean of yf
            Cxy =  (xf-meanxf)*(yf-meanyf)'/(Ne-1);                        % Cross-covariance between xf and yf
            Cyy =  (yf-meanyf)*(yf-meanyf)'/(Ne-1);                        % Auto-covariance of yf
            K = Cxy/(Cyy+Cd);                                              % Kalman gain
            if t < T*2/3, l = 1; else, l = sign(2*rand-1); end             % After 2/3 of total iterations, use +/- updates to maintain detailed balance
            dX(i,1:Nx) = l*K*(Obs+randn(Ny,1).*sd_kf-Y(i,:)');             % Calculate the Kalman update
        end
        Xp(i,:) = X(i,1:Nx) + dX(i,1:Nx);                                  % Compute the ith proposal
        Xp(i,:) = Boundary_handling(Xp(i,:),range,'reflect');              % Boundary handling
    end
    [p_Xp,Yp] = pdf(Xp,Obs,sd,range);                                      % Calculate the log-density and model responses
    if strcmp(method,'snooker')                                            % Snooker correction: non-symmetry proposal distribution
        alfa_sn = (sum((Xp - Z(R(1,1:N),1:Nx)).^2,2)./sum((X - ...
            Z(R(1,1:N),1:Nx)).^2,2)).^((Nx-1)/2);
    else                                                                   % Otherwise no correction needed
        alfa_sn = ones(N,1);
    end
    for i = 1:N                                                            % Accept/reject proposals
        p_acc(i) = min(1,alfa_sn(i)*exp(p_Xp(i,1)-P(i,1)));                % Compute acceptance probability
        if p_acc(i) > rand                                                 % p_acc(i) larger than U[0,1]?
            X(i,:) = Xp(i,:); Y(i,:) = Yp(i,:); P(i,1) = p_Xp(i,1);        % True: Accept proposal
        else
            dX(i,1:Nx) = 0;                                                % Set jump back to zero for p_CR
        end
        J(id(i)) = J(id(i)) + sum((dX(i,1:Nx)./std_X).^2);                 % Update jump distance crossover idx
        n_id(id(i)) = n_id(id(i)) + 1;                                     % How many times idx crossover used
    end
    if (mod(t,k) == 0)                                                     % Check whether to append X to archive Z
        Z(m+1:m+N,1:Nx+1) = [X P]; m = m + N;                              % Append current values to Z after k generations
        if t < T/10
            p_CR = J./n_id; p_CR = p_CR/sum(p_CR);                         % Update selection probability crossover
        end
    end
    f = fopen('x.bin','a+','n'); fwrite(f,X','double'); fclose(f);         % Store model parameters of new population
    f = fopen('y.bin','a+','n'); fwrite(f,Y','double'); fclose(f);         % Store model responses of new population
    f = fopen('p.bin','a+','n'); fwrite(f,P','double'); fclose(f);         % Store log-densities of new population
    if mod(t,ceil(T/50)) == 0
        save DREAM_KZS.mat                                                 % Save intermediate results
    end
end

end

function x = Boundary_handling(x,range,flag)
% Function to check whether parameter values remain within prior bounds

% First determine the size of x
[m,~] = size(x);
% Now replicate min and max
min_d = repmat(range(:,1)',m,1); max_d = repmat(range(:,2)',m,1);
% Now find which elements of x are smaller than their respective bound
[ii_low] = find(x < min_d);
% Now find which elements of x are larger than their respective bound
[ii_up] = find(x > max_d);

switch flag
    case 'reflect'  % reflection
        % reflect in min
        x(ii_low)= 2 * min_d(ii_low) - x(ii_low);
        % reflect in max
        x(ii_up)= 2 * max_d(ii_up) - x(ii_up);
    case 'bound'    % set to bound
        % set lower values to min
        x(ii_low)= min_d(ii_low);
        % set upper values to max
        x(ii_up)= max_d(ii_up);
    case 'fold'     % folNparg
        % Fold parameter space lower values
        x(ii_low) = max_d(ii_low) - ( min_d(ii_low) - x(ii_low) );
        % Fold parameter space upper values
        x(ii_up) = min_d(ii_up) + ( x(ii_up) - max_d(ii_up) );
        % Now double check in case elements are still out of bound -- this is
        % theoretically possible if values are very small or large
    otherwise
        disp('do not know this boundary handling option! - treat as unbounded parameter space')
end

% Now double check if all elements are within bounds
% If parameter so far out of space then possible that reflection or
% folding still creates point that is out of bound - this is very
% rare but needs explicit consideration
if strcmp(flag,'reflect') || strcmp(flag,'fold')
    % Lower bound
    [ii_low] = find(x < min_d); x(ii_low) = min_d(ii_low) + rand(size(ii_low)).* ( max_d(ii_low) - min_d(ii_low) );
    % Upper bound
    [ii_up]  = find(x > max_d); x(ii_up)  = min_d(ii_up)  + rand(size(ii_up)).*  ( max_d(ii_up)  - min_d(ii_up)  );
end

end