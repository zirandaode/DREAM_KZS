function [log_post,Y] = pdf(X,Obs,sd,range)

% Change pdf.m and log_prior_pdf.m if needed

N = size(X,1);
log_lik = nan(N,1);
Y = nan(N,length(Obs));

log_prior = log_prior_pdf(X,range);

parfor i = 1:N
    Y(i,:) = forwardmodel(X(i,:),i);
    log_lik(i,1) = log_lik_func(X(i,:),Y(i,:),Obs,sd)
end

log_post = log_prior + log_lik;

end

function log_lik = log_lik_func(x,y,Obs,sd)

if sum( size(Obs) == size(y) ) == 0
    y = y';
end

Err = Obs - y;

if isempty(sd)
    sd = error_model(x,Obs);
end

log_lik = - ( length(Err) / 2) * log(2 * pi) - sum ( log( sd ) ) - ...
    1/2 * sum ( ( Err./sd ).^2);

end

function log_prior = log_prior_pdf(X,range)

% Calculate the log pdf of X in the prior distribution

N = size(X,1);
Nx = size(X,2);
log_prior = nan(N,1);

for i = 1:N
    temp = nan(Nx,1);
    if Nx == 108
        for j = 1:8
            temp(j,1) = log(1/(range(j,2)-range(j,1)));
        end
        for j = 9:Nx
            temp(j,1) = log(normpdf(X(i,j),0,1));
        end
    else
        for j = 1:8
            temp(j,1) = log(1/(range(j,2)-range(j,1)));
        end
        for j = 9:Nx-2
            temp(j,1) = log(normpdf(X(i,j),0,1));
        end
        temp(Nx-1,1) = log(1/(range(Nx-1,2)-range(Nx-1,1)));
        temp(Nx,1) = log(1/(range(Nx,2)-range(Nx,1)));
    end
    log_prior(i,1) = sum(temp);
end

end


