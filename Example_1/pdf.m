function [log_post,Y] = pdf(X,Obs,sd,range)

% Change pdf.m and log_prior_pdf.m if needed 

N = size(X,1);
log_lik = nan(N,1);
Y = nan(N,length(Obs));

log_prior = log_prior_pdf(X,range);

for i = 1:N
    Y(i,:) = forwardmodel(X(i,:));
    log_lik(i,1) = log_lik_func(X(i,:),Y(i,:),Obs,sd);
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
xmin = range(:,1)';
xmax = range(:,2)';

for i = 1:N
    temp = nan(Nx,1);
    for j = 1:Nx
        temp(j,1) = log(1/(xmax(j)-xmin(j)));
    end
    log_prior(i,1) = sum(temp);
end

end

