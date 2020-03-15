function X = prior(N,Nx,range)


X = nan(N,Nx);
if Nx == 120
    for i = 1:N
        X(i,:) = randn(1,Nx);
    end
else
    for i = 1:N
        X(i,1:Nx-1) = randn(1,Nx-1);
        X(i,Nx) = unifrnd(range(Nx,1),range(Nx,2));
    end
end

end