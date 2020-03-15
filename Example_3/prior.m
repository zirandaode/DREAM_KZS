function X = prior(N,Nx,range)


X = nan(N,Nx);
if Nx == 108
    for i = 1:N
        for j = 1:8
            X(i,j) = unifrnd(range(j,1),range(j,2));
        end
        X(i,9:Nx) = randn(1,100);
    end
else
    for i = 1:N
        for j = 1:8
            X(i,j) = unifrnd(range(j,1),range(j,2));
        end
        X(i,9:Nx-2) = randn(1,100);
        X(i,Nx-1) = unifrnd(range(Nx-1,1),range(Nx-1,2));
        X(i,Nx) = unifrnd(range(Nx,1),range(Nx,2));
    end
end

end