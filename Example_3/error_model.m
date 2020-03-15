function sd = error_model(x,Obs)


sd = [x(end-1)*ones(135,1);x(end)*ones(15,1)];

end

