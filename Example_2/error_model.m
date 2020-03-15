function sd = error_model(x,Obs)


sd = x(end-1) + x(end)*Obs;

end

