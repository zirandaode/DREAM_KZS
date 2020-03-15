function [] = pre_settings()

%% Generate files for the KL expansion
MeanY = -6.5;
kl_num = 40;
[V,Fi] = KLexpansion();
save data.mat Fi V MeanY kl_num

%% Generate N_obs observation locations

Lx = 35;Ly = 40;
exampledir = [pwd,'\example'];
ii = 1;
delete([exampledir,'\parallel_',num2str(ii),'\three_layers.hob']);

filein=[exampledir,'\parallel_',num2str(ii),'\three_layers_ref.hob'];
fileout=[exampledir,'\parallel_',num2str(ii),'\three_layers.hob'];

fidin = fopen(filein,'r');
fidout = fopen(fileout,'a+');

N_obs = 81;
% Nt = (Lx-2)*(Ly-2);
% [y,x] = meshgrid(2:Ly-1,2:Lx-1);
% x = x(:);
% y = y(:);
% rs = randsample(1:Nt,N_obs);
% x = x(rs,1);
% y = y(rs,1);
[y,x] = meshgrid(3:4:38,5:3:31);
x = x(:);
y = y(:);
save location.mat x y

% Copy the first line from the template to three_layers.mlt
writeline(filein,fileout,1);

fprintf(fidout,'%1d %s\r\n',N_obs*3,'0 2 14 999.0  # Data Set 1: NH MOBS MAXM IUHOBSV HOBDRY');

% Copy the third line from the template to three_layers.mlt
writeline(filein,fileout,3);

% Write the N_obs locations into three_layers.hob
for i = 1:N_obs
    fprintf(fidout,'%s %1d %2d %2d %s\r\n',['obs1_',num2str(i,'%.2d')],1,y(i),x(i),...
        '1 1.0  0.0000   0.0000  0.0000');
end

for i = 1:N_obs
    fprintf(fidout,'%s %1d %2d %2d %s\r\n',['obs2_',num2str(i,'%.2d')],2,y(i),x(i),...
        '1 1.0  0.0000   0.0000  0.0000');
end

for i = 1:N_obs
    fprintf(fidout,'%s %1d %2d %2d %s\r\n',['obs3_',num2str(i,'%.2d')],3,y(i),x(i),...
        '1 1.0  0.0000   0.0000  0.0000');
end

fclose(fidin);
fclose(fidout);

end

function [V,Fi] = KLexpansion()

% Field parameters
Lx = 35*1.5;Ly = 40*2;Lz = Ly;         % the size of domian
Ncol = 35;Nrow = 40;Nlay = Nrow; % the discretization of domain
CorLengthx = 25*1.5;CorLengthy = 30*2;CorLengthz = CorLengthy;% correlation length in each dimension
VarF = 0.5;                      % variance of field
dx = Lx/(Ncol-1)*ones(1,Ncol);
dy = Lz/(Nrow-1)*ones(1,Nrow);
dz = Lz/(Nrow-1)*ones(1,Nrow);   % the size of cells

% Settings for KL expansion
Npd = 1000;                      % the total number of KL terms
Npr = 1;
Nroot = 40;                      %超越方程的个数
Np = 40;                         % how many terms used

% construct the KL expansion files
zz = 0;
for k = 1:Nlay
    yy = 0;
    for j = 1:Nrow
        xx = 0;
        for i=1:Ncol
            N=i+(j-1)*Ncol+(k-1)*Ncol*Nrow;
            xa(N) = xx; %#ok<*SAGROW>
            ya(N) = yy;
            za(N) = zz;
            xx = xx+dx(i);
        end
        yy = yy+dy(j);
    end
    zz = zz+dz(k);
end

% eigenpairAna_2D
x0RLn = 0.1;
dXRLn = 1.0e-5;
z0RLn = 0.1;
dZRLn = 1.0e-5;
Rx = CorLengthx/Lx;
Rz = CorLengthz/Lz;
[wx,~] = Findwxz(Nroot,x0RLn,dXRLn,1.0,Rx);
[wz,~] = Findwxz(Nroot,z0RLn,dZRLn,1.0,Rz);

for i=1:Nroot
    wx(i) = wx(i)/Lx;
    wz(i) = wz(i)/Lz;
end

RLn = zeros(4,Nroot^2);
for i = 1:Nroot
    for j = 1:Nroot
        k = (i-1)*Nroot+j;
        RLn(1,k) = wx(i);
        RLn(2,k) = wz(j);
        RLn(3,k) = 4.0*CorLengthx*CorLengthz*VarF/...
            ((CorLengthx*wx(i))^2+1.)/((CorLengthz*wz(j))^2+1.);
        if(k > 1)
            for kk = k:-1:2
                if(RLn(3,kk)>RLn(3,kk-1))
                    r1 = RLn(1,kk);
                    r2 = RLn(2,kk);
                    r3 = RLn(3,kk);
                    RLn(1,kk) = RLn(1,kk-1);
                    RLn(2,kk) = RLn(2,kk-1);
                    RLn(3,kk) = RLn(3,kk-1);
                    RLn(1,kk-1) = r1;
                    RLn(2,kk-1) = r2;
                    RLn(3,kk-1) = r3;
                end
            end
        end
    end
end

Fi = zeros(Npd,Ncol*Nlay);
for k = 1:Np
    z = 0;
    for j = 1:Nrow
        x = 0;
        for i = 1:Ncol
            N = (j-1)*Ncol+i;
            r1 = CorLengthx*RLn(1,k);
            r2 = CorLengthz*RLn(2,k);
            r4 = RLn(1,k)*x;
            r5 = RLn(2,k)*z;
            r3 = (r1*cos(r4)+sin(r4))*(r2*cos(r5)+sin(r5))/...
                sqrt(((r1*r1+1.)*Lx/2.+CorLengthx)*...
                ((r2*r2+1.)*Lz/2.+CorLengthz));
            Fi(k,N) = sqrt(RLn(3,k))*r3;
            x = x+dx(i);
        end
        z = z+dz(j);
    end
end

Eigenv = zeros(1,Npd);
for i = 1:Npd
    Eigenv(i) = RLn(3,i);
end
V = sum(Eigenv(1:Np))/(VarF*Lx*Lz);

Fi = Fi';

end

function f = fout( x,C,L )

f = (C^2*x^2-1)*sin(x*L)-2*C*x*cos(x*L);

end

function [x,y] = Findwxz(Nroot,x0,dx,L,CorLength)

accuricywn = 1e-5;
Nr = 1;
x1 = x0;
f0 = fout(x0,CorLength,L);

while(Nr <= Nroot)
    x1 = x1+dx;
    f1 = fout(x1,CorLength,L);
    if (f1 == 0)
        x(Nr) = x1; %#ok<*AGROW>
        y(Nr) = fout(x(Nr),CorLength,L);
        Nr = Nr+1;
        x0 = x1+0.5*dx;
        f0 = fout(x0,CorLength,L);
    elseif(f0*f1 < 0)
        x3 = x1;
        while(abs(x0-x1) > accuricywn)
            x2 = (x0+x1)*0.5;
            f2 = fout(x2,CorLength,L);
            if(f2 == 0)
                x(Nr) = x2;
                y(Nr) = f2;
                Nr = Nr+1;
                x0 = x1+0.5*dx;
                f0 = fout(x0,CorLength,L);
                % go to 10?
            elseif(f1*f2 < 0)
                x0 = x2;
                f0 = f2; %#ok<*NASGU>
            else
                x1 = x2;
                f1 = f2;
            end
        end
        x(Nr) = x1-0.5*abs(x1-x0);
        y(Nr) = fout(x(Nr),CorLength,L);
        Nr = Nr+1;
        % 10 continue?
        x0 = x3;
        f0 = fout(x3,CorLength,L);
    end
end

end




