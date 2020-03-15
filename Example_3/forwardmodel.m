function  [yout]=forwardmodel( x1,ii )

if nargin<2
    ii = 1;
end

currentdir = pwd;
exampledir = [currentdir,'\example'];

XY = importdata('XY.dat');
timesteps = [ 4 5 6  7 8 9 10 11 12];

nrow = 41;
ncol = 81;
nc=nrow*ncol;
nobs = max(size(XY));

xx = x1(1);
yy = x1(2);
Ss = x1(3:8);
klterm = x1(9:108);

[dx1,dx2,nx] = getint(20,0.25,xx );    
icol=zeros(81,1);
icol(1:nx)=dx1;
icol(nx+1:end)=dx2;

[dy1,dy2,ny] = getint(10,0.25,yy );    
irow=zeros(41,1);
irow(1:ny)=dy1;
irow(ny+1:end)=dy2;

modifydis(ii,exampledir,icol,irow);
modifybtn(ii,exampledir,icol,irow,timesteps);

tt = zeros(18,1);
tt(1) = ny; % for irow
tt(2) = nx; % for icol
tt(3) = Ss(1);
tt(4) = Ss(2);
tt(5) = Ss(3);
tt(6) = Ss(4);
tt(7) = Ss(5);
tt(8) = Ss(6);
tt(9:18) = 0; %#ok<NASGU>
delete([exampledir,'\parallel_',num2str(ii),'\mt.dat']);
save([exampledir,'\parallel_',num2str(ii),'\mt.dat'],'tt','-ascii');

load fun.mat
ty = ones(41*81,1)*2;
for j = 1:41*81;
    for k = 1:length(klterm);
        ty(j) = ty(j) + klterm(k) * fun(j,k);
    end;
end; 
tyy = reshape(ty,81,41);
tyy=tyy';
yy=exp(tyy);
yyy = reshape(yy,3321,1);
delete([exampledir,'\parallel_',num2str(ii),'\Cond_01.dat']);
delete([exampledir,'\parallel_',num2str(ii),'\Cond.dat']);
dlmwrite([exampledir,'\parallel_',num2str(ii),'\Cond_01.dat'],yy,'delimiter', '', 'precision', '%10.4f','newline', 'pc');
dlmwrite([exampledir,'\parallel_',num2str(ii),'\Cond.dat'],yyy,'delimiter', '', 'precision', '%10.4f','newline', 'pc');

cd([exampledir,'\parallel_',num2str(ii)]);
delete([exampledir,'\parallel_',num2str(ii),'\conc.dat']);
delete([exampledir,'\parallel_',num2str(ii),'\head.dat']);
system('MT3D_mat');
cd(currentdir);

cout = importdata([exampledir,'\parallel_',num2str(ii),'\conc.dat']);
hout = importdata([exampledir,'\parallel_',num2str(ii),'\head.dat']);

for j = 1:length(timesteps);            
    for k = 1:nc;             
        conc(k,j) = cout((j-1)*nc+k,3);        
    end;    
end; 

for j = 1:length(timesteps);  
    tpv = zeros(1,nc);    
    tpv(1:nc) = conc(:,j);    
    tpcon = mfdata(tpv,nrow,ncol);    
    for k = 1:nobs    
        y(nobs*j-nobs+k,1) = getobsval(dx1,dx2,nx,ncol,dy1,dy2,ny,nrow,tpcon,XY(1,k),XY(2,k));        
    end;    
end;

tpv = zeros(1,nc);    
tpv(1:nc) = hout(:,3);    
tphead = mfdata(tpv,nrow,ncol); 
for k = 1:nobs    
    h(k,1) = getobsval(dx1,dx2,nx,ncol,dy1,dy2,ny,nrow,tphead,XY(1,k),XY(2,k));        
end;

yout = [y;h];

end


function modifydis(ii,filedir,icol,irow)

delete ([filedir,'\parallel_',num2str(ii),'\MODFO.DIS']);
filein=[filedir,'\parallel_',num2str(ii),'\MODFO_as.DIS'];
fileout=[filedir,'\parallel_',num2str(ii),'\MODFO.DIS'];
fidin=fopen(filein,'r');
fidout=fopen(fileout,'a+');

for i=1:3
    writeline(filein,fileout,i);
end

for i=1:8
    for j=1:9
        fprintf(fidout,'%10.6f',icol(i*10+j-10));
    end
    fprintf(fidout,'%10.6f\r\n',icol(i*10));
end

fprintf(fidout,'%10.6f\r\n',icol(81));

for i=13
    writeline(filein,fileout,i);
end

for i=1:4
    for j=1:9
        fprintf(fidout,'%10.6f',irow(i*10+j-10));
    end
    fprintf(fidout,'%10.6f\r\n',irow(i*10));
end

fprintf(fidout,'%10.6f\r\n',irow(41));

for i=19:21
    writeline(filein,fileout,i);
end

fclose(fidin);
fclose(fidout);

end


function modifybtn(ii,filedir,icol,irow,timesteps)

delete ([filedir,'\parallel_',num2str(ii),'\modfo.btn']);
filein= [filedir,'\parallel_',num2str(ii),'\modfo_as.btn'];
fileout= [filedir,'\parallel_',num2str(ii),'\modfo.btn'];
fidin=fopen(filein,'r');
fidout=fopen(fileout,'a+');

for i=1:7
    writeline(filein,fileout,i);
end

for i=1:80
    fprintf(fidout,'%10.6f',icol(i));
end

fprintf(fidout,'%10.6f\r\n',icol(81));

for i=9
    writeline(filein,fileout,i);
end

for i=1:40
    fprintf(fidout,'%10.6f',irow(i));
end

fprintf(fidout,'%10.6f\r\n',irow(41));

for i=11:17
    writeline(filein,fileout,i);
end

fprintf(fidout,'%10.0f\r\n',length(timesteps));

for i=1:length(timesteps)-1
    fprintf(fidout,'%10.0f',timesteps(i));
end

fprintf(fidout,'%10.0f\r\n',timesteps(length(timesteps)));

for i=20:55
    writeline(filein,fileout,i);
end

fclose(fidin);
fclose(fidout);

end


function dataout = writeline(filein,fileout,line)

fidin=fopen(filein,'r');
fidout=fopen(fileout,'a+');
nline=0;
while ~feof(fidin)
    tline=fgetl(fidin); 
    nline=nline+1;
    if nline==line
        fprintf(fidout,'%s\r\n',tline);
        dataout=tline;
    end
end
fclose(fidin);
fclose(fidout);

end


function [dx1,dx2,n1] = getint(L,dx,x)

N = round(L/dx);
n1 = floor(x/(dx))+1;
dx1 = x./(n1-1);
n2 = N + 1 - n1;
dx2 = (L - (n1-0.5).*dx1)./(n2-0.5);

end


function obsv = getobsval(dx1,dx2,nx1,ncx,dy1,dy2,ny1,ncy,conc,xobs,yobs)

[m,n]=size(conc);
if (m~=ncy || n ~= ncx)
    printf('error, check dim');
    pause;
end
xx = [ dx1*ones(1,nx1-1) 0.5*(dx1+dx2) dx2*ones(1,ncx-nx1-1)];
xc = [0 cumsum(xx)]; 
yy = [ dy1*ones(1,ny1-1) 0.5*(dy1+dy2) dy2*ones(1,ncy-ny1-1)];
yc = [0 cumsum(yy)]; 
obsv = interp2(xc,yc,conc,xobs,yobs,'spline');

end


function y = mfdata(x,m,n)
% x is a vector with length m*n, this function transform the vector to a
% (m,n) matrix
y = zeros(m,n);

for k = 1:m;
    for l = 1:n;
        y(k,l)=x((k-1)*n+l);
    end;
end;

end



