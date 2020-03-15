function y = forwardmodel(x,ii)

if nargin < 2
    ii = 1;
end

[a,~] = size(x);
if a == 1
    x = x';
end

load data.mat

maindir = pwd;
exampledir = [maindir,'\example'];

% Delete the previous results
delete([exampledir,'\parallel_',num2str(ii),'\three_layers.hbs']);

% Generate the K field for the first layer
y1 = MeanY + Fi(:,1:kl_num)*x(1:kl_num,1);   % logK1
k1 = exp(y1);
k1 = reshape(k1,35,40);k1 = k1';

% Generate the K field for the second layer
y2 = MeanY + Fi(:,1:kl_num)*x(kl_num+1:kl_num*2,1); % logK2
k2 = exp(y2);
k2 = reshape(k2,35,40);k2 = k2';

% Generate the K field for the third layer
y3 = MeanY + Fi(:,1:kl_num)*x(kl_num*2+1:kl_num*3,1); % logK3
k3 = exp(y3);
k3 = reshape(k3,35,40);k3 = k3';

% Save the two K fields into the .mlt file
modifymlt(ii,exampledir,k1,k2,k3)

% Run the model
cd([exampledir,'\parallel_',num2str(ii)]);
system('modflow.bat');
cd(maindir);

% Import the results
fid = fopen([exampledir,'\parallel_',num2str(ii),'\three_layers.hbs']);
y = textscan(fid,'%f %f %s','headerlines',1);
y = y{1};
fclose(fid);

end

function modifymlt(ii,filedir,k1,k2,k3)

% Delete the previous .mlt file
delete ([filedir,'\parallel_',num2str(ii),'\three_layers.mlt']);

% three_layers_ref.mlt is the template
filein=[filedir,'\parallel_',num2str(ii),'\three_layers_ref.mlt'];
fileout=[filedir,'\parallel_',num2str(ii),'\three_layers.mlt'];
fidin=fopen(filein,'r');
fidout=fopen(fileout,'a+');

% Copy the first 4 lines from the template to three_layers.mlt
for i = 1:4
    writeline(filein,fileout,i);
end

% Save the conductivity field for the first layer into three_layers.mlt
for i = 1:40
    for j = 1:35
        if j == 35
            fprintf(fidout,'%14.7e\r\n',k1(i,j));
        else
            fprintf(fidout,'%14.7e ',k1(i,j));
        end
    end
end

% Copy lines 45-46 from the template to three_layers.mlt
for i = 45:46
    writeline(filein,fileout,i);
end

% Save the conductivity field for the second layer into three_layers.mlt
for i = 1:40
    for j = 1:35
        if j == 35
            fprintf(fidout,'%14.7e\r\n',k2(i,j));
        else
            fprintf(fidout,'%14.7e ',k2(i,j));
        end
    end
end

% Copy lines 87-88 from the template to three_layers.mlt
for i = 87:88
    writeline(filein,fileout,i);
end

% Save the conductivity field for the third layer into three_layers.mlt
for i = 1:40
    for j = 1:35
        if j == 35
            fprintf(fidout,'%14.7e\r\n',k3(i,j));
        else
            fprintf(fidout,'%14.7e ',k3(i,j));
        end
    end
end


fclose(fidin);
fclose(fidout);

end

