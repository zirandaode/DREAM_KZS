function []= copyexample(NumPoints,flag)
str = 'parallel_';
if nargin < 2  
    parfor i=2:NumPoints
        copyfile([str,'1'],[str,num2str(i)]);
    end    
elseif flag == -1
    parfor i = 2:NumPoints
        rmdir([str,num2str(i)],'s');
    end
end