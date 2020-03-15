function []= copyexample(NumPoints,flag)
if nargin < 2  
    parfor i=2:NumPoints
        copyfile('parallel_1',['parallel_',num2str(i)]);
    end    
elseif flag == -1
    parfor i = 2:NumPoints
        rmdir(['parallel_',num2str(i)],'s');
    end
end