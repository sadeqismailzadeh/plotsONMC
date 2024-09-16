%% load data
clc
clear all
path = ['E:\visualize\NewAnalysis\Lambda vs Distance\'];
% methods = ["Dynamic thinning", "Fukui-Todo"];

% i = 1;
% for method = methods
%     MethodName = sprintf('%s',method);
%     data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
%     error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);
%     i = i + 1;
% end


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_LambdaVsDistance'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;

%%---------------------------------------------------------------------------------
%% plot

data_path = [path, 'info.h5'];
Jmax = h5read(data_path, '/JmaxByDistance');
R0j = h5read(data_path, '/R0j');

beta = 1/0.849;
lambda = 4 * beta * Jmax;
lambdaCumul = cumsum(lambda);

timeFunction = zeros(size(lambdaCumul));

for i = 1:length(timeFunction)
    timeFunction(i) = i + (lambdaCumul(end)- lambdaCumul(i)) * 100 * 1.5;
end

[~ , idx] = min(timeFunction);





