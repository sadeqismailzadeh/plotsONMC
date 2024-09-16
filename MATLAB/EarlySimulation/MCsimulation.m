%  Monte carlo simulation
clc
clear all;
addpath(pwd);
addpath([pwd,'\mc']);
addpath([pwd,'\lattice2d']);
% run('largecell2D.m')
tic();
largecell2D


toc();
%%
tic();
% T_E_C_op = [];
number_of_runs = 20;
getdata_T_steps = 1; % to not get data  in every temperature

for run = 1: number_of_runs
%     largecell2Ddisorder
    
    %  dipoles is 3 * n stores dipole moments' vectors.
    %  i.e dipoles(5,1) is x component of 5th dipole.
    dipoles = rand(3*n, 1) - 0.5;
    
    %  normalize magnitude of each moment to 1
    for i = 1:n
        dipoles(3*i-2: 3*i) = dipoles(3*i-2: 3*i)  / norm(dipoles(3*i-2: 3*i));
    end
%     dipoles
    
    datarun = [];
%     for T = 3: -0.05 :0.05
%         dipoles = stabilize(dipoles, KpcSym, T, 30);
%         [dipoles, meanE, C, op, X] = getStatisticalData(dipoles, KpcSym, KpcUT, T, 500);
%         datarun = [datarun; T, meanE, C, op, X];
%         [T meanE, C, op, X];
%         T
%     end
    
    T = 2;
    number_T_runs = 0
%     for T=3.5:-0.1:0
    while T > 0.02
        number_T_runs = number_T_runs  + 1;
        dipoles = stabilize(dipoles, KpcSym, T, 5000);
        
%         if mod(number_T_runs, getdata_T_steps) == 0
%             [dipoles, meanE, C, op, X] = getStatisticalData(dipoles, KpcSym, KpcUT, T, 10000, 2);
%             datarun = [datarun; T, meanE, C, op, X];
%             disp('here')
%         end
%         [T meanE, C, op, X]
        [run T]
%         if T > 0.2
%             T = T * 0.93;
%         else
%             T = T - 0.025;
%         end
        T= T-0.02;
            
    end
    
    if run == 1
        T_E_C_op = datarun;
    else
        T_E_C_op = T_E_C_op + datarun;
    end
    
end

T_E_C_op = T_E_C_op / number_of_runs;
t=toc()


%% save date
date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
folder_address = ['../results/' date_time_string];
[status, msg, msgID] = mkdir(folder_address)
save([folder_address, '/all_variables.mat'])

%% visualize
MCvisualize

%% time
disp(['elapsed time elapseed is ', datestr(datenum(0,0,0,0,0,t),'HH:MM:SS')]) 


