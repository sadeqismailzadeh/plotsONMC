%% load data
% TODO for error estimate sigma / sqrt(N - 1)
clear all
clc
addpath(pwd);
addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);

path = ['E:/visualize/inputHistogram/histogram/'];
date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_Histogram'];
[status, msg, msgID] = mkdir(save_address)

Llist = readmatrix([path, 'list.txt'])

Nlist = Llist.^2;
N = length(Nlist);

% data = readmatrix([path 'MeanX.csv']);
% error = readmatrix([path 'errorX.csv']);
myxLabel = '$T$';
for i=1:N
    EnsembleIndex = num2str(Llist(i));
    data_path = [path, EnsembleIndex , '/'];

    % [status, msg, msgID] = mkdir(folder_address)

%     dataRegular(:,:,i) = dlmread([data_path 'MeanX.txt']);
%     errorRegular(:,:,i) = readmatrix([data_path 'errorX.txt']);
end

dpi = 1200;

% data = dlmread([path 'MeanX.txt']);
% error = dlmread([path 'errorX.txt']);

% orderParam = data(:,1);
% myxLabel = '$B\mu /\Lambda $';


% error = error / sqrt(10); % for results before 2022.08.25


% locations = readmatrix([path, '0/locations.txt']);
% DistanceVec = readmatrix([path, 'DistanceVec.txt']);
% % Plabel = readmatrix([path, 'Plabel.txt']);
% Mean_GVec = readmatrix([path, 'Mean_GVec.txt']);
% Mean_GConVec = readmatrix([path, 'Mean_GConVec.txt']);
% error_GVec = readmatrix([path, 'error_GVec.txt']);
% error_GConVec = readmatrix([path, 'error_GConVec.txt']);
%%--------------------------------------------------------------------------------------
%% save address (for loaded data)
date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['d:/visualize/output/', date_time_string , '_Histogram'];
[status, msg, msgID] = mkdir(save_address)

%%---------------------------------------------------------------------------------------
%% read size
Ebase = 200000;  % remember to change this!!!!!!!!!
dT = 0.001;
% T = 0.845:dT:0.865;
T = [0.840:0.0002:0.8458, 0.846:0.0001:0.855, 0.8552:0.0002:0.86, 0.8605:0.0005:1];
disp(['size(T) = ',num2str(size(T))]);
Nensemble = 12;  % remember to change this!!!!!!!!!
T0L = zeros(1,length(Llist));

E1TEnsL = zeros(length(T),Nensemble,length(Llist));
E2TEnsL = zeros(length(T),Nensemble,length(Llist));
E3TEnsL = zeros(length(T),Nensemble,length(Llist));
E4TEnsL = zeros(length(T),Nensemble,length(Llist));
m1TEnsL = zeros(length(T),Nensemble,length(Llist));
m2TEnsL = zeros(length(T),Nensemble,length(Llist));
m3TEnsL = zeros(length(T),Nensemble,length(Llist));
m4TEnsL = zeros(length(T),Nensemble,length(Llist));
m1ETEnsL = zeros(length(T),Nensemble,length(Llist));
m2ETEnsL = zeros(length(T),Nensemble,length(Llist));
m3ETEnsL = zeros(length(T),Nensemble,length(Llist));
m4ETEnsL = zeros(length(T),Nensemble,length(Llist));

tic;
for k = 1:length(Llist)
    sizeIndex = num2str(Llist(k))
    data = readmatrix([path, sizeIndex,  '/', 'MeanX.csv']);
    Etime = [];
    mtime = [];
    for i=1:Nensemble
        EnsembleIndex = num2str(i-1);
        data_path = [path, sizeIndex, '/'];
        EtimeTemporary = h5read([data_path, EnsembleIndex, '/EnsembleResults.h5'],['/TimeSeries/Energy']);
        mtimeTemporary = h5read([data_path, EnsembleIndex, '/EnsembleResults.h5'],['/TimeSeries/Magnetization']);
        % 50 must be changed 5 or 4 because tau is calculated based on 1 monte carlo step
        % but data is stored every 11 monte carlo steps
        % for testing I could use 50 because it accord to tau of big sizes which are 10000

        % mtimeTemporary = mtimeTemporary(1:5:end);
        % EtimeTemporary = EtimeTemporary(1:5:end);
        if Llist(k) < 32  % TODO this
            mtimeTemporary = mtimeTemporary(1:5:end);
            EtimeTemporary = EtimeTemporary(1:5:end);
        elseif Llist(k) < 64
            mtimeTemporary = mtimeTemporary(1:15:end);
            EtimeTemporary = EtimeTemporary(1:15:end);
        else
            mtimeTemporary = mtimeTemporary(1:50:end);
            EtimeTemporary = EtimeTemporary(1:50:end);
        end

        Etime(:,i) = EtimeTemporary;
        mtime(:,i) = mtimeTemporary;

    end
    Etime = double(Etime);
    mtime = double(mtime);
    Etime = Etime + Ebase;
    T0 = data(1)
    % T(end+1) = T0   or somthing like this.
    T0L(k) = T0;
    beta0 = 1/T0;

    onesVector = ones(size(Etime(:,1)));

    E1L = zeros(size(T));
    E2L = zeros(size(T));
    E3L = zeros(size(T));
    E4L = zeros(size(T));
    m1L = zeros(size(T));
    m2L = zeros(size(T));
    m3L = zeros(size(T));
    m4L = zeros(size(T));
    m1EL = zeros(size(T));
    m2EL = zeros(size(T));
    m3EL = zeros(size(T));
    m4EL = zeros(size(T));


    for i=1:Nensemble
        disp(['L = ',num2str(Llist(k)), ', Ensemble ', num2str(i)])
        E1ens = zeros(size(T));
        E2ens = zeros(size(T));
        E3ens = zeros(size(T));
        E4ens = zeros(size(T));
        m1ens = zeros(size(T));
        m2ens = zeros(size(T));
        m3ens = zeros(size(T));
        m4ens = zeros(size(T));
        m1Eens = zeros(size(T));
        m2Eens = zeros(size(T));
        m3Eens = zeros(size(T));
        m4Eens = zeros(size(T));
        for j = 1:length(T)
            betaj = 1/T(j);
            mDb = -(betaj - beta0);
            dlog = logsumexp(Etime(:,i),mDb,onesVector);
            E1ens(j) = exp(logsumexp(Etime(:,i),mDb,Etime(:,i))    - dlog);
            E2ens(j) = exp(logsumexp(Etime(:,i),mDb,Etime(:,i).^2) - dlog);
            E3ens(j) = exp(logsumexp(Etime(:,i),mDb,Etime(:,i).^3) - dlog);
            E4ens(j) = exp(logsumexp(Etime(:,i),mDb,Etime(:,i).^4) - dlog);
            m1ens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i))    - dlog);
            m2ens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i).^2) - dlog);
            m3ens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i).^3) - dlog);
            m4ens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i).^4) - dlog);
            m1Eens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i)    .*Etime(:,i)) - dlog);
            m2Eens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i).^2 .*Etime(:,i)) - dlog);
            m3Eens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i).^3 .*Etime(:,i)) - dlog);
            m4Eens(j) = exp(logsumexp(Etime(:,i),mDb,mtime(:,i).^4 .*Etime(:,i)) - dlog);
        end
        E1L = E1L + E1ens;
        E2L = E2L + E2ens;
        E3L = E3L + E3ens;
        E4L = E4L + E4ens;
        m1L = m1L + m1ens;
        m2L = m2L + m2ens;
        m3L = m3L + m3ens;
        m4L = m4L + m4ens;
        m1EL = m1EL + m1Eens;
        m2EL = m2EL + m2Eens;
        m3EL = m3EL + m3Eens;
        m4EL = m4EL + m4Eens;

        E1TEnsL(:,i,k) = E1ens;
        E2TEnsL(:,i,k) = E2ens;
        E3TEnsL(:,i,k) = E3ens;
        E4TEnsL(:,i,k) = E4ens;
        m1TEnsL(:,i,k) = m1ens;
        m2TEnsL(:,i,k) = m2ens;
        m3TEnsL(:,i,k) = m3ens;
        m4TEnsL(:,i,k) = m4ens;
        m1ETEnsL(:,i,k) = m1Eens;
        m2ETEnsL(:,i,k) = m2Eens;
        m3ETEnsL(:,i,k) = m3Eens;
        m4ETEnsL(:,i,k) = m4Eens;
    end
    toc;
end

Etime = 0;
mtime = 0;

save([save_address '/', 'variables.mat'])



%%---------------------------------------------------------------------------------
%% make params Ensemble
Cpp = zeros(size(E1TEnsL));
Xpp = zeros(size(E1TEnsL));
Binder2 = zeros(size(E1TEnsL));
Binder4 = zeros(size(E1TEnsL));
BinderE4 = zeros(size(E1TEnsL));
Binder4Der = zeros(size(E1TEnsL));
Binder2Der = zeros(size(E1TEnsL));
mDer = zeros(size(E1TEnsL));
DK2 = zeros(size(E1TEnsL));
DK3 = zeros(size(E1TEnsL));
DK4 = zeros(size(E1TEnsL));
lnm1Der = zeros(size(E1TEnsL));
lnm2Der = zeros(size(E1TEnsL));
lnm3Der = zeros(size(E1TEnsL));
lnm4Der = zeros(size(E1TEnsL));
V = zeros(size(E1TEnsL,1),size(E1TEnsL,2),6);
U3 = zeros(size(E1TEnsL));
U4 = zeros(size(E1TEnsL));
Epp = zeros(size(E1TEnsL));
mpp = zeros(size(E1TEnsL));

beta = 1 ./ T;

for Ti = 1:length(T)
    for Ensi = 1:Nensemble
        for Li=1:length(Llist)
            E1 = E1TEnsL(Ti,Ensi,Li) / Nlist(Li);
            E2 = E2TEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 2;
            E3 = E3TEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 3;
            E4 = E4TEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 4;
            m1 = m1TEnsL(Ti,Ensi,Li) / Nlist(Li);
            m2 = m2TEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 2;
            m3 = m3TEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 3;
            m4 = m4TEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 4;
            m1E = m1ETEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 2;
            m2E = m2ETEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 3;
            m3E = m3ETEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 4;
            m4E = m4ETEnsL(Ti,Ensi,Li) / Nlist(Li) ^ 5;
            Temp = T(Ti);
            InvTemp = 1/ Temp;
            
            Cpp(Ti,Ensi,Li) = T(Ti)^(-2) * (E2-E1^2) * Nlist(Li);
            Xpp(Ti,Ensi,Li) = T(Ti)^(-1) * (m2-m1^2) * Nlist(Li);
            Binder2(Ti,Ensi,Li) = 1 - m2/(3 * m1^2);
            Binder4(Ti,Ensi,Li) = 1 - m4/(3 * m2^2);
            BinderE4(Ti,Ensi,Li) = 1 - (E4 -4*E1*E3+ 6*E1^2 * E2 - 3*E1^4)/ ...
                                   (3*(E2 - E1^2)^2);
            % Binder4Der(Ti,Ensi,Li) = (-1/(3*m2^2) * (m4E - m4*E1) + ...
            %     (2*m4E)/(3*m2^3) * (m2E - m2 * E1)) * (-1* T(Ti)^(-2)) ;
            % Binder4Der(Ti,Ensi,Li) = (1 - Binder4(Ti,Ensi,Li)) * ...
            %     (E1 - 2*m2E/m2 + m4E/m4)* Nlist(Li);
            m1der = -(m1E - m1 * E1) * Nlist(Li);
            m2der = -(m2E - m2 * E1) * Nlist(Li);
            m3der = -(m3E - m3 * E1) * Nlist(Li);
            m4der = -(m4E - m4 * E1) * Nlist(Li);
            mDer(Ti,Ensi,Li) = m1der;
            DK2(Ti,Ensi,Li) =  (m2der - 2*m1*m1der);
            DK3(Ti,Ensi,Li) =  (m3der - 3*(m2der*m1 + m2*m1der) ...
                +6*m1der*(m1^2));
            DK4(Ti,Ensi,Li) =  m4der -4*(m3der*m1 + m3*m1der) - 3*m2*m2der ...
                + 12*(m2der*m1^2 + m2*2*m1*m1der) - 6*4*(m1^3)*m1der;
            Binder4Der(Ti,Ensi,Li) = Nlist(Li) * (1 - Binder4(Ti,Ensi,Li)) * ...
                                     (E1 - 2*m2E/m2 + m4E/m4);
            Binder2Der(Ti,Ensi,Li) = Nlist(Li) * (1 - Binder2(Ti,Ensi,Li)) * ...
                                     (E1 - 2*m1E/m1 + m2E/m2);
            Epp(Ti,Ensi,Li) = (E1-Ebase/Nlist(Li));

            lnm1Der(Ti,Ensi,Li) = -(m1E/m1 - E1) * Nlist(Li);
            lnm2Der(Ti,Ensi,Li) = -(m2E/m2 - E1) * Nlist(Li);
            lnm3Der(Ti,Ensi,Li) = -(m3E/m3 - E1) * Nlist(Li);
            lnm4Der(Ti,Ensi,Li) = -(m4E/m4 - E1) * Nlist(Li);
            mpp(Ti,Ensi,Li) = m1;
            bm1 = log(m1der);
            bm2 = log(m2der);
            bm3 = log(m3der);
            bm4 = log(m4der);
            V(Ti,Ensi,Li,1) = 4*bm3 - 3*bm4;
            V(Ti,Ensi,Li,2) = 2*bm2 - bm4;
            V(Ti,Ensi,Li,3) = 3*bm2 - 2*bm3;
            V(Ti,Ensi,Li,4) = (4*bm1 - bm4)/3;
            V(Ti,Ensi,Li,5) = (3*bm1 - bm3)/2;
            V(Ti,Ensi,Li,6) = 2*bm1 - bm2;

            U3(Ti,Ensi,Li) = (m3-3*m2*m1 +2*m1^3) /(m1*(m2 - m1^2));
            U4(Ti,Ensi,Li) = (m4 - 4*m3*m1 - 3*m2^2 + 12*m2*(m1^2) - 6*(m1)^4) / (3*m2^2);
            
        end
    end
end

CppMean = squeeze(mean(Cpp,2));
XppMean = squeeze(mean(Xpp,2));
Binder2Mean = squeeze(mean(Binder2,2));
Binder4Mean = squeeze(mean(Binder4,2));
BinderE4Mean = squeeze(mean(BinderE4,2));
Binder4DerMean = squeeze(mean(Binder4Der,2));
Binder2DerMean = squeeze(mean(Binder2Der,2));
mDerMean = squeeze(mean(mDer,2));
DK2Mean = squeeze(mean(DK2,2));
DK3Mean = squeeze(mean(DK3,2));
DK4Mean = squeeze(mean(DK4,2));
lnm1DerMean = squeeze(mean(lnm1Der,2));
lnm2DerMean = squeeze(mean(lnm2Der,2));
lnm3DerMean = squeeze(mean(lnm3Der,2));
lnm4DerMean = squeeze(mean(lnm4Der,2));
VMean = squeeze(mean(V,2));
U3Mean = squeeze(mean(U3,2));
U4Mean = squeeze(mean(U4,2));
EppMean = squeeze(mean(Epp,2));
mppMean = squeeze(mean(mpp,2));

CppErr = squeeze(std(Cpp,0,2)) / sqrt(Nensemble);
XppErr = squeeze(std(Xpp,0,2)) / sqrt(Nensemble);
Binder2Err = squeeze(std(Binder2,0,2)) / sqrt(Nensemble);
Binder4Err = squeeze(std(Binder4,0,2)) / sqrt(Nensemble);
BinderE4Err = squeeze(std(BinderE4,0,2)) / sqrt(Nensemble);
Binder4DerErr = squeeze(std(Binder4Der,0,2)) / sqrt(Nensemble);
Binder2DerErr = squeeze(std(Binder2Der,0,2)) / sqrt(Nensemble);
mDerErr = squeeze(std(mDer,0,2)) / sqrt(Nensemble);
DK2Err = squeeze(std(DK2,0,2)) / sqrt(Nensemble);
DK3Err = squeeze(std(DK3,0,2)) / sqrt(Nensemble);
DK4Err = squeeze(std(DK4,0,2)) / sqrt(Nensemble);
lnm1DerErr = squeeze(std(lnm1Der,0,2)) / sqrt(Nensemble);
lnm2DerErr = squeeze(std(lnm2Der,0,2)) / sqrt(Nensemble);
lnm3DerErr = squeeze(std(lnm3Der,0,2)) / sqrt(Nensemble);
lnm4DerErr = squeeze(std(lnm4Der,0,2)) / sqrt(Nensemble);
VErr = squeeze(std(V,0,2)) / sqrt(Nensemble);
U3Err = squeeze(std(U3,0,2)) / sqrt(Nensemble);
U4Err = squeeze(std(U4,0,2)) / sqrt(Nensemble);
EppErr = squeeze(std(Epp,0,2)) / sqrt(Nensemble);
mppErr = squeeze(std(mpp,0,2)) / sqrt(Nensemble);

% for Li=1:length(Llist)
%     ValidDeltaT(Li,1) = Llist(Li);
%     ValidDeltaT(Li,2) = T0L(Li);
%     % ValidHalfRange = 1*T0L(Li)/sqrt(data(4,Li)*Nlist(Li)); %% gives wrong results when merging ensembles
%     % ValidDeltaT(Li,3) = ValidHalfRange;
%     ValidDeltaT(Li,4) = T0L(Li) - ValidHalfRange;
%     ValidDeltaT(Li,5) = T0L(Li) + ValidHalfRange;
% end

mppMean(1)
%%---------------------------------------------------------------------------------
%% E

clf('reset')
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
           'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
for i = 1:N
%     plot(T, Epp(:,i), 'Marker', '.','MarkerSize', 15);
    errorbar(T, EppMean(:,i), EppErr(:,i), 'Marker', markers{i},'MarkerSize', 15);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel '$E/ \Lambda$'
xlabel(myxLabel)
legendStrings = "L = " + string(Llist);
legend(legendStrings, 'Location','southeast','Interpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = jet(N+1);
cmap = turbo(N+1);
% cmap = linspecer(N);
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'energy.png']) 
% export_fig   %% presentation/article quality figures. alternatively use exportgraphics

%%---------------------------------------------------------------------------------
%%  m
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
%     plot(T, mpp(:,i), 'Marker', '.','MarkerSize', 15);
    errorbar(T, mppMean(:,i), mppErr(:,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel(myxLabel)
ylabel '$|\textbf{m}|$'
legendStrings = "L = " + string(Llist);
legend(legendStrings, 'Location','southwest','Interpreter','latex','fontsize',8)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
cmap = jet(N+1);
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'm.png']);

%%---------------------------------------------------------------------------------
%% C
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
SizeRange = 1:N;
for i = SizeRange
%     plot(T, Cpp(:,i), 'Marker', '.','MarkerSize', 15);
    errorbar(T, CppMean(:,i), CppErr(:,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
xlabel(myxLabel);
ylabel '$c$'
legendStrings = "L = " + string(Llist(SizeRange));
legend(legendStrings,'Interpreter','latex');
set(gca,'fontsize',12) ;
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
cmap = jet(N+1);
ax.ColorOrder = cmap;
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',18)  
% axes('Position',[.65 .65 .25 .25]) % begin inset
% box on
% scatter(T_E_m_c_X(:,1), T_E_m_c_X(:,2), 15, 'filled');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',14) 
% ylabel '$E$'
% xlabel '$T$'
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'C.png']) % end inset

%%---------------------------------------------------------------------------------
%% chi
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
SizeRange = 1:N;
for i = SizeRange
    % plot(T, XppMean(:,i), 'Marker', '.','MarkerSize', 15);
    errorbar(T, XppMean(:,i), XppErr(:,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
xlabel(myxLabel)
ylabel '$\chi$'
legendStrings = "L = " + string(Llist(SizeRange))
legend(legendStrings,'Interpreter','latex','fontsize',8)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
cmap = jet(N+1);
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'X.png'])
%%---------------------------------------------------------------------------------
%% Binder
% clf('reset')


SizeRange = 18:1:21;
Trange = 20:5:110;

% fig = figure;
% fig.Position  = [0 0 300 200];
colordiff=0;
cmap = linspecer(length(SizeRange));
markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
                'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
% markerShapes = { 's', 'o', '^', 'v', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
mi=0;
% axes(axs(1))
for i = SizeRange
    mi = mi + 1;
%     SplineVal = spline(orderParam(:,:,i), data(:,6,i),xRange);
%     errorbar(orderParam(:,:,i), data(:,6,i), error(:,6,i), '.', xRange, SplineVal);
%     f = fit(orderParam(:,:,i),data(:,6,i),'smoothingspline');
    h = errorbar(T(Trange), BinderE4Mean(Trange,i), BinderE4Err(Trange,i), 'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
         'Color',cmap(mi+colordiff,:), 'LineWidth', 0.8);
    % h = plot(T(Trange), Binder4Mean(Trange,i));
    % h = plot(T(Trange), Binder4Mean(Trange,i), 'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
    %      'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h.MarkerFaceColor = h.Color;
    % plot(T(Trange), Binder4Mean(Trange,i), 'Marker', markers{i},'MarkerSize', 5);
    hold on
    
end

xlim([0.845 0.853])
% xlim([0.845 0.853])
yl = ylim;
xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/4) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/4) :xl(2));
xticks([0.845, 0.849, 0.853])
% axis tight
xlabel(myxLabel)
ylabel '$U_4$'

myFontSize =13;
set(gca,'fontsize',myFontSize) 

legendStrings = "$L = $" + string(Llist(SizeRange))
lgd = legend(legendStrings, 'Location','southwest','Interpreter','latex', 'FontSize', myFontSize-1);
set(lgd,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;

set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% pbaspect([2 1 1])

% axes('Position',[0.2 .65 .25 .25]) % begin inset
% box on
% for i = 1:N
%     plot(orderParam(:,:,i), data(:,7,i), '--.','MarkerSize', 9);
%     hold on
% end
% xlim([0.665 0.675])
% ylim([1.04 1.2])
% % ylim([1 1.2])
% % xticks(0.64:0.02:0.68)
% % xlabel '$T$'
% % ylabel '$Binder parameter$'
% ax = gca;
% ax.TickLength(1) = 0.04;
% % ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.25:ax.YAxis.Limits(2);
% % ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):0.5:ax.XAxis.Limits(2);
% set(gca,'XMinorTick','on','YMinorTick','on')
% % set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',9)  % end inset
% set(gca,'TickLabelInterpreter','latex');
% set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'Binder.png'])
exportgraphics(gcf,[save_address '/',  'Binder' '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  'Binder' '.eps'])
%%---------------------------------------------------------------------------------
%% Binder Derivative

clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
% SizeRange = 1:N;
SizeRange = 14;
for i = SizeRange
    % errorbar(T, Binder2DerMean(:,i), Binder2DerErr(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, Binder4DerMean(:,i), Binder4DerErr(:,i), 'Marker', '.','MarkerSize', 15);
    % plot(T, Binder2DerMean(:,i), 'Marker', '.','MarkerSize', 15);
    errorbar(T, DK3Mean(:,i), DK3Err(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, DK2Mean(:,i), DK2Err(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(beta, lnm1DerMean(:,i), lnm1DerErr(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, lnm2DerMean(:,i), lnm2DerErr(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, lnm3DerMean(:,i), lnm3DerErr(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, lnm4DerMean(:,i), lnm4DerErr(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, mDerMean(:,i), mDerErr(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, U3Mean(:,i), U3Err(:,i), 'Marker', '.','MarkerSize', 15);
    % errorbar(T, U4Mean(:,i), U4Err(:,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
xlabel(myxLabel)
ylabel 'Binder Derivative'
legendStrings = "L = " + string(Llist(SizeRange))

set(gca,'fontsize',13) 
legend(legendStrings,'Interpreter','latex', 'FontSize', myFontSize-1)
set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'BinderDeriv.png'])
%%---------------------------------------------------------------------------------
%% FSS

% Tc = data(1);
Tc = 0.8495;
reducedTemp = T / Tc -1;
L = sqrt(Nlist);
gamma = 1.65;
nu = 1/(1.023);
beta = 0.15;
% Binder = data(:,7,:);
Xtilde =zeros(size(Xpp));
x =zeros(size(Xpp));
SizeRange = 3:8;
for i=SizeRange 
%     Xtilde(:,i) = L(i)^(-gamma/nu) * Xipp(:,i);
    Xtilde(:,i) = L(i)^(beta/nu) * mpp(:,i);
    x(:,i) = L(i)^(1/nu) * reducedTemp;
end

clf('reset')
for i = SizeRange
    plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 10);
%     plot(x(:,i), Binder(:,i), '.','MarkerSize', 10);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
xlim([-0.5 0.5])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel '$E/ \Lambda$'
xlabel(myxLabel)
legendStrings = "N = " + string(Nlist);
% legend(legendStrings, 'Location','southeast','Interpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

%%---------------------------------------------------------------------------------
%% chi max vs L, power fit
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
SizeRange = 1:21;
% InvNu = 0.98;
sizesInvExpNu = Llist(SizeRange);
maxidxL = 0;
% nuInv = 1/0.95338
for i = 1:length(sizesInvExpNu)
%     TT = 0;
    % chidiff = 0;
    % [maxchi, maxidx] = max(XppMean(:,i));
    [maxchi, maxidx] = max(lnm1DerMean(:,i));
    % [maxchi, maxidx] = max(Binder4DerMean(:,i));
    % [maxchi, maxidx] = max(mDerMean(:,i));
    maxidxL(i) = maxchi;
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizesInvExpNu = sizesInvExpNu';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
plot(sizesInvExpNu, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(min(sizesInvExpNu) - min(sizesInvExpNu)/10, max(sizesInvExpNu) + max(sizesInvExpNu)/10, 500);
p = polyfit(log(sizesInvExpNu),log(maxidxL),1);
yFit = polyval(p, log(xFit));
plot(xFit, exp(yFit) ,'r' , 'LineWidth', 1 );

mdl = fitlm((sizesInvExpNu),log(maxidxL));
R2 = mdl.Rsquared.Ordinary
Error = mdl.Coefficients.SE(1)

% 
% xl = xlim;
% yl = ylim;
% xt = 0.2 * (xl(2)-xl(1)) + xl(1);
% yt = 0.2 * (yl(2)-yl(1)) + yl(1);
% caption = sprintf('$\\chi \\sim  L^{%.3f}$', p(1));
% h = text(xt, yt, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
% set(h,'Rotation',30);

axis padded
xlabel '$L$'
ylabel '$\chi_{\mathrm{max}}$'
ylim([1 (max(maxidxL)+ max(maxidxL)/3)])
set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exponent = p(1)
% exponent = p(1)
saveas(gcf, [save_address '/', 'XvsN power fit.png'])

%%---------------------------------------------------------------------------------
%% Tempereture range for valid chi max
TchiMax = zeros(length(Llist),8);
TchiMaxRelDistanceFromT0 = zeros(length(Llist),8);
TchiMax(:,1) = Llist;
TchiMaxRelDistanceFromT0(:,1) = Llist;

[~, maxidx] = max(XppMean,[],1);
TchiMax(:,2) = T(maxidx);
TchiMaxRelDistanceFromT0(:,2) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

[~, maxidx] = max(lnm1DerMean,[],1);
TchiMax(:,3) = T(maxidx);
TchiMaxRelDistanceFromT0(:,3) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

[~, maxidx] = max(lnm2DerMean,[],1);
TchiMax(:,4) = T(maxidx);
TchiMaxRelDistanceFromT0(:,4) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

[~, maxidx] = max(lnm3DerMean,[],1);
TchiMax(:,5) = T(maxidx);
TchiMaxRelDistanceFromT0(:,5) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

[~, maxidx] = max(lnm4DerMean,[],1);
TchiMax(:,6) = T(maxidx);
TchiMaxRelDistanceFromT0(:,6) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

[~, maxidx] = max(mDerMean,[],1);
TchiMax(:,7) = T(maxidx);
TchiMaxRelDistanceFromT0(:,7) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

[~, maxidx] = max(CppMean,[],1);
TchiMax(:,8) = T(maxidx);
TchiMaxRelDistanceFromT0(:,8) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

[~, maxidx] = min(DK2Mean,[],1);
TchiMax(:,9) = T(maxidx);
TchiMaxRelDistanceFromT0(:,9) = (T(maxidx) - T0L) ./ValidDeltaT(:,3)';

disp("done");

%%---------------------------------------------------------------------------------
%% find critical point by Tc (L) plot
% TExtreames = [min(equal_within_error), (max(equal_within_error))];
% InvNuExtreames = [InvNuMean(min(equal_within_error)), InvNuMean(max(equal_within_error))];
Tindices = 1:length(T); 
% Tindices = 40:155;
% Tindices = 1:200; 
TRange =  T(Tindices);

index = 1;
SizeRange = 10:N;
InvNu = 1/0.96;
% InvNu = InvNuExtreames(index);
LExpMinuxInvNu = Llist(SizeRange) .^ (-InvNu);
sizesInvExpNu = LExpMinuxInvNu;
maxidxL = [];
Tmax = [];
TmaxErr = []; 
j = 1;
for i = SizeRange
%     TT = 0;
% %     chidiff = 0;
    % [~, maxidx] = max(XppMean(Tindices,i));
    % [~, maxidx] = max(lnm2DerMean(Tindices,i));
    % [~, maxidx] = max(DK2Mean(:,i));
    % [~, maxidx] = max(DK3Mean(:,i));
    % [~, maxidx] = max(DK4Mean(Tindices,i));
    % [~, maxidx] = max(U4Mean(Tindices,i));
    % [~, maxidx] = min(U3Mean(Tindices,i));
    % [~, maxidx] = max(mDerMean(Tindices,i));
    % [~, maxidx] = max(CppMean(:,i));
    % [~, maxidx] = findpeaks(Binder4DerMean(Tindices,i)); 
    % maxidx = maxidx(1);
    % maxidxL(j) = TRange(maxidx);

    % x_xerr = xAtMaxWithError(T, XppMean(:,i), XppErr(:,i));
    % x_xerr = xAtMaxWithError(T, CppMean(:,i), CppErr(:,i));
    % x_xerr = xAtMaxWithError(T, DK2Mean(:,i), DK2Err(:,i));
    % x_xerr = xAtMaxWithError(T, -DK2Mean(:,i), DK2Err(:,i));
    % x_xerr = xAtMaxWithError(T, -U3Mean(:,i), U3Err(:,i));

    % x_xerr = xAtMaxWithError3(T, DK2(:, :, i));
    % x_xerr = xAtMaxWithError3(T, -DK2(:, :, i));
    % x_xerr = xAtMaxWithError3(T, DK3(:, :, i));
    % x_xerr = xAtMaxWithError3(T, -DK3(:, :, i));
    % x_xerr = xAtMaxWithError3(T, DK4(:, :, i));
    % x_xerr = xAtMaxWithError3(T, -U3(:, :, i));
    % x_xerr = xAtMaxWithError3(T, U4(:, :, i));

    % x_xerr = xAtMaxWithError3(T, Xpp(:, :, i));
    % x_xerr = xAtMaxWithError3(T, Cpp(:, :, i));
    % x_xerr = xAtMaxWithError3(T, mDer(:, :, i));
    % x_xerr = xAtMaxWithError3(T, lnm4Der(:, :, i));
    x_xerr = xAtMaxWithError3(T, Binder4Der(:, :, i));
    


    Tmax(j) = x_xerr(1);
    TmaxErr(j) = x_xerr(2);
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
    j = j + 1;
end
sizesInvExpNu = sizesInvExpNu';
maxidxL = maxidxL';
maxidxL = Tmax;
% f = fit(sizes,maxidxL,'exp1')
errorbar(sizesInvExpNu, Tmax, TmaxErr, 's', 'Marker', '.','MarkerSize', 15);
% plot(sizesInvExpNu, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(0, max(sizesInvExpNu) + max(sizesInvExpNu)/10, 500);
p = polyfit((sizesInvExpNu),(maxidxL),1);
yFit = polyval(p, (xFit));
plot(xFit, (yFit) ,'r' , 'LineWidth', 1 );

mdl = fitlm((sizesInvExpNu),(maxidxL));
R2 = mdl.Rsquared.Ordinary
Error = mdl.Coefficients.SE(1)


% xl = xlim;
% yl = ylim;
% xt = 0.2 * (xl(2)-xl(1)) + xl(1);
% yt = 0.2 * (yl(2)-yl(1)) + yl(1);
% caption = sprintf('$\\chi \\sim  L^{%.3f}$', p(1));
% h = text(xt, yt, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
% set(h,'Rotation',30);

% axis padded
xlabel '$L^{-\nu}$'
ylabel '$\chi_{\mathrm{max}}$'
% ylim([1 (max(maxidxL)+ max(maxidxL)/3)])
% set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
disp(['Tc = ', num2str(polyval(p, 0))]);
exponent = p(1) * InvNu
% exponent = p(1)
% saveas(gcf, [save_address '/', 'Tc by L Exp -InvNu fit.png'])


%%---------------------------------------------------------------------------------
%% find critical point by Tc (L) plot multiple
% TExtreames = [min(equal_within_error), (max(equal_within_error))];
% InvNuExtreames = [InvNuMean(min(equal_within_error)), InvNuMean(max(equal_within_error))];
Tindices = 1:length(T); 
% Tindices = 40:155;
% Tindices = 1:200; 
TRange =  T(Tindices);

markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
                'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
% markerShapes = { 's', 'o', '^', 'v', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
mi=0;
cmap = linspecer(8);

index = 1;
SizeRange = 5:N;
InvNu = 0.95;
% InvNu = InvNuExtreames(index);
LExpMinuxInvNu = Llist(SizeRange) .^ (-InvNu);
sizesInvExpNu = LExpMinuxInvNu;
maxidxL = [];
Tmax = [];
TmaxErr = []; 
x_xerr = [];
j = 1;
for i = SizeRange
    k = 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, DK2(:, :, i)); k = k + 1;
    x_xerr(j, :, k)  = xAtMaxWithError3(T, -DK2(:, :, i)); k = k + 1;
    % x_xerr(j, :, k)  = xAtMaxWithError3(T, DK3(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, -DK3(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, DK4(:, :, i)); k = k + 1;
    x_xerr(j, :, k) = xAtMaxWithError3(T, -U3(:, :, i));  k = k + 1;
    x_xerr(j, :, k) = xAtMaxWithError3(T, U4(:, :, i)); k = k + 1;

    % x_xerr(j, :, k)= xAtMaxWithError3(T, Xpp(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, Cpp(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, mDer(:, :, i)); k = k + 1;

    % x_xerr(j, :, k) = xAtMaxWithError3(T, lnm1Der(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, Binder4Der(:, :, i)); k = k + 1;

    j = j + 1;
end
sizesInvExpNu = sizesInvExpNu';
maxidxL = maxidxL';
maxidxL = Tmax;
Tmax = squeeze(x_xerr(:,1,:));
TmaxErr = squeeze(x_xerr(:,2,:));


for k = 1:size(Tmax,2)
    mi = mi + 1;
    h=errorbar(sizesInvExpNu, Tmax(:, k), TmaxErr(:, k), 'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
     'Color',cmap(mi,:), 'LineWidth', 0.5, 'LineStyle','none');
    h.MarkerFaceColor = h.Color;
    hold on

    maxidxL = Tmax(:, k);
    xFit = linspace(0, max(sizesInvExpNu) + max(sizesInvExpNu)/10, 500);
    p = polyfit((sizesInvExpNu),(maxidxL),1);
    yFit = polyval(p, (xFit));
    plot(xFit, (yFit) ,'Color',cmap(mi,:) , 'LineWidth', 0.8, 'HandleVisibility', 'off');
end



% index = 1;
SizeRange = 8:N;
% InvNu = 1/0.96;
% InvNu = InvNuExtreames(index);
LExpMinuxInvNu = Llist(SizeRange) .^ (-InvNu);
sizesInvExpNu = LExpMinuxInvNu;
maxidxL = [];
Tmax = [];
TmaxErr = []; 
x_xerr = [];
j = 1;
for i = SizeRange
    k = 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, DK2(:, :, i)); k = k + 1;
    % x_xerr(j, :, k)  = xAtMaxWithError3(T, -DK2(:, :, i)); k = k + 1;
    x_xerr(j, :, k) = xAtMaxWithError3(T, DK3(:, :, i)); k = k + 1;
    x_xerr(j, :, k) = xAtMaxWithError3(T, -DK3(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, DK4(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, -U3(:, :, i));  k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, U4(:, :, i)); k = k + 1;

    x_xerr(j, :, k)= xAtMaxWithError3(T, Xpp(:, :, i)); k = k + 1;
    x_xerr(j, :, k) = xAtMaxWithError3(T, Cpp(:, :, i)); k = k + 1;
    x_xerr(j, :, k) = xAtMaxWithError3(T, mDer(:, :, i)); k = k + 1;

    % x_xerr(j, :, k) = xAtMaxWithError3(T, lnm1Der(:, :, i)); k = k + 1;
    % x_xerr(j, :, k) = xAtMaxWithError3(T, Binder4Der(:, :, i)); k = k + 1;

    j = j + 1;
end
sizesInvExpNu = sizesInvExpNu';
maxidxL = maxidxL';
maxidxL = Tmax;
Tmax = squeeze(x_xerr(:,1,:));
TmaxErr = squeeze(x_xerr(:,2,:));

for k = 1:size(Tmax,2)
    mi = mi + 1;
    errorbar(sizesInvExpNu, Tmax(:, k), TmaxErr(:, k), 'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
     'Color',cmap(mi,:), 'LineWidth', 0.8, 'LineStyle','none');
    % h.MarkerFaceColor = h.Color;
    hold on
    maxidxL = Tmax(:, k);
    xFit = linspace(0, max(sizesInvExpNu) + max(sizesInvExpNu)/10, 500);
    p = polyfit((sizesInvExpNu),(maxidxL),1);
    yFit = polyval(p, (xFit));
    plot(xFit, (yFit) ,'Color',cmap(mi,:) , 'LineWidth',  0.8, 'HandleVisibility', 'off');
end


legendStrings = ["$D_{K_2}(\mathrm{min})$", "$P_3$", "$P_4$", "$D_{K_3}(\mathrm{max})$", ...
                 "$D_{K_3}(\mathrm{min})$", "$\chi$", "$c$", "$dm/d\beta$"];
legend(legendStrings, 'Location','best', ...
    'Interpreter','latex')

% xl = xlim;
% yl = ylim;
% xt = 0.2 * (xl(2)-xl(1)) + xl(1);
% yt = 0.2 * (yl(2)-yl(1)) + yl(1);
% caption = sprintf('$\\chi \\sim  L^{%.3f}$', p(1));
% h = text(xt, yt, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
% set(h,'Rotation',30);

% axis padded
xlabel '$L^{-1/\nu}$'
ylabel '$T_{\mathrm{c}}(L)$'
% ylim([1 (max(maxidxL)+ max(maxidxL)/3)])
% set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% disp(['Tc = ', num2str(polyval(p, 0))]);
% exponent = p(1) * InvNu
% exponent = p(1)
saveas(gcf, [save_address '/', 'find critical point by Tc (L).png'])


%%---------------------------------------------------------------------------------
%% find critical point by Tc (L) multiple
% TExtreames = [min(equal_within_error), (max(equal_within_error))];
% InvNuExtreames = [InvNuMean(min(equal_within_error)), InvNuMean(max(equal_within_error))];

SizeRange1 = 5:N;
SizeRange2 = 8:N;
InvNu = 0.97;
 
Nboot = 100;
Nensemble = size(V, 2);
Tc = [];
for iboot = 1:Nboot

    BootEnsembles = randi(Nensemble, Nensemble, 1);
    BootEnsembles = BootEnsembles';

    SizeRange = SizeRange1;
    sizesInvExpNu = Llist(SizeRange) .^ (-InvNu);
    xAtMax = [];
    j = 1;
    mi = 0;
    for i = SizeRange
        k = 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, -DK2(:, :, i), BootEnsembles); k = k + 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, -U3(:, :, i), BootEnsembles); k = k + 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, U4(:, :, i), BootEnsembles); k = k + 1;
        j = j + 1;
    end
    Tmax = xAtMax;
    for k = 1:size(Tmax,2)
        mi = mi + 1;
        p = polyfit((sizesInvExpNu),(Tmax(:, k)),1);
        Tc(mi, iboot) = polyval(p, 0);
    end

    % ---------- second part ------------
    SizeRange = SizeRange2;
    sizesInvExpNu = Llist(SizeRange) .^ (-InvNu);
    xAtMax = [];
    j = 1;
    for i = SizeRange
        k = 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, DK3(:, :, i), BootEnsembles); k = k + 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, -DK3(:, :, i), BootEnsembles); k = k + 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, Xpp(:, :, i), BootEnsembles); k = k + 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, Cpp(:, :, i), BootEnsembles); k = k + 1;
        xAtMax(j , k)  = xAtMaxInBootEns(T, mDer(:, :, i), BootEnsembles); k = k + 1;
        j = j + 1;
    end
    Tmax = xAtMax;
    for k = 1:size(Tmax,2)
        mi = mi + 1;
        p = polyfit((sizesInvExpNu),(Tmax(:, k)),1);
        Tc(mi, iboot) = polyval(p, 0);
    end
end
meanTc = mean(Tc,2)
stdTc = std(Tc,0,2)
varTc = var(Tc,0,2);
InvSigma2 = sum(1./varTc);
Sigma2 = 1 / InvSigma2;
UnWeightedAverageTcAll = mean(meanTc)
stdUnWeightedAverageTcAll = std(meanTc)
WeightedAverageTcAll = sum(meanTc ./ varTc) * Sigma2
stdWeightedAverageTcAll = sqrt(Sigma2)
%%---------------------------------------------------------------------------------
%% chi vs L at Tc, power fit
clf('reset')
% TExtreames = [min(equal_within_error), (max(equal_within_error))];
% InvNuExtreames = [InvNuMean(min(equal_within_error)), InvNuMean(max(equal_within_error))];

% TcIndex = TExtreames(1);
TcIndex = 69;
InvNu = InvNuMean(TcIndex);


% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
SizeRange = 1:21;
sizesInvExpNu = Llist(SizeRange);
maxidxL = 0;
erridxL = 0;

% TcIndex = 40;
% InvNu = 0.99;
for i = 1:length(sizesInvExpNu)
%     TT = 0;
    chidiff = 0;
    % ChiChritical = XppMean(TcIndex,i);
    ChiChritical = mppMean(TcIndex,i);
    % ChiChritical = EppMean(TcIndex,i);
    % ChiChritical = CppMean(TcIndex,i);
    % ChiChritical = lnm1DerMean(TcIndex,i);
    % ChiChritical = mDerMean(TcIndex,i);
    % ChiChritical = Binder2DerMean(TcIndex,i);

%     [maxchi, maxidx] = max(lnm1Der(:,i));
%     [maxchi, maxidx] = max(mDer(:,i));
    maxidxL(i) = ChiChritical;
    erridxL(i) = mppErr(TcIndex,i);
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizesInvExpNu = sizesInvExpNu';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
% plot(sizesInvExpNu, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
errorbar(sizesInvExpNu, maxidxL, erridxL', 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(min(sizesInvExpNu) - min(sizesInvExpNu)/10, max(sizesInvExpNu) + max(sizesInvExpNu)/10, 500);
p = polyfit(log(sizesInvExpNu),log(maxidxL),1);
yFit = polyval(p, log(xFit));
plot(xFit, exp(yFit) ,'r' , 'LineWidth', 1 );

mdl = fitlm(log(sizesInvExpNu),log(maxidxL));
R2 = mdl.Rsquared.Ordinary

% axis padded
xlabel '$L$'
ylabel '$\chi$ at $T_c$'; caption = sprintf('$\\chi \\sim  L^{%.3f}$', p(1));
% ylabel '$m$ at $T_c$'; caption = sprintf('$m \\sim  L^{%.3f}$', p(1));
% ylim([1 (max(maxidxL)+ max(maxidxL)/3)])
set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 

xl = xlim;
yl = ylim;
xt = 0.1 * (xl(2)-xl(1)) + xl(1);
yt = 0.2 * (yl(2)-yl(1)) + yl(1);
h = text(xt, yt, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
set(h,'Rotation',30);

set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) ;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exponent = p(1) / InvNu
% exponent = p(1)

saveas(gcf, [save_address '/', 'XvsN power fit.png'])

%%---------------------------------------------------------------------------------
%% nu by V experiment
Vindices = 1:6;
SizeRange= 1:N;
Ti = 45;

fig = figure();
% clf('reset')
ax = axes(fig);
cmap = linspecer(length(Vindices));
% cmap = turbo(length(Vindices));
ax.ColorOrder = cmap;
mi = 1;
slopes =  [];
for Vindex = Vindices
    h = errorbar((Llist(SizeRange)), VMean(Ti,SizeRange,Vindex), VErr(Ti,SizeRange,Vindex),...
         'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 5, ...
         'LineStyle', 'none');
    % h = plot((Llist(SizeRange)), VMean(Ti,SizeRange,Vindex),...
    %  'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 6, ...
    %  'LineStyle', 'none');
    h.MarkerFaceColor = h.Color;
    hold on
    x = (Llist(SizeRange));
    y = VMean(Ti,SizeRange,Vindex);
    % Calculate the coefficients of the fit
    coeffs = polyfit(log(x), y, 1); % 1 for a linear fit
    % Generate a range of x values
    xFit = linspace(10, max(x) + max(x)/5, 1000);
    
    % Apply the coefficients to the xFit data
    yFit = polyval(coeffs, log(xFit));
    plot(xFit, yFit, '-', 'HandleVisibility', 'off','Color',cmap(mi,:)); 
    slopes(Vindex) = coeffs(1);
    mi = mi + 1;
end
disp('slopes = ')
disp(slopes')
mdl = fitlm(log(Llist(SizeRange)),VMean(Ti,SizeRange,Vindex));
R2 =  mdl.Rsquared.Ordinary
Temperature = T(Ti)

% ylim([-0.5 3.5])
% xlim([10 144])

set(gca, 'XScale', 'log')   %using plot
xlabel '$L$'
ylabel '$V_i$' 
title(['$T = $ ', num2str(T(Ti))])
legendStrings = "$V_" + string(Vindices) + "$";
legend(legendStrings,'Location','northwest','Interpreter','latex')
set(gca,'fontsize',14) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'V_i.png'])
saveas(gcf, [save_address '/', 'V_i at T = ', num2str(T(Ti)), '.png'])

%%---------------------------------------------------------------------------------
%% nu by V Ensemble
% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
Tindices = 1:length(T);
InvNuByV = zeros(6,Nensemble);
SizeRange= 1:19;
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);
InvNuByVErr = zeros(6,1);
InvNuByVMean = zeros(6,1);
k = 1;
for Ti = Tindices
    for Ensi = 1:Nensemble
        for i = 1:6
            [p,S] = polyfit(log(Llist(SizeRange)), V(Ti,Ensi,SizeRange,i),1);
            InvNuByV(i,Ensi) = p(1);
            
            % InvNuByV(i,Ensi) = mdl.Coefficients.Estimate(2);
        end
    end
    InvNuByVMean = squeeze(mean(InvNuByV,2));
    InvNuByVErr = squeeze(std(InvNuByV,0,2)) / sqrt(Nensemble);
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    k = k + 1;
end

[~, idx] = min(InvNuStd);
idxOfTc = Tindices(idx);
TCritical = T(idxOfTc);
InvNu = InvNuMean(idxOfTc);
nuInv = InvNu;
disp(['idxOfTc = ' num2str(idxOfTc)]);
disp(['TCritical = ' num2str(TCritical)]);
disp(['InvNu = ' num2str(InvNu)]);
disp(['nu = ' num2str(1/InvNu)]);
disp(['std = ' num2str(InvNuStd(idxOfTc))]);
fprintf('\n');


%%---------------------------------------------------------------------------------
%% nu by V Ensemble no Bins
% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
Tindices = 1:length(T);
InvNuByV = zeros(6);
SizeRange= 1:21;
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);
InvNuByVErr = zeros(6,1);
InvNuByVMean = zeros(6,1);
k = 1;
for Ti = Tindices
    for i = 1:6
        [p,S] = polyfit(log(Llist(SizeRange)), VMean(Ti,SizeRange,i),1);
        InvNuByV(i) = p(1);
        mdl = fitlm(log(Llist(SizeRange)), squeeze(V(Ti,Ensi,SizeRange,i)));
        InvNuByVMean(i) = mdl.Coefficients.Estimate(2);
        InvNuByVErr(i) = mdl.Coefficients.SE(2);
        
    end
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    k = k + 1;
end

[~, idx] = min(InvNuStd);
idxOfTc = Tindices(idx);
TCritical = T(idxOfTc);
InvNu = InvNuMean(idxOfTc);
nuInv = InvNu;
disp(['idxOfTc = ' num2str(idxOfTc)]);
disp(['TCritical = ' num2str(TCritical)]);
disp(['InvNu = ' num2str(InvNu)]);
disp(['nu = ' num2str(1/InvNu)]);
disp(['std = ' num2str(InvNuStd(idxOfTc))]);
fprintf('\n');

%%---------------------------------------------------------------------------------
%% nu by V Ensemble no Bins, Tc and nu for max L
% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
Tindices = 1:length(T);
LmaxRange = 6:N;
Lmax_Tc_InvNu = zeros(length(LmaxRange), 4);
Li = 1;
for Lmax = LmaxRange
SizeRange= 1:Lmax;
InvNuByV = zeros(6);
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);
InvNuErrMean = zeros(length(Tindices),1);
InvNuByVErr = zeros(6,1);
InvNuByVMean = zeros(6,1);
k = 1;
for Ti = Tindices
    for i = 1:6
        [p,S] = polyfit(log(Llist(SizeRange)), VMean(Ti,SizeRange,i),1);
        InvNuByV(i) = p(1);
        mdl = fitlm(log(Llist(SizeRange)), squeeze(V(Ti,Ensi,SizeRange,i)));
        InvNuByVMean(i) = mdl.Coefficients.Estimate(2);
        InvNuByVErr(i) = mdl.Coefficients.SE(2);
        
    end
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    InvNuErrMean(k) = mean(InvNuByVErr);
    k = k + 1;
end

[~, idx] = min(InvNuStd);
idxOfTc = Tindices(idx);
TCritical = T(idxOfTc);
InvNu = InvNuMean(idxOfTc);
nuInv = InvNu;
InVnuErr = InvNuErrMean(idxOfTc);
Lmax_Tc_InvNu(Li, :) = [Llist(Lmax), TCritical, InvNu, InVnuErr];
Li = Li + 1;
end
%%---------------------------------------------------------------------------------
%% Lmax vs Tc

plot(Lmax_Tc_InvNu(:,1), Lmax_Tc_InvNu(:,2), 'Marker', '.','MarkerSize', 15)

yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel("$L_{\mbox{max}}$")
ylabel '$T_c$'
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'Lmax vs Tc.png']) 

%%---------------------------------------------------------------------------------
%% Lmax vs nu

errorbar(Lmax_Tc_InvNu(:,1), Lmax_Tc_InvNu(:,3), Lmax_Tc_InvNu(:,4),...
    'Marker', '.','MarkerSize', 15)

yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel("$L_{\mbox{max}}$")
ylabel '$1/\nu$'
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'Lmax vs nu.png']) 
%%--------------------------------------------------------------------------------
%% Equal Nu By V Temperatue Range no Bins
Tindices = 1:length(T);
SizeRange = 1:19;

% Initialize array to keep track of temperature points where conditions are met
equal_within_error = [];

InvNuRange = [];
InvNuByV = zeros(6);
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);
InvNuByVErr = zeros(6,1);
InvNuByVMean = zeros(6,1);

k = 1;
for Ti = Tindices
    % Evaluate variables at this temperature
    for i = 1:6
        [p,S] = polyfit(log(Llist(SizeRange)), VMean(Ti,SizeRange,i),1);
        InvNuByV(i) = p(1);
        mdl = fitlm(log(Llist(SizeRange)), squeeze(V(Ti,Ensi,SizeRange,i)));
        InvNuByVMean(i) = mdl.Coefficients.Estimate(2);
        InvNuByVErr(i) = mdl.Coefficients.SE(2);
        
    end
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    k = k + 1;
    variable_values = InvNuByVMean;
    errors = InvNuByVErr;
    % Now check if all variables are within error bounds of each other
    equal = true; % Assume they are equal to start
    for i = 1:length(variable_values)
        for j = i+1:length(variable_values)
            if abs(variable_values(i) - variable_values(j)) > (errors(i) + errors(j))
                equal = false; % If any are not within error, they are not equal
                break; % No need to check further if one inequality is found
            end
        end
        if ~equal
            break; % Exit the outer loop as well
        end
    end
    
    % If all variables are equal within the error, add the temperature to the list
    if equal
        equal_within_error = [equal_within_error, Ti];
        InvNuRange = [InvNuRange, mean(InvNuByVMean)];
    end
end
% Display the temperature range where variables are equal within errors
if ~isempty(equal_within_error)
    fprintf('Temperature range = \n %f to %f.\n', ...
            T(min(equal_within_error)), T(max(equal_within_error)));
    fprintf('respective temperature indices = \n %i to %i\n', ...
            (min(equal_within_error)), (max(equal_within_error)));
    fprintf('Index range width = \n %i \n', ...
            max(equal_within_error) - min(equal_within_error));
    fprintf('InvNu range = \n %f to %f.\n', ...
            InvNuMean(min(equal_within_error)), InvNuMean(max(equal_within_error)))
else
    disp('No temperature range found where all variables are equal');
end
fprintf('\n')
%%--------------------------------------------------------------------------------
%% Equal Nu By V Temperatue Range (bins)
% Tindices = 1:length(T);
Tindices = 2:120;
TRange =  T(Tindices);
SizeRange = 1:21;

% Initialize array to keep track of temperature points where conditions are met
equal_within_error = [];

InvNuRange = [];
InvNuByV = zeros(6,Nensemble);
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);
InvNuByVErr = zeros(6,1);
InvNuByVMean = zeros(6,1);

k = 1;
for Ti = Tindices
    % Evaluate variables at this temperature
    for Ensi = 1:Nensemble
        for i = 1:6
            [p,S] = polyfit(log(Llist(SizeRange)), V(Ti,Ensi,SizeRange,i),1);
            InvNuByV(i,Ensi) = p(1);
            % InvNuByV(i,Ensi) = mdl.Coefficients.Estimate(2);
        end
    end
    InvNuByVMean = squeeze(mean(InvNuByV,2));
    InvNuByVErr = squeeze(std(InvNuByV,0,2)) / sqrt(Nensemble);
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    k = k + 1;
    variable_values = InvNuByVMean;
    errors = InvNuByVErr;
    % Now check if all variables are within error bounds of each other
    equal = true; % Assume they are equal to start
    for i = 1:length(variable_values)
        for j = i+1:length(variable_values)
            if abs(variable_values(i) - variable_values(j)) > (errors(i) + errors(j))
                equal = false; % If any are not within error, they are not equal
                break; % No need to check further if one inequality is found
            end
        end
        if ~equal
            break; % Exit the outer loop as well
        end
    end
    
    % If all variables are equal within the error, add the temperature to the list
    if equal
        equal_within_error = [equal_within_error, Ti];
        InvNuRange = [InvNuRange, mean(InvNuByVMean)];
    end
end
equal_within_error_unmodified = equal_within_error;
equal_within_error = equal_within_error - Tindices(1);
% Display the temperature range where variables are equal within errors
if ~isempty(equal_within_error)
    fprintf('Temperature range = \n %f to %f.\n', ...
            TRange(min(equal_within_error)), TRange(max(equal_within_error)));
    fprintf('respective temperature indices = \n %i to %i\n', ...
            (min(equal_within_error_unmodified)), (max(equal_within_error_unmodified)));
    fprintf('Index range width = \n %i \n', ...
            max(equal_within_error) - min(equal_within_error));
    fprintf('InvNu range = \n %f to %f.\n', ...
            InvNuMean(min(equal_within_error)), InvNuMean(max(equal_within_error)))
    fprintf('nu range = \n %f to %f.\n', ...
        1/InvNuMean(max(equal_within_error)), 1/InvNuMean(min(equal_within_error)))
else
    disp('No temperature range found where all variables are equal');
end
sort(equal_within_error);
IndexTc = ceil(mean(equal_within_error));
Tc = TRange(IndexTc);
deltaT = TRange(max(equal_within_error)) - T(min(equal_within_error));
InvNu = InvNuMean(IndexTc);
deltaInvNu = InvNuMean(max(equal_within_error)) - InvNuMean(min(equal_within_error));
nu = 1/InvNuMean(IndexTc);
deltaNu = 1/InvNuMean(min(equal_within_error)) - 1/InvNuMean(max(equal_within_error));
disp(['Tc = ' fmtMeanUnc(Tc, deltaT/2)]);
disp(['InvNu = ' fmtMeanUnc(InvNu, deltaInvNu/2)]);
disp(['nu = ' fmtMeanUnc(nu, deltaNu/2), newline]);

disp(['Tc index = ' num2str(IndexTc)]);
disp(['Tc (no error) = ' num2str(Tc)]);
disp(['InvNu (no error) = ' num2str(InvNu)]);
disp(['nu (no error) = ' num2str(nu)]);

fprintf('\n')

%%--------------------------------------------------------------------------------
%% Equal Nu By V Temperatue Range (bootstrap)
% Tindices = 1:length(T);
Tindices = 30:110;
TRange =  T(Tindices);
SizeRange = 1:21;

Nboot = 500;

% Initialize array to keep track of temperature points where conditions are met
equal_within_error = [];

InvNuRange = [];
% InvNuByV = zeros(6,Nensemble);
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);
% InvNuByVErr = zeros(6,1);
% InvNuByVMean = zeros(6,1);

k = 1;
for Ti = Tindices
    Nensemble = size(V, 2);
    InvNuByV = zeros(6,Nboot);
    for iboot = 1:Nboot
        ensembles = randi(Nensemble, Nensemble, 1);
        ensembles = ensembles';
        VBootTot = zeros(size(squeeze(V(Ti,1,SizeRange,:))));
        for Ensi = ensembles
            VBootTot = VBootTot + squeeze(V(Ti,Ensi,SizeRange,:));
        end
        VBootMean = squeeze(VBootTot) / Nensemble;
        for i = 1:6
            [p,S] = polyfit(log(Llist(SizeRange)), VBootMean(SizeRange,i),1);
            InvNuByV(i, iboot) = p(1);
            % InvNuByV(i,Ensi) = mdl.Coefficients.Estimate(2);
        end
    end
    InvNuByVMean = squeeze(mean(InvNuByV,2));
    InvNuByVErr = squeeze(std(InvNuByV,0,2));
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    k = k + 1;
    variable_values = InvNuByVMean;
    errors = InvNuByVErr;




    % Now check if all variables are within error bounds of each other
    equal = true; % Assume they are equal to start
    for i = 1:length(variable_values)
        for j = i+1:length(variable_values)
            if abs(variable_values(i) - variable_values(j)) > (errors(i) + errors(j))
                equal = false; % If any are not within error, they are not equal
                break; % No need to check further if one inequality is found
            end
        end
        if ~equal
            break; % Exit the outer loop as well
        end
    end
    
    % If all variables are equal within the error, add the temperature to the list
    if equal
        equal_within_error = [equal_within_error, Ti];
        InvNuRange = [InvNuRange, mean(InvNuByVMean)];
    end
end
equal_within_error_unmodified = equal_within_error;
equal_within_error = equal_within_error - Tindices(1);
% Display the temperature range where variables are equal within errors
if ~isempty(equal_within_error)
    fprintf('Temperature range = \n %f to %f.\n', ...
            TRange(min(equal_within_error)), TRange(max(equal_within_error)));
    fprintf('respective temperature indices = \n %i to %i\n', ...
            (min(equal_within_error_unmodified)), (max(equal_within_error_unmodified)));
    fprintf('Index range width = \n %i \n', ...
            max(equal_within_error) - min(equal_within_error));
    fprintf('InvNu range = \n %f to %f.\n', ...
            InvNuMean(min(equal_within_error)), InvNuMean(max(equal_within_error)))
    fprintf('nu range = \n %f to %f.\n', ...
        1/InvNuMean(max(equal_within_error)), 1/InvNuMean(min(equal_within_error)))
else
    disp('No temperature range found where all variables are equal');
end
sort(equal_within_error);
IndexTc = ceil(mean(equal_within_error));
Tc = TRange(IndexTc);
deltaT = TRange(max(equal_within_error)) - T(min(equal_within_error));
InvNu = InvNuMean(IndexTc);
deltaInvNu = InvNuMean(max(equal_within_error)) - InvNuMean(min(equal_within_error));
nu = 1/InvNuMean(IndexTc);
deltaNu = 1/InvNuMean(min(equal_within_error)) - 1/InvNuMean(max(equal_within_error));
disp(['Tc = ' fmtMeanUnc(Tc, deltaT/2)]);
disp(['InvNu = ' fmtMeanUnc(InvNu, deltaInvNu/2)]);
disp(['nu = ' fmtMeanUnc(nu, deltaNu/2), newline]);

disp(['Tc index = ' num2str(IndexTc + Tindices(1))]);
disp(['Tc (no error) = ' num2str(Tc)]);
disp(['InvNu (no error) = ' num2str(InvNu)]);
disp(['nu (no error) = ' num2str(nu)]);

fprintf('\n')
%%---------------------------------------------------------------------------------
%% nu by V Ensemble plot
% clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
indexT = 1;
maxidxL = 0;
% Tindices = [26, 36, 41:2:48, 53, 63];
Tindices = 32:8:64;
InvNuByV = zeros(6,Nensemble);
InvNuByVMean = zeros(6,1);
InvNuByVErr = zeros(6,1);
randColor = rand(2,1);
SizeRange= 1:21;
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);

fig = figure();
ax = axes(fig);
cmap = linspecer(length(Tindices)+1);
% cmap = turbo(length(SizeRange)+1);
ax.ColorOrder = cmap;
k = 1;
mi = 1;
for Ti = Tindices
    for Ensi = 1:Nensemble
        for i = 1:6
            [p,S] = polyfit(log(Llist(SizeRange)), V(Ti,Ensi,SizeRange,i),1);
            % mdl = fitlm(log(Llist(SizeRange)), squeeze(V(Ti,Ensi,SizeRange,i)));
    %         sss= table2array(mdl.Coefficients);
    %         nuByVError(i) =  mdl.Coefficients.SE(2);
    %         nuByVError(i) = S.normr;
            InvNuByV(i,Ensi) = p(1);  % why real part???
            % InvNuByV(i,Ensi) = mdl.Coefficients.Estimate(2);
        end
    end
    InvNuByVMean = squeeze(mean(InvNuByV,2));
    InvNuByVErr = squeeze(std(InvNuByV,0,2)) / sqrt(Nensemble);
    h = errorbar(1:6,InvNuByVMean,InvNuByVErr, ...
     'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 6, ...
     'LineStyle', 'none');
    % h = plot(1:6,InvNuByVMean, ...
    %  'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 6, ...
    %  'LineStyle', 'none');
    if mi < 10
        h.MarkerFaceColor = h.Color;
    end
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    X = 6.6;
    Y = InvNuByVMean;
    text(X(end), Y(end), num2str(T(Ti)));
    k = k + 1;
    mi = mi + 1;
    hold on
end
% ylim([0.86 1.06])
yl = ylim;
text(X(end), (yl(2)*(1.004)), '\quad $T$');
legendStrings = "T = " + string(T(Tindices));
% legend(legendStrings,'Interpreter','latex','Location','northeastoutside');
xlim([0.5 6.5]);
xlabel '$V_i$';
ylabel '$1/ \nu$';
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) ;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'nuByV, Error by Bins.png'])

%%---------------------------------------------------------------------------------
%% nu by V Ensemble plot no Bins
% clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
indexT = 1;
maxidxL = 0;
% Tindices = [26, 36, 41:2:48, 53, 63];
Tindices = 16:8:72;
InvNuByV = zeros(6);
SizeRange= 1:N;
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);
InvNuByVErr = zeros(6,1);
InvNuByVMean = zeros(6,1);
randColor = rand(2,1);

fig = figure();
ax = axes(fig);
cmap = linspecer(length(Tindices)+1);
% cmap = turbo(length(SizeRange)+1);
ax.ColorOrder = cmap;
k = 1;
mi = 1;

for Ti = Tindices
    for i = 1:6
        [p,S] = polyfit(log(Llist(SizeRange)), VMean(Ti,SizeRange,i),1);
        InvNuByV(i) = p(1);
        mdl = fitlm(log(Llist(SizeRange)), squeeze(V(Ti,Ensi,SizeRange,i)));
        InvNuByVMean(i) = mdl.Coefficients.Estimate(2);
        InvNuByVErr(i) = mdl.Coefficients.SE(2);
        
    end
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    h = errorbar(1:6,InvNuByVMean, InvNuByVErr, ...
     'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 6, ...
     'LineStyle', 'none');
    % h = plot(1:6,InvNuByVMean, ...
    %  'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 6, ...
    %  'LineStyle', 'none');
    if mi < 10
        h.MarkerFaceColor = h.Color;
    end
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    X = 6.6;
    Y = InvNuByVMean;
    text(X(end), Y(end), num2str(T(Ti)));
    k = k + 1;
    mi = mi + 1;
    hold on
end
% ylim([0.86 1.06])
yl = ylim;
text(X(end), (yl(2)*(1.004)), '\quad $T$');
legendStrings = "T = " + string(T(Tindices));
% legend(legendStrings,'Interpreter','latex','Location','northeastoutside');
xlim([0.5 6.5]);
xlabel '$V_i$';
ylabel '$1/ \nu$';

ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) ;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'nuByV, Error by SE of regression.png'])

%%---------------------------------------------------------------------------------
%% nu by V Ensemble plot (bootstrap)
% clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
% invsize = 1 ./ Nlist;
indexT = 1;
maxidxL = 0;
% Tindices = [26, 36, 41:2:48, 53, 63];
Tindices = 50:4:80;
InvNuByV = zeros(6,Nensemble);
InvNuByVMean = zeros(6,1);
InvNuByVErr = zeros(6,1);
randColor = rand(2,1);
SizeRange= 1:21;
InvNuStd = zeros(length(Tindices),1);
InvNuMean = zeros(length(Tindices),1);

fig = figure();
ax = axes(fig);
cmap = linspecer(length(Tindices)+1);
% cmap = turbo(length(SizeRange)+1);
ax.ColorOrder = cmap;
k = 1;
mi = 1;
for Ti = Tindices
    % for Ensi = 1:Nensemble
    %     for i = 1:6
    %         [p,S] = polyfit(log(Llist(SizeRange)), V(Ti,Ensi,SizeRange,i),1);
    %         % mdl = fitlm(log(Llist(SizeRange)), squeeze(V(Ti,Ensi,SizeRange,i)));
    % %         sss= table2array(mdl.Coefficients);
    % %         nuByVError(i) =  mdl.Coefficients.SE(2);
    % %         nuByVError(i) = S.normr;
    %         InvNuByV(i,Ensi) = p(1);  % why real part???
    %         % InvNuByV(i,Ensi) = mdl.Coefficients.Estimate(2);
    %     end
    % end
    % InvNuByVMean = squeeze(mean(InvNuByV,2));
    % InvNuByVErr = squeeze(std(InvNuByV,0,2)) / sqrt(Nensemble);

    Nensemble = size(V, 2);
    Nboot = 100;
    InvNuByV = zeros(6,Nboot);
    for k = 1:Nboot
        ensembles = randi(Nensemble, Nensemble, 1);
        ensembles = ensembles';
        VBootTot = zeros(size(squeeze(V(Ti,1,SizeRange,:))));
        for Ensi = ensembles
            VBootTot = VBootTot + squeeze(V(Ti,Ensi,SizeRange,:));
        end
        VBootMean = squeeze(VBootTot) / Nensemble;
        for i = 1:6
            [p,S] = polyfit(log(Llist(SizeRange)), VBootMean(SizeRange,i),1);
            InvNuByV(i, k) = p(1);
            % InvNuByV(i,Ensi) = mdl.Coefficients.Estimate(2);
        end
    end
    InvNuByVMean = squeeze(mean(InvNuByV,2));
    InvNuByVErr = squeeze(std(InvNuByV,0,2));

    h = errorbar(1:6,InvNuByVMean,InvNuByVErr, ...
     'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 6, ...
     'LineStyle', 'none');
    % h = plot(1:6,InvNuByVMean, ...
    %  'Color',cmap(mi,:), 'Marker', markers{mi},'MarkerSize', 6, ...
    %  'LineStyle', 'none');
    if mi < 10
        h.MarkerFaceColor = h.Color;
    end
    InvNuStd(k) = std(InvNuByVMean);
    InvNuMean(k) = mean(InvNuByVMean);
    X = 6.6;
    Y = InvNuByVMean;
    text(X(end), Y(end), num2str(T(Ti)));
    k = k + 1;
    mi = mi + 1;
    hold on
end
% ylim([0.86 1.06])
yl = ylim;
text(X(end), (yl(2)*(1.004)), '\quad $T$');
legendStrings = "T = " + string(T(Tindices));
% legend(legendStrings,'Interpreter','latex','Location','northeastoutside');
xlim([0.5 6.5]);
xlabel '$V_i$';
ylabel '$1/ \nu$';
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) ;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'nuByV, Error by Bins.png'])
%%---------------------------------------------------------------------------------------
%% 
N = 12544;
T0 = data(1);
T = 0.6805;
beta0 = 1/T0;
beta = 1/T;
mDb = -(beta - beta0);
dT = 0.0001;

onesVector = ones(size(Etime(:,1)));
Etime = Etime - min(Etime,[],'all');
tic();

[dlog,denum] = logsumexp(Etime(:,1),0,onesVector);
[numlog]  = logsumexp(Etime(:,1),0,Etime(:,1));
E1L = exp(numlog - dlog);
[numlog]  = logsumexp(Etime(:,1),0,Etime(:,1).^2);
E2L = exp(numlog - dlog);
C = beta^2 * (E2L-E1L^2);
ValidDeltaT = T0/sqrt(C);
TRange =[T0-1*ValidDeltaT, T0+1*ValidDeltaT]


T = TRange(1):dT:TRange(2);
% T = 0.65:dT:0.73;
E1L = zeros(size(T));
E2L = zeros(size(T));
C = zeros(size(T));

for i = 1:length(T)
    beta = 1/T(i);
    mDb = -(beta - beta0);
    dlog = logsumexp(Etime(:,1),mDb,onesVector);
    E1L(i)  = exp(logsumexp(Etime(:,1),mDb,mtime(:,1))    - dlog);
    E2L(i) = exp(logsumexp(Etime(:,1),mDb,mtime(:,1).^2) - dlog);    
end
Cpp = T.^(-1) .* (E2L-E1L.^2) ./ N;

plot(T,Cpp, 'Marker', '.','MarkerSize', 7)

%% binning
T0 = data(1);
T = 0.67;
beta0 = 1/T0;
beta = 1/T;
mDb = -(beta - beta0);

[Nbin,EbinEdge] = histcounts (Etime(:,1),'NumBins',1000);
Ebin = (EbinEdge(1:end-1) + EbinEdge(2:end))/2;
WE = Nbin.*exp(mDb*Ebin);

plot(Ebin, WE)


toc()

%%---------------------------------------------------------------------------------
%% FSS Binder

Tc = 0.844;
reducedTemp = data(:,1,:) / Tc -1;
L = sqrt(listN);
for i=1:N
    Xtilde(:,i) = data(:,7,i);  %Binder
    x(:,i) = L(i)^(1) * reducedTemp(:,:,i);
end

clf('reset')
for i = 1:N
    plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 12);
%     plot(x(:,i), data(:,7,i), '.','MarkerSize', 10);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
% xlim([-5 5])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));

% ylabel '$m L^{\beta / \nu}$'
% ylabel '$\chi L^{-\gamma / \nu}$'
ylabel '$U_4$'
xlabel('$t L^{1/ \nu}$')
legendStrings = "L = " + string(listL);
legend(legendStrings, 'Location','southeast','Interpreter','latex');
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
cmap = turbo(N+1);
ax.ColorOrder = cmap;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'Fss binder.png'],'Resolution',dpi)


%%---------------------------------------------------------------------------------
%% FSS chi
dpi = 600;
fontSize = 14;
Tc = TCritical;
reducedTemp = (T / Tc) -1;
L = Llist;
x = [];
Xtilde = [];
sizeRange = 1:15;
for i=sizeRange
    Xtilde(:,i) = L(i)^(-7/4) .* XppMean(:, i);  
    x(:,i) = L(i)^(1) * reducedTemp;
end

clf('reset')
for i = sizeRange
    plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 12);
%     plot(x(:,i), data(:,7,i), '.','MarkerSize', 10);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
% xlim([-5 5])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',fontSize) 
% ylabel '$m L^{\beta / \nu}$'
ylabel '$\chi L^{-\gamma / \nu}$'
% ylabel '$U_4$'
xlabel('$t L^{1/ \nu}$')
legendStrings = "L = " + string(Llist);
% legend(legendStrings, 'Location','northwestoutside','Interpreter','latex');
% legend(legendStrings, 'Location','northwest','Interpreter','latex');
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
cmap = turbo(N+1);
ax.ColorOrder = cmap;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'Fss chi.png'],'Resolution',dpi)


%%---------------------------------------------------------------------------------
%% FSS m

dpi = 600;
fontSize = 14;
Tc = TCritical;
reducedTemp = (T / Tc) -1;
L = Llist;
x = [];
Xtilde = [];
sizeRange = 10:15;
for i=sizeRange
    Xtilde(:,i) = L(i)^(1/8) .* mppMean(:, i);  
    x(:,i) = L(i)^(1) * reducedTemp;
end

clf('reset')
for i = sizeRange
    plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 12);
%     plot(x(:,i), data(:,7,i), '.','MarkerSize', 10);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
% xlim([-5 5])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',fontSize) 
ylabel '$m L^{\beta / \nu}$'
% ylabel '$\chi L^{-\gamma / \nu}$'
% ylabel '$U_4$'
xlabel('$t L^{1/ \nu}$')
% legendStrings = "L = " + string(listL);
% legend(legendStrings, 'Location','southwest','Interpreter','latex');
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
cmap = turbo(N+1);
ax.ColorOrder = cmap;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'Fss m.png'],'Resolution',dpi)

%%---------------------------------------------------------------------------------
%% test
i = 10;
% y = XppMean(:,i);
% yerr = XppErr(:,i);
y = Binder4DerMean(:,i);
yerr = Binder4DerErr(:,i);

x = T;
[maxY, maxIndex] = max(y);
xAtMax = x(maxIndex);

% upperCrossings = find(y + yerr >= maxY);
upperCrossings = find(y + yerr >= maxY - yerr(maxIndex));

upperCrossings = extractSequence(upperCrossings, maxIndex);
xErrorLow = x(upperCrossings(1));
xErrorHigh = x(upperCrossings(end));

fprintf('Max y at x = %.4f (%.4f - %.4f)\n', xAtMax, xErrorLow, xErrorHigh);
xAtMaxErr = max(xErrorHigh - xAtMax, xAtMax - xErrorLow)
% x_err = xAtMaxWithError(T, Binder4DerMean(:,i), Binder4DerErr(:,i))
% x_err = xAtMaxWithError(T, DK3Mean(:,i), DK3Err(:,i))
% x_err = xAtMaxWithError(T, -DK3Mean(:,i), DK3Err(:,i))
% x_err = xAtMaxWithError(T, DK2Mean(:,i), DK2Err(:,i))
% x_err = xAtMaxWithError2(T, DK2(:, :, i))
% x_err = xAtMaxWithError3(T, DK2(:, :, i))

% x_err = xAtMaxWithError(T, -DK2Mean(:,i), DK2Err(:,i))
% x_err = xAtMaxWithError2(T, -DK2(:, :, i))
% x_err = xAtMaxWithError3(T, -DK2(:, :, i))

x_err = xAtMaxWithError(T, Binder4DerMean(:,i), Binder4DerErr(:,i))
x_err = xAtMaxWithError2(T, Binder4Der(:, :, i))
x_err = xAtMaxWithError3(T, Binder4Der(:, :, i))
%%---------------------------------------------------------------------------------
%% functions


function x_xErr = xAtMaxWithError(x, y, yerr)
    % [maxY, maxIndex] = max(y);
    [maxY, maxIndex] = findpeaks(y);
    maxY = maxY(1);
    maxIndex = maxIndex(1);
    xAtMax = x(maxIndex);
    
    upperCrossings = find(y + yerr >= maxY);
    % upperCrossings = find(y + yerr >= maxY - yerr(maxIndex));
    
    upperCrossings = extractSequence(upperCrossings, maxIndex);
    xErrorLow = x(upperCrossings(1));
    xErrorHigh = x(upperCrossings(end));

    xAtMaxErr = max(xErrorHigh - xAtMax, xAtMax - xErrorLow);
    x_xErr = [xAtMax, xAtMaxErr];
end

function x_xErr = xAtMaxWithError2(x, yEns)
    % [maxY, maxIndex] = max(y);
    xAtMax = [];
    Nensemble = size(yEns, 2);
    for  i = 1:Nensemble
        y = squeeze(yEns(:, i));
        [maxY, maxIndex] = findpeaks(y);
        maxY = maxY(1);
        maxIndex = maxIndex(1);
        xAtMax(i) = x(maxIndex);
    end
    x_xErr = [];
    x_xErr(1) = mean(xAtMax);
    x_xErr(2) = std(xAtMax) / sqrt(length(xAtMax));
end

function x_xErr = xAtMaxWithError3(x, yEns)
    % [maxY, maxIndex] = max(y);
    xAtMax = [];
    Nensemble = size(yEns, 2);
    Nboot = 100;
    for k = 1:Nboot
        ensembles = randi(Nensemble, Nensemble, 1);
        ensembles = ensembles';
        y = zeros(size(squeeze(yEns(:, 1))));
        for l = ensembles
            y = y + squeeze(yEns(:, l));
        end
        [maxY, maxIndex] = findpeaks(y);
        maxY = maxY(1);
        maxIndex = maxIndex(1);
        xAtMax(k) = x(maxIndex);
    end
    x_xErr = [];
    x_xErr(1) = mean(xAtMax);
    % x_xErr(2) = std(xAtMax) / sqrt(Nensemble - 1);
    x_xErr(2) = std(xAtMax);
end


function x_mean_Ens = xAtMaxInBootEns(x, yEns, BootEnsembles)
    y = zeros(size(squeeze(yEns(:, 1))));
    for l = BootEnsembles
        y = y + squeeze(yEns(:, l));
    end
    [maxY, maxIndex] = findpeaks(y);
    maxY = maxY(1);
    maxIndex = maxIndex(1);
    x_mean_Ens = x(maxIndex);
end

function BootAveraged = BootAverage(Param, BootEnsembles)
    % [maxY, maxIndex] = max(y);
    xAtMax = [];
    Nensemble = size(yEns, 2);
    Nboot = 100;
    for k = 1:Nboot
        ensembles = randi(Nensemble, Nensemble, 1);
        ensembles = ensembles';
        y = zeros(size(squeeze(yEns(:, 1))));
        for l = BootEnsembles
            y = y + squeeze(yEns(:, l));
        end
        [maxY, maxIndex] = findpeaks(y);
        maxY = maxY(1);
        maxIndex = maxIndex(1);
        xAtMax(k) = x(maxIndex);
    end
    x_xErr = [];
    x_xErr(1) = mean(xAtMax);
    % x_xErr(2) = std(xAtMax) / sqrt(Nensemble - 1);
    x_xErr(2) = std(xAtMax);

end

function sequence = extractSequence(arr, value)
    % Find the index of the given value
    idx = find(arr == value);
    
    if isempty(idx)
        sequence = [];
        return;
    end
    
    % Find the start and end of the sequence
    start = idx;
    while start > 1 && arr(start-1) == arr(start) - 1
        start = start - 1;
    end
    
    endIdx = idx;
    while endIdx < length(arr) && arr(endIdx+1) == arr(endIdx) + 1
        endIdx = endIdx + 1;
    end
    
    % Extract the sequence
    sequence = arr(start:endIdx);
end


