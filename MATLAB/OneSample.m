%% load data
clear
clc
addpath(pwd);
addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);
date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_Sample'];
[status, msg, msgID] = mkdir(save_address)

dpi = 600;
%% set subfigs
fig = figure;
fig.Position  = [0 0 500 200];
% tiledlayout(4,2)
axs = tight_subplot(1,2,[.05 .15],[.25 .07],[.13 .03])
%%

SystemSizes = [32,64, 96]
% MethodNames=["16","32","64", "96"]
legendStrings=[];
% fig = figure;
% fig.Position  = [0 0 250 180];
colordiff=2;
cmap = linspecer(length(SystemSizes)+colordiff);
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
lineStyles = {'-', '--', '-.'};
mi=0;
axes(axs(2))
for SystemSize = SystemSizes
mi = mi+1;
        
    path = ['E:/visualize/inputSample/' num2str(SystemSize) '/']
    
    
    data = readmatrix([path 'MeanX.csv']);
    error = readmatrix([path 'errorX.csv']);
    % data = dlmread([path 'MeanX.txt']);
    % error = dlmread([path 'errorX.txt']);
    % orderParam = data(:,1);
    % myxLabel = '$B\mu /\Lambda $';
    % myxLabel = '$k_{\mathrm{B}} T/\Lambda$';
    
    % error = error / sqrt(10); % for results before 2022.08.25
    
    
    % locations = readmatrix([path, '0/locations.txt']);
    % DistanceVec = readmatrix([path, 'DistanceVec.txt']);
    % Plabel = readmatrix([path, 'Plabel.txt']);
    % Mean_GVec = readmatrix([path, 'Mean_GVec.txt']);
    % Mean_GConVec = readmatrix([path, 'Mean_GConVec.txt']);
    % error_GVec = readmatrix([path, 'error_GVec.txt']);
    % error_GConVec = readmatrix([path, 'error_GConVec.txt']);
    % Xtime = readmatrix([path, '0/Xtime.txt']);
    % XtimeDiff = readmatrix([path, '0/XtimeDiff.txt']);
    % XtimeNormalized = readmatrix([path, '0/XtimeNormalized.txt']);
    % Etime = readmatrix([path, '9/EtimeSeries.csv']);
    % mtime = readmatrix([path, '9/mtimeSeries.csv']);
    % Etime = h5read([path '0/EnsembleResults.h5'],['/timeSeries/EtimeSeries']);
    % mtime = h5read([path '0/EnsembleResults.h5'],['/timeSeries/mtimeSeries']);
    Etime = h5read([path '0/EnsembleResults.h5'],['/TimeSeries/Energy']);
    mtime = h5read([path '0/EnsembleResults.h5'],['/TimeSeries/Magnetization']);
    % Etime = h5read([path '0/EnsembleResults.h5'],['/TimeSeries/Energy2step']);
    % mtime = h5read([path '0/EnsembleResults.h5'],['/TimeSeries/Magnetization2step']);
    % mtime = readmatrix([path, '2/mtimeSeries.csv']);
    % mtimeNormalized = readmatrix([path, '0/mtimeNormalized.txt']);
    N = SystemSize^2;
    Nensemble = 8;
    XtimeNensemble=[];
    tauMean = 0;
    X = 0;
    C = 0;
    tauVec=[];
    for i=1:Nensemble
        disp(['Ensemble ', num2str(i)])
        EnsembleIndex = num2str(i-1);
        data_path = [path, EnsembleIndex , '/'];
        % Etime = readmatrix([data_path 'EtimeSeries.csv']);
        % mtime = readmatrix([data_path 'mtimeSeries.csv']);
        Etime = h5read([path EnsembleIndex '/EnsembleResults.h5'],['/TimeSeries/Energy']);
        mtime = h5read([path EnsembleIndex '/EnsembleResults.h5'],['/TimeSeries/Magnetization']);
        Etime = double(Etime);
        mtime = double(mtime);
        mtimeNormalized = mtime - mean(mtime);  % mean(mtime,2) is column vector containing the mean of each row
        mw = fft(mtimeNormalized);
        Xw = abs(mw).^2;
        Xtime = ifft(Xw);
        XtimeNensemble(:,i) = Xtime;
        tau = 1;
        for jj = 1:30
            tau = sum(Xtime(1:floor(6*tau)))/Xtime(1);
        %     disp(['i = ',num2str(i), ', tau = ', num2str(tau)])
        end
        % disp(['tau = ', num2str(tau)])
        tau=tau +0.5;
        mtime2 = mtime.^2;
        Etime2 = Etime.^2;
        X = X + (mean(mtime2) - mean(mtime)^2)/(data(1,1)   *N);
        C = C + (mean(Etime2) - mean(Etime)^2)/((data(1,1)^2)*N);
        tauVec(i) = tau;
    end
    X = X /Nensemble
    C = C/ Nensemble
    tauVec = tauVec *11;
    tauMean=mean(tauVec)
    tauSTE=std(tauVec)/ (sqrt(Nensemble)) 
    XtimeMean = mean(XtimeNensemble,2);
    XtimeErr = std(XtimeNensemble, 0,2) / (sqrt(Nensemble));
    XtimeNormMean = XtimeMean ./ XtimeMean(1);
    XtimeNormErr = XtimeErr ./ XtimeMean(1);



    % Indices to keep (where you do not want NaN) for error
    keep_indices = [10, 50, 100, 200, 300]; % Example indices
    
    keep_indices = keep_indices(1:mi+1);


    % Create a logical mask of the same size as the data
    mask = false(size(XtimeNormErr));
    
    % Set the mask to true at the indices to keep
    mask(keep_indices) = true;
    
    % Assign NaN to all elements not in the keep_indices
    XtimeNormErr(~mask) = NaN;




    % negativeIndex = find(XtimeNormMean < 0, 1);
    % if ~isempty(negativeIndex)
    %     % Set all elements after the negative index to zero
    %     XtimeNormMean(negativeIndex:end) = 1e-30;
    % end
    % 
    % %%---------------------------------------------------------------------------------
    % %% X(t)/X(0) fig
    % clf('reset')
    Step = 1:1:400;
    Time = Step * 11;

    h = errorbar(Time, XtimeNormMean(Step), XtimeNormErr(Step) , 'LineStyle', lineStyles{mi}, ...
            'Color',cmap(mi+colordiff,:), 'LineWidth', 1.2);

    % h = plot(Time, XtimeNormMean(Step), 'LineStyle', lineStyles{mi}, ...
    %         'Color',cmap(mi+colordiff,:), 'LineWidth', 2);
    h.MarkerFaceColor = h.Color;
    hold on

    % errorbar(Time(selected_err_idx), XtimeNormMean(selected_err_idx), ...
    %          XtimeNormErr(selected_err_idx), 'k', 'none' , 'LineWidth', 0.6);

    legendStrings = [legendStrings, string(['$L =' sprintf('%i',SystemSize) , '$' ])];
    
    % end

end
myFontSize =13;
set(gca,'fontsize',myFontSize) 
ylabel ('$\phi_M(t) / \phi_M(0)$', 'FontSize', myFontSize)

xlabel('$t \; [\mathrm{MCS}]$', 'FontSize', myFontSize)
xlim([0 5000])
ylim([1e-2 1])

set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
%     ylim([1e-1 1])
% ylim([0 1])
%     legendStrings = "$T = $ " + string(data(TempNumArr,1));
% legendStrings = "$T = $ " + string(data(1));
leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', myFontSize-1);
set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'YScale', 'log') 
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% tau = sum(XtimeNormMean(1:50)) * 11
exportgraphics(gcf,[save_address '/',  'X(t)divX(0)' '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  'X(t)divX(0)' '.eps'])

%%---------------------------------------------------------------------------------
 %% one ensemble investigation 
Etime = double(Etime);
mtime = double(mtime);
mtimeNormalized = mtime - mean(mtime);  % mean(mtime,2) is column vector containing the mean of each row
mw = fft(mtimeNormalized);
Xw = abs(mw).^2;
Xtime = ifft(Xw);
tau = 1;
for i = 1:30
    tau = sum(Xtime(1:floor(6*tau)))/Xtime(1);
%     disp(['i = ',num2str(i), ', tau = ', num2str(tau)])
end
% disp(['tau = ', num2str(tau)])
mtime2 = mtime.^2;
Etime2 = Etime.^2;
X = (mean(mtime2) - mean(mtime)^2)/(data(1,1)   *N);
C = (mean(Etime2) - mean(Etime)^2)/((data(1,1)^2)*N);
%% fft 
mtimeNormalized = mtime - mean(mtime);  % mean(mtime,2) is column vector containing the mean of each row
mw = fft(mtimeNormalized);
Xw = abs(mw).^2;
Xtime = ifft(Xw);
% Xtime = Xtime1';

%% visualize statistics
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

%%---------------------------------------------------------------------------------
%% X(t)
clf('reset')
Step = 1:20;
% TempNumArr = 2:4:size(data,1);

i = 1;
% for TempNum = TempNumArr
    plot(Step, Xtime(1:length(Step)), 'Marker', '.','MarkerSize', 12);
%     plot(Step, Xtime(TempNum, 1:length(Step)), 'Marker', '.','MarkerSize', 12);
%     plot(Step, Xtime(1:length(Step)), 'Marker', '.','MarkerSize', 12);
    hold on
% end
set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
%     ylim([1e-3 1])
set(gca,'fontsize',12) 
ylabel '$\chi(t)$'
xlabel('$t$');
legendStrings = "$T = $ " + string(data(TempNumArr,1));
legend(legendStrings, 'Location','northeast','Interpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on');
%     set(gca, 'YScale', 'log') 
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'X(t).png'])

%%---------------------------------------------------------------------------------
%% D(t)
clf('reset')
Step = 1:230;
TempNumArr = 2:4:size(data,1);
i = 1;
for TempNum = TempNumArr
    plot(Step, XtimeDiff(TempNum, 1:length(Step)), 'Marker', '.','MarkerSize', 12);
%     plot(Step, Xtime(1:length(Step)), 'Marker', '.','MarkerSize', 12);
    hold on
end
set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
%     ylim([1e-3 1])
set(gca,'fontsize',12) 
ylabel '$D(t)$'
xlabel('$t$');
legendStrings = "$T = $ " + string(data(TempNumArr,1));
legend(legendStrings, 'Location','northwest','Interpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on');
%     set(gca, 'YScale', 'log') 
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'XDiff(t).png'])

%%---------------------------------------------------------------------------------
%% D(t) movie
clf('reset')
Step = 1:1000;
TempNumArr = 1:1:size(data,1);
i = 1;
for TempNum = TempNumArr
    clf('reset')
    plot(Step, XtimeDiff(TempNum, 1:length(Step)), 'Marker', '.','MarkerSize', 12);
    set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
%     ylim([1e-3 1])
    set(gca,'fontsize',12) 
    ylabel '$D(t)$'
    xlabel('$t$');
    legendStrings = "$T = $ " + string(data(TempNum,1));
    legend(legendStrings, 'Location','southwest','Interpreter','latex');
    title(['$T = $ ', num2str(data(TempNum,1), '%4.3f')])
    set(gca,'XMinorTick','on','YMinorTick','on');
%     set(gca, 'YScale', 'log') 
    ax = gca;
    ax.TickLength(1) = 0.02;
    set(gca,'TickLabelInterpreter','latex');
    set(0,'defaulttextinterpreter','latex');

    F(i) = getframe(gcf);
    i = i + 1;
    saveas(gcf, [save_address '/', 'D(t)_', num2str(data(TempNum,1), '%1.3f') '.png'])
end

writerObj = VideoWriter([save_address, '/D(t)mov.avi'], 'Motion JPEG AVI');
writerObj.FrameRate = 1;
open(writerObj); % open the video writer
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);
close(gcf);

%%---------------------------------------------------------------------------------
%% X(t)/X(0) fig
clf('reset')
Step = 1:84000;
% TempNumArr = 2:4:size(data,1);
i = 41;
XtimeReduced = Xtime(1:length(Step))./Xtime(1);
% for TempNum = TempNumArr
% plot(Step, Xtime(TempNum, 1:length(Step))/Xtime(TempNum, 1), 'Marker', '.','MarkerSize', 12);    
plot(Step, XtimeReduced , 'Marker', '.','MarkerSize', 12);
    hold on
% end

    set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
%     ylim([1e-1 1])
    ylim([0 1])
    set(gca,'fontsize',12) 
    ylabel '$\chi(t) / \chi(0)$'
    xlabel('$t$');
%     legendStrings = "$T = $ " + string(data(TempNumArr,1));
    legendStrings = "$T = $ " + string(data(1));
    legend(legendStrings, 'Location','northeast','Interpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on');
%     set(gca, 'YScale', 'log') 
    ax = gca;
    ax.TickLength(1) = 0.02;
    set(gca,'TickLabelInterpreter','latex');
    set(0,'defaulttextinterpreter','latex');
    saveas(gcf, [save_address '/', 'X(t)divX(0).png'])
   tau = sum(Xtime(1:length(Step)))/Xtime(1)
%%---------------------------------------------------------------------------------
%% evaluating correlation time recursively
    tau = 1
    for i = 1:30
        tau = sum(Xtime(1:floor(6*tau)))/Xtime(1);
        disp(['i = ',num2str(i), ', tau = ', num2str(tau)])
    end
%%---------------------------------------------------------------------------------
%% X(t)/X(0) movie
clf('reset')
Step = 1:100;
TempNumArr = 1:1:size(data,1);
i = 1;
for TempNum = TempNumArr
    clf('reset')
    plot(Step, XtimeNormalized(TempNum, 1:length(Step)), 'Marker', '.','MarkerSize', 12);
    set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
    ylim([1e-3 1])
    set(gca,'fontsize',12) 
    ylabel '$\chi(t) / \chi(0)$'
    xlabel('$t$');
    legendStrings = "$T = $ " + string(data(TempNum,1));
    legend(legendStrings, 'Location','southwest','Interpreter','latex');
    title(['$T = $ ', num2str(data(TempNum,1), '%4.3f')])
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca, 'YScale', 'log') 
    ax = gca;
    ax.TickLength(1) = 0.02;
    set(gca,'TickLabelInterpreter','latex');
    set(0,'defaulttextinterpreter','latex');

    F(i) = getframe(gcf);
    i = i + 1;
%     saveas(gcf, [save_address '/', 'm(t)_', num2str(TempNum, '%1.3f') '.png'])
end

writerObj = VideoWriter([save_address, '/X(t)divX(0).mp4'], 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj); % open the video writer
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);
close(gcf);

%%---------------------------------------------------------------------------------
%% m(t)
clf('reset')
% mtimeNormalized = mtime - mean(mtime,2)
Step = 1:5000;
TempNumArr = 1:1:size(data,1);
i = 1;
for TempNum = TempNumArr
    clf('reset')
    % for TempNum = TempNumArr
        plot(Step, mtimeNormalized(TempNum, 1:length(Step)),'LineWidth',0.7);
    %     hold on
    % end
    set(gca,'TickLabelInterpreter','latex');
%     ylim([0, 1])
    yl = ylim;
    xl = xlim;
    yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
    xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
    set(gca,'fontsize',12) 
    ylabel '$m(t) -\bar{m}$'
    % xlabel(['$T = $ ', num2str(data(TempNum,1), '%4.3f')]);
    xlabel('$t$');
    legendStrings = "$T = $ " + string(data(TempNum,1));
    legend(legendStrings, 'Location','southeast','Interpreter','latex');
    title(['$T = $ ', num2str(data(TempNum,1), '%4.3f')])
    set(gca,'XMinorTick','on','YMinorTick','on');
    ax = gca;
    ax.TickLength(1) = 0.02;
    set(gca,'TickLabelInterpreter','latex');
    set(0,'defaulttextinterpreter','latex');
    F(i) = getframe(gcf);
    i = i + 1;
%     saveas(gcf, [save_address '/', 'm(t)_', num2str(TempNum, '%1.3f') '.png'])
end

writerObj = VideoWriter([save_address, '/m(t)bar_', num2str(length(Step), '%i'), '.mp4'], 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj); % open the video writer
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);
close(gcf);

%%---------------------------------------------------------------------------------
%% m(t) Normalized
clf('reset')
Step = 1:1000;
TempNumArr = 1:1:size(data,1);
for TempNum = TempNumArr
    clf('reset')
    % for TempNum = TempNumArr
        plot(Step, mtimeNormalized(TempNum, 1:length(Step)),'LineWidth',0.7);
    %     hold on
    % end
    set(gca,'TickLabelInterpreter','latex');
    ylim([-0.8, 0.8])
    yl = ylim;
    xl = xlim;
    yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
    xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
    set(gca,'fontsize',12) 
    ylabel '$m(t) - \langle m \rangle$'
    % xlabel(['$T = $ ', num2str(data(TempNum,1), '%4.3f')]);
    xlabel('$t$');
    legendStrings = "$T = $ " + string(data(TempNum,1));
    legend(legendStrings, 'Location','southeast','Interpreter','latex');
    title(['$T = $ ', num2str(data(TempNum,1), '%4.3f')])
    set(gca,'XMinorTick','on','YMinorTick','on');
    ax = gca;
    ax.TickLength(1) = 0.02;
    set(gca,'TickLabelInterpreter','latex');
    set(0,'defaulttextinterpreter','latex');
    saveas(gcf, [save_address '/', 'm(t)_', num2str(TempNum, '%1.3f') '.png'])
end

%%---------------------------------------------------------------------------------
%% G movie
clf('reset')
TempNumArr = 1:1:size(data,1);
mkdir([save_address '/G_images']);
xmax = max(DistanceVec);
i = 1;
for TempNum = TempNumArr
    clf('reset')
%     errorbar(DistanceVec, Mean_GConVec(TempNum, :), error_GVec(TempNum, :), 'Marker', '.','MarkerSize', 14);
%     plot(DistanceVec, Mean_GConVec(TempNum, :), 'Marker', '.','MarkerSize', 14);
    plot(DistanceVec, Mean_GConVec(TempNum, :), '.','MarkerSize', 14);
    set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
    ylim([1e-5 1])
    xlim([0 xmax])
    set(gca,'fontsize',12) 
    ylabel '$G$'
    xlabel('$r$');
    legendStrings = "$T = $ " + string(data(TempNum,1));
    legend(legendStrings, 'Location','southwest','Interpreter','latex');
    title(['$T = $ ', num2str(data(TempNum,1), '%4.3f')])
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca, 'YScale', 'log') 
%     set(gca, 'XScale', 'log') 
    ax = gca;
    ax.TickLength(1) = 0.02;
    set(gca,'TickLabelInterpreter','latex');
    set(0,'defaulttextinterpreter','latex');

    F(i) = getframe(gcf);
    i = i + 1;
    saveas(gcf, [save_address '/G_images/', 'G_', num2str(data(TempNum,1), '%1.3f') '.png'])
end

writerObj = VideoWriter([save_address, '/G_mov.avi'], 'Motion JPEG AVI');
writerObj.FrameRate = 1;
open(writerObj); % open the video writer
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);
close(gcf);

%%---------------------------------------------------------------------------------
%% Snapshot
clf('reset')
i = 1;
% sname = 'snapshot'; mkdir([save_address '/' sname '_images']); fid = fopen([path, '0/' sname '.txt']);
sname = 'SnapshotTransformed'; mkdir([save_address '/' sname '_images']); fid = fopen([path, '0/' sname '.txt']);

NumParticle = size(locations,2);
xmax = max(locations(1,:)) + 1;
xmin = min(locations(1,:)) - 1;
ymax = max(locations(2,:)) + 1;
ymin = min(locations(2,:)) - 1;
lozationZ = zeros(100,1);

gcf;
% set(gcf, 'Position',  [0, 0, 800, 600]);

while true
    clf('reset')
    Temperature = fscanf(fid, '%f', 1); if feof(fid); break; end
    Ux = fscanf(fid, '%f', NumParticle); if feof(fid); break; end
    Uy = fscanf(fid, '%f', NumParticle); if feof(fid); break; end
%     Uz = fscanf(fid, '%f', NumParticle); if feof(fid); break; end
    
    quiver(locations(1,:)',locations(2,:)',Ux, Uy, 0.4, 'LineWidth', 1);
    text(locations(1,:)',locations(2,:)',num2str(Plabel));
%     quiver3(locations(1,:)',locations(2,:)', lozationZ,Ux, Uy, Uz, 1, 'LineWidth', 1);
%     set(gca,'TickLabelInterpreter','latex');
%     yl = ylim;
%     xl = xlim;
%     yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
%     xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
    xlim([xmin xmax])
    ylim([ymin ymax])
    set(gca,'fontsize',12) 
    set(gca,'XColor', 'none','YColor','none')
%     ylabel '$G$'
%     xlabel('$r$');
%     legendStrings = "$T = $ " + string(Temperature);
%     legend(legendStrings, 'Location','southeast','Interpreter','latex');
%     title(['$T = $ ', num2str(Temperature, '%4.3f')])
    title(['T =  ', num2str(Temperature, '%4.3f')])
%     set(gca,'XMinorTick','on','YMinorTick','on');
    
%     ax = gca;
%     ax.TickLength(1) = 0.02;
%     set(gca,'TickLabelInterpreter','latex');
%     set(0,'defaulttextinterpreter','latex');

%     set(gca,'visible','off')

%     axis off

    F(i) = getframe(gcf);
    i = i + 1;
    saveas(gcf, [save_address '/' sname '_images/T=', ...
        num2str(Temperature, '%1.3f'), '_', num2str(i, '%d'), '.png'])
end
fclose(fid);

% % % % % writerObj = VideoWriter([save_address, '/snapshot_mov.avi'], 'Motion JPEG AVI');
% writerObj = VideoWriter([save_address, '/' sname '_mov.mp4'], 'MPEG-4');
% writerObj.FrameRate = 1;
% open(writerObj); % open the video writer
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% close(writerObj);

close(gcf);
disp('ended');

%%---------------------------------------------------------------------------------
%% Particle label

% plot(locations(1,:)',locations(2,:)', '.','MarkerSize', 14);
plot(locations(1,:)',locations(2,:)', '.','MarkerSize', 1);
text(locations(1,:)',locations(2,:)',num2str(Plabel));
set(gca,'XColor', 'none','YColor','none');
axis equal

%% Lattice Plot
plot(locations(1,:)',locations(2,:)', '.','MarkerSize', 14);
axis equal




