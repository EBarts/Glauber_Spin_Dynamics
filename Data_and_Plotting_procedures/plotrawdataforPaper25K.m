%%plot raw data to reproduce T = 25 K data in the paper 
myfigure = figure;
set(gcf,'Renderer','painters')
left_fig_color = [0,0,0];
right_fig_color = left_fig_color;
set(myfigure,'defaultAxesColorOrder',[left_fig_color;right_fig_color]);

myfntsize_label = 8;
myfntsize_gca = 8;

yyaxis left

%axes
hold on
xline(0)
yline(0)
hold off

%Magnetization plots
ylabel('M (\mu_B/f.u.)','Fontsize', myfntsize_label);
filename = 'MvsH_rawtab_T=25K_VirginCurve';
plotMH(filename,"VC");

filename = 'MvsH_rawtab_T=25K_SweepBack';
plotMH(filename,"H_reversal_Back");

filename = 'MvsH_rawtab_T=25K_SweepForward';
plotMH(filename,"H_reversal_Forward");


yyaxis right
ylabel('P (\muC/m^2)','Fontsize', myfntsize_label);

%------- Uncomment the following lines to add STD shading plots ----------
% filename = 'PvsH_rawtab_T=35K_VirginCurve';
% plotPstdH(filename,"VC");
% 
% filename = 'PvsH_rawtab_T=35K_SweepBack';
% plotPstdH(filename,"H_reversal_Back");
% 
% filename = 'PvsH_rawtab_T=35K_SweepForward';
% plotPstdH(filename,"H_reversal_Forward");

%Polarization plots
filename = 'PvsH_rawtab_T=25K_VirginCurve';
plotPH(filename,"VC");

filename = 'PvsH_rawtab_T=25K_SweepBack';
plotPH(filename,"H_reversal_Back");

filename = 'PvsH_rawtab_T=25K_SweepForward';
plotPH(filename,"H_reversal_Forward");


ylim([100 1400]);
set(gca,  'YTick'       , [400, 800 ,1200]);

xticks([-10, -5, 0, 5, 10]);
xlabel('\mu_0H (T)','Fontsize', myfntsize_label);


set(gca,'Fontsize',myfntsize_gca);
set(gca,'LineWidth',0.9)
set(gca,'box','on')

yyaxis left

set(gca,'Fontsize',myfntsize_gca);
set(gca,'Layer','top')
set(gca,'LineWidth',0.7)
set(gca,'box','on')

text(0.05,0.85,'25 K','Units','normalized','FontSize',8)

%save plot as pdf
x0 = 30;
y0 = 10;
width = 7.;%17.0;
height = 4.0;%8.;%7.5;%15.0;

set(gcf,'units','centimeters','position',[x0,y0,width,height])

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
exportgraphics(h,'SweepsForPaperT=25K.pdf')

function plotMH(filename, indicator)
    [Hlist, Magnetization, ~] = read_the_file(filename);
    hold on
    plotHM = plot(Hlist, Magnetization);
    hold off

    plotHM.LineStyle = '--';
    plotHM.LineWidth = 1.5;
    myblue = [2 , 143, 163]/255;
    mygreen = [143 , 232, 10]/255;
    if indicator == "VC"
        plotHM.Color = 'k';

    elseif indicator == "H_reversal_Forward"
        plotHM.Color = mygreen; 
    elseif indicator == "H_reversal_Back"
        plotHM.Color = myblue;
    end
end

function plotPH(filename, indicator)
    [Hlist, Polarization, ~] = read_the_file(filename);
    hold on
    plotHM = plot(Hlist, Polarization);
    hold off
    plotHM.LineStyle = '-';
    plotHM.Marker = '.';
    plotHM.LineWidth = 0.5;
    plotHM.MarkerSize = 11;
    myblue = [2 , 143, 163]/255;
    mygreen = [143 , 232, 10]/255;
    if indicator == "VC"
        plotHM.Color = 'k';

    elseif indicator == "H_reversal_Forward"
        plotHM.Color = mygreen; 
    elseif indicator == "H_reversal_Back"
        plotHM.Color = myblue;
    end
end


function plotPstdH(filename, indicator)
    [Hlist, Polarization, Polarization_std] = read_the_file(filename);
    
    methsmooth = 'sgolay'; 
    Polarization0_Top = smooth(Polarization+Polarization_std/2,methsmooth);
    Polarization0_Bot = smooth(Polarization-Polarization_std/2,methsmooth);
    Polarization0_Top = Polarization0_Top';
    Polarization0_Bot = Polarization0_Bot';
    Hlist1 = Hlist;
    myblue = [2 , 143, 163]/255;
    mygreen = [143 , 232, 10]/255;
    if indicator == "VC"
        hold on
        Ppatch = patch([Hlist1 fliplr(Hlist1)],...
            [Polarization0_Bot fliplr(Polarization0_Top)],...
            [0, 0, 0],'EdgeAlpha',0.01);
        hold off

    elseif indicator == "H_reversal_Forward"
        hold on
        Ppatch = patch([Hlist1 fliplr(Hlist1)],...
            [Polarization0_Bot fliplr(Polarization0_Top)],...
            mygreen,'EdgeAlpha',0.01);
        hold off
    elseif indicator == "H_reversal_Back"
        hold on
        Ppatch = patch([Hlist1 fliplr(Hlist1)],...
            [Polarization0_Bot fliplr(Polarization0_Top)],...
            myblue,'EdgeAlpha',0.01);
        hold off
    end
    Ppatch.FaceAlpha = 1.0;

end

function [Hlist, XValues, XstdValues] = read_the_file(filename)
    A = readmatrix(filename);
    allResults = A'; clear A
    Hlist = squeeze(allResults(1,:));
    XValues = squeeze(allResults(2,:));
    XstdValues =squeeze(allResults(3,:));
    
end


