%---------------------------------------------%
% Plot averaged magnetization and polarization
%---------------------------------------------%

%Line style
pltcolor = [0.0001, 0.4470, 0.7410];
pltstl = '-';
NumberOfSweeps = 24;
figure
plotWeightedAveragesCPP2('FSdata\\FieldSweeps',...
    NumberOfSweeps, pltcolor, pltstl)  ;

%-------uncomment to plot only an individual sweep
%figure
%filename = 'FSdata\\FieldSweeps_0';
%oneplt(filename, pltcolor, pltstl)

return


%% save as pdf
x0 = 30;
y0 = 10;
width = 8.;%17.0;
height = 8.;%7.5;%15.0;

set(gcf,'units','centimeters','position',[x0,y0,width,height])

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
exportgraphics(h,'FSsum2.pdf')


function plotWeightedAveragesCPP2(filename, Nsw, pltcolor, pltstl)    
    filename_0 = [filename, sprintf('_%d',0)];
    A = readmatrix(filename_0);
    allResults = A'; clear A
    Hlist = squeeze(allResults(1,:)); clear allResults
    
    [Mav, MavSigma, Pav, PavSigma, L2av, L2avSigma]...
        = givemeWAverages(filename, Nsw, numel(Hlist));
    %figure plot
    %%Load file
MWav  = 0*Hlist; MWavSigma = 0*Hlist;
L2Wav = 0*Hlist; L2WavSigma= 0*Hlist;
PWav  = 0*Hlist; PWavSigma = 0*Hlist;
for iH = 1:numel(Hlist)
    [MWav(iH), MWavSigma(iH), L2Wav(iH), L2WavSigma(iH), PWav(iH), PWavSigma(iH) ]...
            = givemeWeightedAverages...
            (Mav(:,iH), MavSigma(:,iH), L2av(:,iH),...
            L2avSigma(:,iH), Pav(:,iH), PavSigma(:,iH));
end
%plot weigted averages


myfntsize_label = 9;
myfntsize_gca = 8;
marksize = 4;

%plot data first
subplot(3,1,1)


hold on
%h1 = plot(Hlist,MWav,'LineWidth', 2.0);
errorbar(Hlist,MWav,MWavSigma);
if ~(pltcolor == 0)
h1.Color = pltcolor;
end
h1.LineStyle = pltstl;
hold off

%xlabel('H (T)','Fontsize',myfntsize_label);
ylabel('M (\mu_B per f.u.)','Fontsize',myfntsize_label);

subplot(3,1,2)

hold on

%h3 = plot(Hlist,PWav,'LineWidth', 2.0);
errorbar(Hlist,PWav,PWavSigma);
if ~(pltcolor == 0)
h3.Color = pltcolor;
end
h3.LineStyle = pltstl;
hold off

%xlabel('H (T)','Fontsize',myfntsize_label);
ylabel('<\sigma_is_j>_{\perp}','Fontsize',myfntsize_label);


subplot(3,1,3)

hold on


gperp = 400;
gpar = -800;
PWav0 = gpar *(-0.9);%(-0.82);
errorbar(Hlist,PWav0 + gperp*PWav, gperp*PWavSigma);

%Pscale = 1200/0.8*(0.35) ; 
%PWav0 = 1200/0.8 * (-0.65)* (-0.82);
%errorbar(Hlist,PWav0 + Pscale*PWav,Pscale*PWavSigma);
if ~(pltcolor == 0)
h3.Color = pltcolor;
end
h3.LineStyle = pltstl;
hold off

xlabel('H (T)','Fontsize',myfntsize_label);
ylabel('P (\muC/m^2)','Fontsize',myfntsize_label);


end



function oneplt(filename, pltcolor, pltstl)

    myfntsize_label = 9;

marksize = 4;

A = readmatrix(filename);
allResults = A'; clear A;

Tlist = squeeze(allResults(1,:));

MFM = squeeze(allResults(2,:));
MFMsigma = squeeze(allResults(3,:));

Pz = squeeze(allResults(4,:));
Pzsigma = squeeze(allResults(5,:));


subplot(3,1,1)
hold on
%errorbar(Tlist,MFM,MFMsigma,'k')
h1 = plot(Tlist,MFM);
h1.Color = pltcolor;
h1.LineStyle = pltstl;
hold off

%xlabel('H (T)','Fontsize',myfntsize_label);
ylabel('M (\mu_B per f.u.)','Fontsize',myfntsize_label);



subplot(3,1,2)

hold on
%h3 = errorbar(Tlist,Pz,Pzsigma);
h3 = plot(Tlist,Pz,'--');
h3.Color = pltcolor;
h3.LineStyle = pltstl;
hold off

xlabel('H (T)','Fontsize',myfntsize_label);
ylabel('P (\muC/m^2)','Fontsize',myfntsize_label);

set([h1,h3],'LineWidth',2.0);

end



function [MWav, MWavSigma, L2Wav, L2WavSigma, PWav, PWavSigma ]...
        = givemeWeightedAverages...
        (MavAll, MavSigmaAll, L2avAll, L2avSigmaAll, PavAll, PavSigmaAll)
    %wegited averages: xbar = sum(x_i/sigma_i^2) / sum(1/sigma_i^2)
    %for each H point
    [Nsw,~] = size(MavAll);
    %MWav = zeros(1,NH); L2Wav=MWav; PWav=MWav; 
    %calc magnetization
    MWav = sum(MavAll./MavSigmaAll.^2)/sum(1./MavSigmaAll.^2);
    %size(MWav)
    %sum(MavAll./MavSigmaAll.^2,1)/sum(1./MavSigmaAll.^2,1);
    %size(MWav)
    MWav1 = repmat(MWav,Nsw,1);
    %size(MWav1)
    MWavSigma = sqrt(sum((MavAll - MWav1).^2,1)/(Nsw-1));
    %MWavSigma
    %size(MWavSigma)
    %calc L^2 of Neel vector
    L2avAll = sign(L2avAll).*L2avAll;
    L2Wav = sum(L2avAll./L2avSigmaAll.^2,1)/sum(1./L2avSigmaAll.^2,1);
    L2Wav1 = repmat(L2Wav,Nsw,1);
    L2WavSigma = sqrt(sum((L2avAll - L2Wav1).^2,1)/(Nsw-1));
    
    %calc polarization
    PWav = sum(PavAll./PavSigmaAll.^2,1)/sum(1./PavSigmaAll.^2,1);
    PWav1 = repmat(PWav,Nsw,1);
    PWavSigma = sqrt(sum((PavAll - PWav1).^2,1)/(Nsw-1));
end

function [MFM, MFMsigma, Pz, Pzsigma, AFMF, AFMFsigma]...
    = givemeWAverages(filename, Nsw, NH)
%[Mav, MavSigma, Pav, PavSigma, L2av, L2avSigma]...
%        = givemeWAverages(filename, Nsw);
    MFM = zeros(Nsw, NH); MFMsigma=MFM;
    Pz = MFM; Pzsigma = MFM; AFMF = MFM; AFMFsigma = MFM;
    for i=1:Nsw
        iN = i-1;
        filename_n = [filename, sprintf('_%d',iN)];
        A = readmatrix(filename_n);
        %A = readmatrix('CoolingResultsFS');
        allResults = A'; clear A;

        %Hlist = squeeze(allResults(1,:));
        MFM(i,:) = squeeze(allResults(2,:));
        MFMsigma(i,:) = squeeze(allResults(3,:));

        Pz(i,:) = squeeze(allResults(4,:));
        Pzsigma(i,:) = squeeze(allResults(5,:));

        AFMF(i,:) = squeeze(allResults(6,:));
        AFMFsigma(i,:) = squeeze(allResults(7,:));
    end
    %Acrate = squeeze(allResults(8,:));
end