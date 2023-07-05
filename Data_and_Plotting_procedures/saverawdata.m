%% T = 35K curves

%Sweep Forward
Nsw = 24;
%
filename = 'sim35KFinalFS20x20x20\\FSdata\\FieldSweeps';
[MWav, MWavSigma, PWav1, PWav1Sigma, Hlist]...
    = givemerawdata(filename, Nsw); 

xy = [Hlist(:), MWav(:), MWavSigma(:)];
writematrix(xy,'MvsH_rawtab_T=35K_SweepForward.txt','Delimiter','tab')

xy = [Hlist(:), PWav1(:), PWav1Sigma(:)];
writematrix(xy,'PvsH_rawtab_T=35K_SweepForward.txt','Delimiter','tab')

%Sweep Back
[MWav, MWavSigma, PWav1, PWav1Sigma, Hlist]...
    = givemerawdataBack(filename, Nsw); 

xy = [Hlist(:), MWav(:), MWavSigma(:)];
writematrix(xy,'MvsH_rawtab_T=35K_SweepBack.txt','Delimiter','tab')

xy = [Hlist(:), PWav1(:), PWav1Sigma(:)];
writematrix(xy,'PvsH_rawtab_T=35K_SweepBack.txt','Delimiter','tab')

%Virgin curves
%
filename = 'sim35KFinalFS20x20x20vc\\FSdata\\FieldSweeps';
[MWav, MWavSigma, PWav1, PWav1Sigma, Hlist]...
    = givemerawdata(filename, Nsw); 

xy = [Hlist(:), MWav(:), MWavSigma(:)];
writematrix(xy,'MvsH_rawtab_T=35K_VirginCurve.txt','Delimiter','tab')

xy = [Hlist(:), PWav1(:), PWav1Sigma(:)];
writematrix(xy,'PvsH_rawtab_T=35K_VirginCurve.txt','Delimiter','tab')

%% T = 25K curves

%Sweep Forward
Nsw = 24;
%
filename = 'sim25KFinalFS20x20x20\\FSdata\\FieldSweeps';
[MWav, MWavSigma, PWav1, PWav1Sigma, Hlist]...
    = givemerawdata(filename, Nsw); 

xy = [Hlist(:), MWav(:), MWavSigma(:)];
writematrix(xy,'MvsH_rawtab_T=25K_SweepForward.txt','Delimiter','tab')

xy = [Hlist(:), PWav1(:), PWav1Sigma(:)];
writematrix(xy,'PvsH_rawtab_T=25K_SweepForward.txt','Delimiter','tab')

%Sweep Back
[MWav, MWavSigma, PWav1, PWav1Sigma, Hlist]...
    = givemerawdataBack(filename, Nsw); 

xy = [Hlist(:), MWav(:), MWavSigma(:)];
writematrix(xy,'MvsH_rawtab_T=25K_SweepBack.txt','Delimiter','tab')

xy = [Hlist(:), PWav1(:), PWav1Sigma(:)];
writematrix(xy,'PvsH_rawtab_T=25K_SweepBack.txt','Delimiter','tab')

%Virgin curves
%
filename = 'sim25KFinalFS20x20x20vc\\FSdata\\FieldSweeps';
[MWav, MWavSigma, PWav1, PWav1Sigma, Hlist]...
    = givemerawdata(filename, Nsw); 

xy = [Hlist(:), MWav(:), MWavSigma(:)];
writematrix(xy,'MvsH_rawtab_T=25K_VirginCurve.txt','Delimiter','tab')

xy = [Hlist(:), PWav1(:), PWav1Sigma(:)];
writematrix(xy,'PvsH_rawtab_T=25K_VirginCurve.txt','Delimiter','tab')



return


function [MWav, MWavSigma, PWav1, PWav1Sigma, Hlist]...
    = givemerawdata(filename, Nsw)    
    filename_0 = [filename, sprintf('_%d',0)];
    %filename_0 = 'FieldSweeps_0';
    A = readmatrix(filename_0);
    allResults = A'; clear A
    Hlist = squeeze(allResults(1,:)); clear allResults
    
    %filename 'FieldSweeps'
    [Mav, MavSigma, Pav, PavSigma, L2av, L2avSigma]...
        = givemeWAverages(filename, Nsw, numel(Hlist));
    %size(Mav)
    %sadsad
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
%save weigted averages
PWav0 = 1200/0.8 * (-0.65)* (-0.8);
Pscale = 1200/0.8*(0.35);
PWav1 = PWav0 + Pscale*PWav;
PWav1Sigma = Pscale*PWavSigma;


end


function [MWavBack, MWavBackSigma, PWavBack1, PWavBack1Sigma, HlistBack]...
    = givemerawdataBack(filename, Nsw)    
    filename_0 = [filename, sprintf('_%d',0)];
    %filename_0 = 'FieldSweeps_0';
    A = readmatrix(filename_0);
    allResults = A'; clear A
    Hlist = squeeze(allResults(1,:)); clear allResults
    HlistBack = Hlist*0;
    for iH = 1:numel(Hlist)
        iHb = numel(Hlist) - iH + 1;
        HlistBack(iH) = -Hlist(iHb);
    end
    %filename 'FieldSweeps'
    [Mav, MavSigma, Pav, PavSigma, L2av, L2avSigma]...
        = givemeWAverages(filename, Nsw, numel(Hlist));
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
%save weigted averages
PWav0 = 1200/0.8 * (-0.65)* (-0.8);
Pscale = 1200/0.8*(0.35);
PWav1 = PWav0 + Pscale*PWav;
PWav1Sigma = Pscale*PWavSigma;

MWavBack  = 0*Hlist; MWavBackSigma = 0*Hlist;
PWavBack1  = 0*Hlist; PWavBack1Sigma = 0*Hlist;
    for iH = 1:numel(Hlist)
        iHb = numel(Hlist) - iH + 1;
        MWavBack(iH) = -MWav(iHb);
        MWavBackSigma(iH) = MWavSigma(iHb);
        PWavBack1(iH) = PWav1(iHb);
        PWavBack1Sigma(iH) = PWav1Sigma(iHb);
    end


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