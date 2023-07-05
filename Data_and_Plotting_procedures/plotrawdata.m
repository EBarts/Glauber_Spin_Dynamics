%% plot raw data
figure
filename = 'PvsH_rawtab_T=25K_SweepForward';
plotfromfile(filename);
filename = 'PvsH_rawtab_T=25K_SweepBack';
plotfromfile(filename);
filename = 'PvsH_rawtab_T=25K_VirginCurve';
plotfromfile(filename);

figure
filename = 'PvsH_rawtab_T=35K_SweepForward';
plotfromfile(filename);
filename = 'PvsH_rawtab_T=35K_SweepBack';
plotfromfile(filename);
filename = 'PvsH_rawtab_T=35K_VirginCurve';
plotfromfile(filename);

figure
filename = 'MvsH_rawtab_T=35K_SweepForward';
plotfromfile(filename);
filename = 'MvsH_rawtab_T=35K_SweepBack';
plotfromfile(filename);
filename = 'MvsH_rawtab_T=35K_VirginCurve';
plotfromfile(filename);
function plotfromfile(filename)
    A = readmatrix(filename);
    allResults = A'; clear A
    Hlist = squeeze(allResults(1,:));
    XValues = squeeze(allResults(2,:));
    XstdValues =squeeze(allResults(3,:));
    hold on 
    errorbar(Hlist,XValues,XstdValues);
    hold off
end