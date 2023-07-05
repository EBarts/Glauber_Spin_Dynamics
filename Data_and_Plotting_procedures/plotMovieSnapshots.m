%

x = 1;%0.14;
MA = 4.21; MB = 4.51;

MS = MB - MA*(1-x);

%
ALLSPINS = importdata('SaveSnapshots\\ConfigFS.txt',',');
%
Hsize = 32; %Nstep = 1;%50;
%Hsize = 64; %Nstep = 1;%50;
Ny = 32; Nx = Ny; Nz = Ny;
SzAll = zeros(Hsize,Ny,Nx,Nz);
iH0 = 0;%32;%30; 
for iH=1:Hsize
    for iy = 1:Ny
        for ix = 1:Nx
            for iz = 1:Nz
                i = ix; j = iy - 1; k = iz - 1;
                SzAll(iH, iy, ix, iz )=...
                    ALLSPINS(iH+iH0,i + j*Nx + k*Nx*Ny);
            end
        end
    end
end

% MAKE EVOLUTION SNAPSHOTS

[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz/2);
%for 7 points
fig = figure('Position', [168    321    1584     584]);
set(gcf,'Renderer','opengl')
%set(gcf,'Renderer','painters');
set(gcf,'units','centimeters','position',[15,5, 11.0, 5.0]) %6.0

%Htoplt = [ 1, 8, 11, 14 , 17, 20, 24 ];
%Htoplt = [ 11, 8, 11, 14 , 17, 20, 24 ];
%Htoplt = 12:18; %24:30;18:24;%12:18;
Htoplt = [ 12, 13, 17 , 18, 20, 24, 28 ];


%
for iHplt = 1:numel(Htoplt)
    subplot(1,numel(Htoplt),iHplt)
    iH = Htoplt(iHplt);
    Sz = squeeze(SzAll(iH,:,:,:));
    visualizeSBLB(Sz,X,Y,Z)
    zlim([1 10])
    set(gca,'Zticklabel',[]) 
    set(gca,'Yticklabel',[]) 
    set(gca,'Xticklabel',[]) 
    xlabel('','Fontsize',20);
    ylabel('','Fontsize',20);
    zlabel('','Fontsize',20);
    set(gca,'Fontsize',10);
    set(gca,'LineWidth',1.0)
    set(gca,'box','on')
end

subplot(1,numel(Htoplt),1)
%set(gca,'Zticklabel',[2,4,6,8,10]) 


fname = 'snapAllcppFinal';
print(fig,fname,'-r600','-dpng'); %-r500
clear fname fig



function visualizeSBLB(Sz,X,Y,Z)
    [Ny,Nx,Nz] = size(Sz);
xslice = [];   
yslice = [];%0;%[];
zslice = 1:10;%1:12;%1:1:Nz;%1:Nz;%1:Nz/2;%1:Nz/2;%0;[1,3,4,5];

LBperUC= zeros(size(Z));
MBperUC= zeros(size(Z));

for iz = 1:Nz/2
   LBperUC(:,:,iz) = (Sz(:,:,2*iz-1) - Sz(:,:,2*iz))/2;
   MBperUC(:,:,iz)  = (Sz(:,:,2*iz-1) + Sz(:,:,2*iz))/2;
end
LBperUCplt = LBperUC + 0.2*MBperUC;


hh = slice(X,Y,Z,LBperUCplt,xslice,yslice,zslice);

%shading flat 
shading faceted
set(hh,'edgecolor','black','linewidth',0.1,'edgealpha',0.1); %'facealpha',0,
cmap = zeros(4,3); 

% 1->AFM=-1 ; 2->FM=+1; 3->FM=-1; 4->AFM=1
% white and black + blue and red color pallets
cmap(1,:) = [0, 154, 222]/255; % AFM down
cmap(2,:) = [247, 240, 232]/255;% FM up
cmap(3,:) = [0,0,0]/255;% FM down
cmap(4,:) = [255, 31, 91]/255; % AFM up

colormap(cmap);

xlabel('a','Fontsize',20);
ylabel('b','Fontsize',20);
xticks([0 60])
yticks([0 60])
zticks([2 4 6 8 10])

set(gca,'Fontsize',10);
view(-20,4)

end