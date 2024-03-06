%% Plot basis vectors for LF_EIT_63

load('/cubric/scratch/c1616132/ROMEG_R123456_CUSTOM5/ROM/Results/ROM/other/LF_EIT_63.mat')
model = '/home/c1616132/Documents/PhD/ROMEG/Models/Real/head_model.mat';

%% Make all the plots


%basis = [2 4 7 9 10 11 12 13 14 15];
basis = [11 4 2];
%basis = [8 7 6 5 4 3 2 1];
%basis = 1;
plotting = PlottingClass();

for i = basis
    plotting.plotHeadModel('model',model,'axis','y','cut',0.1025,'fill',ROM.V(1:end-ROM.L,i))
    caxis([-5e-3,5e-3])
    view(0,0)
    label = ['Basis', char(10), 'Vector', char(10), num2str(i)];
    %text(0,0,0.165,label,'FontWeight','bold','FontSize',15,'HorizontalAlignment','left')
    colorbar('FontSize',14)
end


%% Make into one plot

figlist=get(groot,'Children');

figure; tcl=tiledlayout(2,2,'TileSpacing',"tight");


for i = 1:numel(figlist)
    figure(figlist(i));
    ax=gca;
    c=colorbar('FontSize',20);
    ax.Parent=tcl;
    ax.Layout.Tile=i+1;

    if (i==1)
        c.Parent = tcl;
        c.Layout.Tile = 'East';
        c.Ticks=[-5e-3,0,5e-3];
        c.TickLabels = {'Min','','Max'};
    end
end