%% Plot the main manuscript figure

samples=1:10;

c = [1,2,3,4,5,6]; c2 = [1,2,3];
snap = SnapShotClass('range',100,'sample_num',samples,'tissue',[3]);
snap = snap.readResults(c,c2);
tissues = {[1,2,3],[3],[1]};

%% Create figures

labels = {'3 Layers','Spongiform','Scalp'};
ypos = [6e-3,1.5e-2,2e-4];
xpos = 280;

for i=1:length(tissues)
    snap.tissue = tissues{i};
    snap = snap.processResults(c2);
    snap.plotSnapshots();
    xlim([0,300])
    lgd = legend('Location','southwest');
    text(xpos,ypos(i),labels{i},'FontWeight','bold','FontSize',14,'HorizontalAlignment','right')
end

%% Set all into one parent figure

figlist=get(groot,'Children');

figure; tcl=tiledlayout(3,1,'TileSpacing',"none");


for i = 1:numel(figlist)
    figure(figlist(i));
    ax=gca;
    lgd = legend();
    if (i==1)
        ylim(ax,[1e-4,4e-1])
        ax.YTick = [1e-4,1e-3,1e-2,1e-1];%ax.YTick(1:end-1);
    end
    if (i==2)
        ylim(ax,[9e-3,2])
        ax.YTick = [1e-2,1e-1,1e0];%ax.YTick(1:end-1);
    end
    if (i==3)
        ax.YTick = [1e-2,1e-1,1e0];
        ylim(ax,[2e-3,3])
    end
    if ~(i==numel(figlist))
        ax.XTickLabels = [];
        %ax.XTick = [];
    end
ax.Parent=tcl;
ax.Layout.Tile=i;

if (i==2)
    lgd.Parent = tcl;
end
end