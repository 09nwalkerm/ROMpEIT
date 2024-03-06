%% MakeThickness.m script

tree = getenv("ROMEG");
load([tree '/Models/Real/head_model.mat'])

thickness = skull_thickness(p,t);

save([tree '/Models/Real/thickness.mat'],'thickness')