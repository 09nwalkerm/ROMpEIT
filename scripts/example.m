%% Example script

% run set_env.sh first


mu_min=[0.15,0.001,1.1,0.05,5]; % minimum conductivities
mu_max=[0.66,0.06,2.3,1.00,5]; % maximum conductivities
model = '/home/c1616132/Documents/PhD/ROMEG/Models/Spherical/head_model.mat';

OrderedModelClass.patterns('model',model,'electrodes',[63,64,65])

GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'nic',40,...
    'Nmax',10,'current',0.02e-3,'debug',true)

GenMeasurements()

GenInverse()