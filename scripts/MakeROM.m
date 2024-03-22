%% Generate a RBModel

% Script to determine which layers are best to estimate by the reduction in
% overall estimation error.

%% Prelims
tree = getenv("ROMEG");
model = [tree '/models/Real/head_model_theta.mat']; % path to model to use

%data = getenv("ROMEG_DATA");
%jobid = getenv("SLURM_JOB_ID");
%logger = log4m.getLogger([data '/logs/' jobid '.log']);
%logger.setLogLevel(logger.INFO); % set to logger.OFF for only slurm log output
%logger.setCommandWindowLevel(logger.INFO); % set to logger.OFF for only log file input

% Conds for spherical head model, comment out as necessary
%mu_max=[0.66,0.060,2.3,1.00,5]; % maximum conductivities
%mu_min=[0.15,0.001,1.1,0.05,5]; % minimum conductivities

% Conds for real head model, comment out as necessary
%mu_min=[0.15,0.001,0.001,1.1,0.05,0.05,5];
%mu_max=[0.66,0.006,0.060,2.3,1.00,0.65,5];

%mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
%mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

mu_min = [0.303,0.002,0.002,1.450,0.268,0.092,5];
mu_max = [0.444,0.043,0.043,1.794,0.508,0.177,5];

%mu_min = [0.303,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,1.794,0.508,0.177,5];

%% Create sink file
% Type `help OrderedModelClass.patterns` for help, options and examples

%OrderedModelClass.patterns('model',model)

%% Run functions to make RB model
% Type `help GenRBModel` for help, options and examples

GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'tolGREEDY',5e-10,'Nmax',100,...
    'use_sinks',true,'nic',150,'anis_rad',2,'anis_tan',3,'angles',true,'debug',true,...
    'Cluster',true,'current',0.02e-3,'use_FOM',true)