%% Bound Script

% Make a folder in storage pool called ROMEG_*_Bound, where * is a letter
% (either R or S) followed by a series of numbers to indicate layers
% (e.g.R123). 

%% Prelims

tree = getenv("ROMEG");
model = [tree '/Models/Real/head_model.mat']; % path to model to use

data = getenv("ROMEG_DATA");
jobid = getenv("SLURM_JOB_ID");
logger = log4m.getLogger([data '/logs/' jobid '.log']);
logger.setLogLevel(logger.INFO); % set to logger.OFF for only slurm log output
logger.setCommandWindowLevel(logger.INFO); % set to logger.OFF for only log file input

% Conds for spherical head model, comment out as necessary
%mu_max=[0.66,0.06,2.3,1.00,1.00,5]; % maximum conductivities
%mu_min=[0.15,0.001,1.1,0.05,0.05,5]; % minimum conductivities

% Conds for real head model, comment out as necessary
%mu_min=[0.15,0.001,0.001,1.1,0.05,0.05,5];
%mu_max=[0.66,0.006,0.060,2.3,1.00,0.65,5];

mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

num_samples=23;
num_start=1;
max_snap = 100;

%% Make measurements
% Type `help GenMeasurements` for examples, help and options

%GenMeasurements('model',model,'sample_num',num_start:num_samples,'mu_min',mu_min,...
%   'mu_max',mu_max,'noise',0.82e-6)

%% Run inverse solutions and error estimator values
% Type `help GenBound` for help and options

GenBound('sample_num',num_start:num_samples,'max_snap',100)

% Once ready, use the command line in MATLAB to load the bound.mat file 
% saved in the ROM folder and then run `BC.plotBound()`