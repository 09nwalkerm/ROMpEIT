%% All Plots Script

% This script generates, processes and saves all the required data for the
% sensitivity analysis, snapshot plot and bound plot in the paper. To run
% this script it is reccommended you use the SubmitScript.sh sbatch script
% so this can run on a local cluster. Prior to running this you will need
% to setup the environment using the setup/setenv script and set the
% ROMEG_DATA environment variable to a directory with an RBModel in it.
% 
% See the MakeROM script to create the RBModel needed for this script.

%% Prelims

tree = getenv("ROMEG");
model = [tree '/Models/Real/head_model_theta.mat']; % path to model to use

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

mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

% Anisotropic values in the skull
%mu_min = [0.303,0.002,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,0.043,1.794,0.508,0.177,5];

% For merged skull and spongiform bone
mu_min = [0.303,0.002,1.450,0.268,0.092,5];
mu_max = [0.444,0.043,1.794,0.508,0.177,5];

num_samples=50;
num_start=1;

c = [1 2 3 4 5];
c2 = [1 2 3];
%layer_names = {'Scalp','Compact Skull','Spongiform bone','CSF','Grey Matter','White Matter'};
layer_names = {'Scalp','Radial Skull','Tangential Skull','CSF','Grey Matter','White Matter'};
%layer_names = {'Scalp','Skull','CSF','Brain'};

%% Make measurements
% Type `help GenMeasurements` for examples, help and options

GenMeasurements('model',model,'sample_num',num_start:num_samples,'mu_min',mu_min,...
  'mu_max',mu_max,'noise',0.82e-6,'anis_rad',2,'anis_tan',3,'angles',true,'ratio',3)

%% Run inverse solutions and error estimator values
% Type `help GenBound` for help and options

GenBound('sample_num',num_start:num_samples,'max_snap',100)

%% Running the inverse problem for sensitivity
% Type `help GenInverse` for help and options 

GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
  'active_layers',c,'sensitivity',true,'use_sinks',true,'use_noise',true)

sense = SenseClass('layers',c,'layer_names',layer_names,'plot','Box','num_samples',num_start:num_samples);
sense = sense.readResults(c);
sense = sense.processResults(c);
sense = sense.saveSense();

%% Running the inverse snapshots for the snap plot
% Type `help GenInverse` for help and options

GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
 'snaps',true,'active_layers',c,'use_sinks',true,'fix_conds',true,'use_noise',true)

GenInverse('model',model,'TRAD',true,'sample_num',num_start:num_samples,...
  'snaps',true,'active_layers',c2,'use_sinks',true,'use_noise',true)

% snap = SnapShotClass('range',100,'num_samples',num_samples,'tissue',[1 2 3]);
% snap = snap.readResults(c,c2);
% snap = snap.processResults(c2);
% snap = snap.saveSnap();

%% Plotting

% Type `GenPlot('sense',true,'snap',true,'bound',true)` to see the plots