%% GenScript
% A script to run the full pipeline from generating a reduced order model
% from a head model to plotting a ROM vs traditional method comparison. To
% modify this pipeline please look at the documentation for each function
% to see the customizable options. To view the docs please use the 'help'
% command folled by the function of interest e.g 'help GenRBModel'.
%
% Please ensure that prior to starting matlab you have run the set_env.sh
% script (using 'source set_env.sh') while being in the setup/ folder. 
% This allows the appropriate environment variables to be setup.

%% Prelims
tree = getenv("ROMEG");
model = [tree '/Models/Spherical/head_model.mat']; % path to model to use

% Conds for spherical head model, comment out as necessary
mu_max=[0.66,0.060,2.3,1.00,5]; % maximum conductivities
mu_min=[0.15,0.001,1.1,0.05,5]; % minimum conductivities

% Conds for real head model, comment out as necessary
%mu_min=[0.15,0.001,0.001,1.1,0.05,0.05,5];
%mu_max=[0.66,0.006,0.060,2.3,1.00,0.65,5];

num_samples = 1;
c = [1,2];
%% Create sink file

%OrderedModelClass.patterns('model',model,'num_sinks',5,'elec_height',0.05,'electrodes',[7,65,62,19,52,56,10,103,16,110])

%% Make the measurements
% Type `help GenMeasurements` for examples, help and options

GenMeasurements('model',model,'num_samples',num_samples,'mu_min',mu_min,'mu_max',mu_max,'noise',0.82e-6)

%% Running the inverse snapshots for the snap plot
% Type `help GenInverse` for help and options

GenInverse('model',model,'ROM',true,'Cluster',true,'TRAD',true,'num_samples',num_samples,...
    'snaps',true,'active_layers',c,'use_sinks',true,'fix_conds',true)

%% Plotting

snap = SnapShotClass('range',200,'num_samples',num_samples,'split_conds',true);
snap = snap.readResults(c);
snap = snap.processResults();
snap = snap.saveSnap();

% Type `GenPlot('snap',true)` for plot