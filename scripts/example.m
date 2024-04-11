%% Example script
% This script offers a lightweight example of the ROMpEIT toolbox and the
% general workflow to get some results quickly. For more information and
% a link to the documentation please see the README in the root of the
% repo. This example is designed to be run on a single PC as it is
% illustrative. The real benefits of ROMpEIT become more apparent when the
% head models are more complex, there more compartments and more injection 
% patterns. However, this requires access to a compute cluster.

%% Set environment: run set_env.sh first
% Please make a directory for the results and then run the set_env.sh
% script OR run the set_env.sh script and leave the prompt blank for the
% test results and logs to go into the results folder.

%% setting parameters
mu_min=[0.303,0.002,1.7,0.3,5]; % minimum conductivities
mu_max=[0.444,0.043,1.7,0.3,5]; % maximum conductivities

% setting the path to the model in the repo
% feel free to change to your own model
tree = getenv("ROMEG");
model = [tree '/models/Spherical/head_model.mat'];

%% Make the elctrode injection pattern
% Creating and saving a file called sinks.mat that describes the injection 
% and extraction electrodes used to train ROM. This creates only 3 patterns
% for training where the injections are electrodes 63, 64 and 65. All the
% sinks are set to the final electrode (120 in this case).

OrderedModelClass.patterns('model',model,'electrodes',[63,64,65])

%% Make the Reduced Order Model
% Creating the Reduced Order Model by training within the min and max
% conductivity ranges. The 'nic' value is the number of bound support
% points and Nmax is the number of snapshots. We'll set debug to true so
% that we can check the logs if something goes wrong.

GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'nic',10,...
    'Nmax',10,'current',0.02e-3,'debug',true)

%% Generate some noisey synthetic measurements
% We need some sythetic measurements to test against. This uses the model
% to simulate the behaviour of the current and then adds noise to the
% measurements to simulate real world use.

GenMeasurements('model',model,'noise',0.82e-6,'synth_cond',...
    [0.4,0.009,1.7,0.3,5],'debug',true)

%% Run the Inverse Problem for ROMpEIT.
% Finally we can use our Reduced Model to run the inverse problem for pEIT.
% We have to tell the function to use the same electrode patterns with
% 'use_sinks' and also for it to read the noisey measurements with
% 'use_noise'. Notice that it completes almost instantly. We also specify
% the active layers to estimate, which is set to the first 2 compartments.

GenInverse('model',model,'ROM',true,'current',0.02e-3,...
    'use_sinks',true,'use_noise',true,'debug',true,...
    'active_layers',[1 2],'ground',1)

% Try adding the "'simultaneous',true" pair to the GenInverse function to
% see all the electrode pairs be estimated together in one optimisation.