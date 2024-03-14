import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.selectors.HasParameter;

% setting up environment
mkdir("results")
mkdir("results/logs")
mkdir("results/ROM")
addpath(genpath('functions'))
addpath(genpath('external'))
addpath(genpath('models'))
setenv("ROMEG_LOGS",'results/logs')
setenv("ROMEG_DATA",'results')
setenv("ROMEG_TOP",'results/ROM')

s = HasParameter('Property','model_shape','Name','Spherical') & ...
    HasParameter('Property','model_type','Name','head_model') & ...
    HasParameter('Property','mu_lims','Name','iso');

FOMSuite = TestSuite.fromClass(?FOMTest,s);

RBFSuite = TestSuite.fromClass(?RBFTest);

ROMSuite = TestSuite.fromClass(?ROMTest);

s = HasParameter('Property','snaps','Value',0) & ...
    HasParameter('Property','use_sinks','Value',1);

InverseSuite = TestSuite.fromClass(?InverseTest,s);

FullSuite = [FOMSuite RBFSuite ROMSuite InverseSuite];

results = FullSuite.run;

assertSuccess(results)
