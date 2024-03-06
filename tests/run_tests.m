import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.selectors.HasParameter;

mkdir("results")

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
