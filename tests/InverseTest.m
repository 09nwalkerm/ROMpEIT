classdef InverseTest < matlab.unittest.TestCase
    
    properties
        IROM
    end
    
    properties (ClassSetupParameter)
        snaps = {1,0};
        use_sinks = {1,0};
        min_snap = {5};
    end
    
    properties (MethodSetupParameter)
    end
    
    properties (TestParameter)
    end
    
    methods (TestClassSetup)
        function classSetup(INVTest,snaps,min_snap,use_sinks)
            load('tests/RBModel.mat','RBModel')
            load('tests/sinks.mat','sinks')
            INVTest.IROM = InverseROMClass('injection',1,'sinks',sinks,'snaps',snaps, ...
                'use_sinks',use_sinks,'min_snap',min_snap,'eL',120);
            INVTest.IROM.LF = RBModel.LF;
        end
    end
    
    methods (TestMethodSetup)
        function methodSetup(INVTest)
            INVTest.IROM.lf = INVTest.IROM.LF{1}.non_active;
            INVTest.IROM.te = INVTest.IROM.LF{1}.active;
            load('tests/measurement.mat','Data')            
            INVTest.IROM.cond_lf = Data.synth_cond(INVTest.IROM.lf);
            INVTest.IROM.lb = INVTest.IROM.LF{1}.mu_min(INVTest.IROM.te);
            INVTest.IROM.ub = INVTest.IROM.LF{1}.mu_max(INVTest.IROM.te);
            INVTest.IROM.x0 = (INVTest.IROM.lb + INVTest.IROM.ub)/2;
            INVTest.IROM.u{INVTest.IROM.injection} = Data.u;
            INVTest.IROM.el_in = INVTest.IROM.injection;
            INVTest.IROM.eL = INVTest.IROM.LF{1}.L;
            INVTest.IROM.synth_cond = Data.synth_cond;
        end
    end
    
    methods (Test,ParameterCombination='sequential')
        function InvROM(INVTest)
            INVTest.IROM = INVTest.IROM.opt();
            load('tests/estimate.mat','estimate')
            verifyEqual(INVTest,INVTest.IROM.estimate,estimate,"RelTol",ones(size(estimate)))
        end
    end
end