classdef RBFTest < matlab.unittest.TestCase
    
    properties
        FOM
        FOM_tmp
        beta
        beta_int
    end
    
    properties (TestParameter)
        mu_test = {[0.3,0.02,1.9,0.3],[0.6,0.002,1.6,0.3]};
    end
    
    methods (TestClassSetup)
        function loadFOM(obj)
            load("tests/FOM.mat",'FOM')
            obj.FOM = FOM;
            obj.FOM = obj.FOM.startLogger();
            obj.FOM.logger.info('loadFOM','Logger on')
        end
    end
    
    methods (Test,ParameterCombination='sequential')
        function Eigtest(obj)
            %obj.FOM.startLogger();
            %obj.FOM.logger.info('Eigtest','Logger on')
            [obj.beta,~]=femeg_ROM_RBF_offline_dual_iter(obj.FOM,1);
            disp(['Beta from Eig func -- ' num2str(obj.beta)])
            verifyEqual(obj,obj.beta,obj.FOM.betaa(1),"AbsTol",1e-6)
            disp('Eigenvalue calculation verified')
        end
        function Coefstest(obj)
            obj.FOM_tmp=femeg_ROM_RRBF(obj.FOM);
            verifyEqual(obj,obj.FOM_tmp.RRBF,obj.FOM.RRBF,"RelTol",0.1*ones(length(obj.FOM.RRBF),1))
            disp('Interpolation coefs verified')
        end
        function Interptest(obj,mu_test)
            obj.beta_int = femeg_ROM_RBF_online(mu_test,obj.FOM);
            disp(['Interpolated beta -- ' num2str(obj.beta_int)])
            verifyGreaterThan(obj,obj.beta_int,1e-6)
            verifyLessThan(obj,obj.beta_int,9e-6)
            disp('Interpolation verified')
        end
    end
end