classdef FOMTest < matlab.unittest.TestCase

    properties
        FOM_test
        FOM_expected
    end
    
    properties (ClassSetupParameter)
        model_shape = {'Spherical','Real'};
        model_type = {'head_model','head_model_theta'};
    end
    
    properties (MethodSetupParameter)
        mu_lims = struct('iso',[0.15,0.001,1.1,0.05,5;0.66,0.06,2.3,1.00,5],...
            'anis',[0.15,0.001,0.001,1.1,0.05,5;0.66,0.06,0.06,2.3,1.00,5]);
    end
    
    properties (TestParameter)
        z = {1};
        q = {3};
    end

    methods (TestClassSetup)
        function classSetup(FOMTest,model_shape,model_type)
            model = ['models/' model_shape '/' model_type '.mat'];
            FOM = FOMClass('model',model,'nic',20);
            FOMTest.FOM_test=FOM;
            FOMTest.FOM_test = FOMTest.FOM_test.assembleFOM();
            load("tests/FOM.mat",'FOM')
            FOMTest.FOM_expected = FOM;
        end
    end
    
    methods (TestClassTeardown)
        function closefig(FOMTest)
        end
        
        function closeotherfig(FOMTest)
        end
    end
    
    methods (TestMethodSetup)
        function methodSetup(FOMTest,mu_lims)
            FOMTest.FOM_test.mu_min = mu_lims(1,:);
            FOMTest.FOM_test.mu_max = mu_lims(2,:);
        end
    end
    
    methods (TestMethodTeardown)
        function closetest(FOMTest)
        end
    end

    methods (Test, ParameterCombination = 'sequential')
        function testAssembleFOM(FOMTest)
            FOMTest.FOM_test = FOMTest.FOM_test.assembleFOM();
            verifyEqual(FOMTest,FOMTest.FOM_test.Qa,size(FOMTest.FOM_test.mu_min,2)-1)
            disp('Verified FOM.assembleFOM()')
        end
        function testCompStiff(FOMTest)
            FOMTest.FOM_test = FOMTest.FOM_test.computeStiff();
            verifyEqual(FOMTest,size(FOMTest.FOM_test.Aq{1}),size(FOMTest.FOM_expected.Aq{1}))
            disp('Verified FOM.computeStiff()')
        end
        function testMuFOM(FOMTest,q)
            FOMTest.FOM_test = FOMTest.FOM_test.computeMuTrain();
            verifyEqual(FOMTest,length(FOMTest.FOM_test.mu_train),length(FOMTest.FOM_expected.mu_train))
            disp('Verified FOM.computeMuTrain()')
        end
        function testXnormFOM(FOMTest,z)
            FOMTest.FOM_test = FOMTest.FOM_test.computeXnorm();
            verifyEqual(FOMTest,size(FOMTest.FOM_test.Xnorm),size(FOMTest.FOM_expected.Xnorm))
            disp('Verified FOM.computeXnorm()')
        end
    end
        
    methods (Test, ParameterCombination = 'pairwise')
        function testcomputeAngles(FOMTest,z,q)
            %write test
        end
    end
end