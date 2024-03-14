classdef ROMTest < matlab.unittest.TestCase
    
    properties
        FOM
        ROM
        ROM_expected
    end
    
    properties (ClassSetupParameter)
        electrode = {1};
        tolGREEDY = {5e-5};
        Nmax = {10};
    end
    
    properties (MethodSetupParameter)
    end
    
    properties (TestParameter)
    end
    
    methods (TestClassSetup)
        function classSetup(ROMTest,electrode,tolGREEDY,Nmax)
            load('tests/FOM.mat','FOM')
            FOM = FOM.startLogger();
            ROMTest.FOM=FOM;
            load('tests/sinks.mat','sinks')
            FOM.sinks = sinks;
            FOM.Fq{1}=sparse([FOM.np+electrode FOM.np+FOM.sinks(electrode,2:end)],ones(1,size(FOM.sinks(electrode,2:end),2)+1),[0.02e-3 -ones(1,size(FOM.sinks(electrode,2:end),2))*(0.02e-3)/size(FOM.sinks(electrode,2:end),2)],FOM.np+FOM.L,1);
            FOM.logger.debug('classSetup','built Fq in FOM object')
            ROMTest.ROM = ROMClass('sinks',sinks,'electrode',electrode,'FOM',FOM,'tolGREEDY',tolGREEDY,'Nmax',Nmax,'debug',true,'cmdLevel',2);
            FOM.logger.debug('classSetup','Initiated ROMClass')
            load('tests/RBModel.mat','RBModel')
            ROMTest.ROM_expected=RBModel;
        end
    end
    
    methods (TestMethodSetup)
        function methodSetup(ROMTest)
        end
    end
    
    methods (Test,ParameterCombination='sequential')
        function Greedytest(ROMTest)
            %ROMTest.ROM = ROMTest.ROM.startLogger();
            ROMTest.ROM.sample_grid = ROMTest.ROM_expected.LF{1}.sample_grid;
            ROMTest.ROM.Cqq = ROMTest.ROM_expected.LF{1}.Cqq;
            ROMTest.ROM = ROMTest.ROM.popFields();
            ROMTest.ROM=ROMTest.ROM.runGreedy(ROMTest.FOM);
            verifyEqual(ROMTest,size(ROMTest.ROM.V,2),size(ROMTest.ROM_expected.LF{1}.V,2))
            verifyLessThan(ROMTest,real(ROMTest.ROM.delta_Mean(end)),1e-1)
            disp('Verified greedy algorithm')
        end
        function testcalcXnorm(ROMTest)
            %TODO
        end
        function testPrelim(ROMTest)
            %TODO
        end
        function testSampleSpace(ROMTest)
            %TODO
        end
    end
end