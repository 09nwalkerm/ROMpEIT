classdef RBModelClass < OrderedModelClass
    properties
        LF
    end
    methods
        function obj = RBModelClass(varargin)
            
            %obj = obj.processArgs(varargin);
            obj@OrderedModelClass(varargin)
            %obj = obj.getTOP();
        end

        function obj = readLF(obj,n)

            for kk=1:n
                try
                    load([obj.top '/Results/ROM/other/LF_EIT_' num2str(kk) '.mat'], 'ROM')
                    if size(ROM.V,1) > ROM.L
                        ROM.V = ROM.V(end-(ROM.L-1):end,:);
                    end
                catch
                    disp(['Could not load LF_EIT_' num2str(kk) '.mat'])
                    continue
                end
                obj.LF{kk}=ROM;
                clear ROM
                disp(['LF_EIT_' num2str(kk) ' loaded of ' num2str(n) '.'])
            end
            
        end

        function saveROM(obj)
        %
        % saveROM(top)
        %
        % Description:
        %   A function to save the reduce order model as RBModel.mat in the
        %   Results/ folder.
        %
        % Arguments:
        %   - top - Path to the top of the ROMEG tree
        %

            RBModel = obj;
            save([obj.top '/Results/ROM/RBModel.mat'],'RBModel')
            fprintf("\n \n The reduced order model has been saved in the Results/ folder as RBModel. \n \n")
            fprintf(" \n To run the inverse problem for this model: type `help GenInverse` into the command line. \n \n")
        end

        function runInverse(obj,varargin)
            GenInverse('ROM',true)
        end
        
    end
end