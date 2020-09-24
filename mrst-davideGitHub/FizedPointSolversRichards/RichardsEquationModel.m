classdef RichardsEquationModel < ReservoirModel
    properties
        forces
        Picard
        Newton
        LScheme
        L_p
        L_c
        Mixed
        
    end
    
    methods
        function model = RichardsEquationModel(G, rock, fluid,varargin)
            model = model@ReservoirModel(G, rock, fluid);
            model.water = true;
            model.oil = false;
            model.gas = false;
            %model.saturationVarNames = {'sw'};

            model = merge_options(model, varargin{:});
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsRichards(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
        end
        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@ReservoirModel(model, name);
            end
        end
    end
end