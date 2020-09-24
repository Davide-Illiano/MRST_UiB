classdef RichardsEquationModelLScheme < RichardsEquationModel
    properties
        L
    end
    
    methods
        function model = RichardsEquationModelLScheme(G, rock, fluid, varargin)
            model = model@RichardsEquationModel(G, rock, fluid);
            
            model = merge_options(model, varargin{:});
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsRichardsLScheme(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
        end

    end
end