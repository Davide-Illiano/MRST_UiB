classdef RichardsEquationModelPicardScheme < RichardsEquationModel
    properties
        
    end
    
    methods
        function model = RichardsEquationModelPicardScheme(G, rock, fluid, varargin)
            model = model@RichardsEquationModel(G, rock, fluid);
            
            model = merge_options(model, varargin{:});
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsRichardsPicardScheme(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
        end

    end
end