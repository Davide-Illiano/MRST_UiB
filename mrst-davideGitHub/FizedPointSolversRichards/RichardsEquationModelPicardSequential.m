classdef RichardsEquationModelPicardSequential < TwoPhaseOilWaterModel
    properties
        polymer
    end
    
    methods
        function model = RichardsEquationModelPicardSequential(G, rock, fluid)
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            model.water = true;
            model.oil = false;
            model.gas = false;
            model.polymer = true;
            model = merge_options(model);
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsRichardsPicardSequential(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
        end

         function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'concentration', model.G.cells.num, 1);
        end

        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, drivingForces);
            if model.polymer
                c     = model.getProp(state, 'c');
            end
        end
        
        %%
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'c', 'concentration'}
                    fn = 'c';
                    index = 1;
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(model, name);
            end
        end

            function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % INPUT:
        %
        % model  - (Base class, automatic)
        %
        % cname  - Name of the component. Must be a property known to the
        %          model itself through getProp/getVariableField.
        %
        % eq     - Equation where the source terms are to be added. Should
        %          be one value per cell in the simulation grid (model.G)
        %          so that the src.sourceCells is meaningful.
        %
        % component - Cell-wise values of the component in question. Used
        %          for outflow source terms only.
        %
        % src    - Source struct containing fields for fluxes etc. Should
        %          be constructed from force and the current reservoir
        %          state by computeSourcesAndBoundaryConditionsAD.
        %
        % force  - Force struct used to produce src. Should contain the
        %          field defining the component in question, so that the
        %          inflow of the component through the boundary condition
        %          or source terms can accurately by estimated.
        if isempty(force)
            return
        end
        c = model.getProp(force, cname);
        cells = src.sourceCells;
        switch lower(cname)
            case {'concentration'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
            otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
        end
        eq(cells) = eq(cells) - qC;
        src.components{end+1} = qC;
    end
        
    end
  
end