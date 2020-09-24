classdef RichardsTransportEquationFixedPointSchemes < ReservoirModel
 %% Function used to solve the coupled probelm Richards and Transport equation. 
 %  If Picard set true we use Modified Picard Method
 %  RichardsTransportEquationModelPicardScheme(G, rock, fluid, 'Picard', 1);
 %
 %  To apply the L scheme specify the 2 constant 'L_p' and 'L_c' 
 %  RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'L_p', x, 'L_c', y)
    
    properties
        % Polymer present
        polymer
        
        forces
        Picard
        Newton
        NewtonDynamic
        LScheme
        L_p
        L_c
        Mixed
        n
        analyticalForce1
        analyticalForce2
        UseNewApp
        Advection
        K
        diffusion
    end

    
    methods
        function model = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, varargin) % no varagin
            model = model@ReservoirModel(G, rock, fluid);
            model.water = true;
            model.oil = false;
            model.gas = false;
          % This is the model parameters for oil/water/polymer
            model.polymer = true;
            
            model = merge_options(model,varargin{:}); 
        end
        
        function forces = getValidDrivingForces(model)
        forces = getValidDrivingForces@ReservoirModel(model);
        
        %% Analytical exteranl forces
        forces.analyticalForce1 = [];
        forces.analyticalForce2 = [];
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            if model.Picard
            [problem, state] = equationsRichardsTransportPicardScheme(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            end
            if model.Newton
             [problem, state] = equationsRichardsTransport(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});   
            end
            if model.Mixed
             [problem, state] = equationsRichardsTransportMixed(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});   
            end
            if model.LScheme
                [problem, state] = equationsRichardsTransportLScheme(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});  
            end
        end
        
         function state = validateState(model, state)
            state = validateState@ReservoirModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'concentration', model.G.cells.num, 1);
        end

        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            [state, report] = updateState@ReservoirModel(model, ...
               state, problem,  dx, drivingForces);
           if model.polymer
                c = model.getProp(state, 'c');
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.polymer
                c     = model.getProp(state, 'c');
            end
        end       
       
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'c', 'concentration'} %, 'polymer'
                    fn = 'c';
                    index = 1;
                case {'cmax'}
                    index = 1;
                    fn = 'cmax';
                case 'qwc'
                    index = 1;
                    fn = 'qWc';
                otherwise
                    
                    [fn, index] = getVariableField@ReservoirModel(model, name);
            end
        end
        
        function names = getComponentNames(model)
            names = getComponentNames@ReservoirModel(model);
            if model.polymer
                names{end+1} = 'concentration';
            end
        end
       
       function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@ReservoirModel(model);
            if model.polymer
                names{end+1} = 'cWells';
                types{end+1} = 'perf';
            end
       end
       
       %
       function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@ReservoirModel(model);
            if model.polymer
                names{end+1} = 'qWc';
            end
       end
       %{
        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ReservoirModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.polymer
                assert(model.water, 'Surfactant injection requires a water phase.');
                f = model.fluid;
                if well.isInjector
                    concWell = model.getProp(well.W, 'c');
                else
                    pix = strcmpi(model.getComponentNames(), 'c');
                    concWell = packed.components{pix};
                end
                qwsft = packed.extravars{strcmpi(packed.extravars_names, 'qwc')};
                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % get volume rate, at
                                            % surface condition.
                cqS = concWell.*cqWs;

                compEqs{end+1} = qwsft - sum(cqWs);
                compSrc{end+1} = cqS;
                compNames{end+1} = 'cWells';              
            end 
    end
       %}
      
       function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ReservoirModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.polymer
                assert(model.water, 'Polymer injection requires a water phase.');
                f = model.fluid;

                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition

                if well.isInjector
                    concWell = model.getProp(well.W, 'c');
                    cqP = concWell.*cqWs;
                else
                    pix = strcmpi(model.getComponentNames(), 'c');
                    concWell = packed.components{pix};

                    a = f.muWMult(f.cmax).^(1-f.mixPar);
                    cbarw     = concWell/f.cmax;

                    % the term (a + (1 - a).*cbarw) account for the
                    % todd-longstaff mixing factor, which model the fact that for
                    % not-fully mixed polymer solution the polymer does not
                    % travel at the same velocity as water. See the governing
                    % equation for polymer (e.g. equationsOilWaterPolymer.m)
                    cqP = concWell.*cqWs./(a + (1-a).*cbarw);
                end

                qwpoly = packed.extravars{strcmpi(packed.extravars_names, 'qwc')};

                compEqs{end+1} = qwpoly - sum(concWell.*cqWs);
                compSrc{end+1} = cqP;
                eqNames{end+1} = 'cWells';
            end
        end
        %}
        %
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
    %% Added diffusion/dispersion contribution
        if model.diffusion == 1
            N = model.G.faces.neighbors(force.face,:);
            BCcells = sum(N, 2);
            nbc = numel(force.face);

            cellToBCMap = sparse((1:nbc)', BCcells, 1, nbc, model.G.cells.num);
            isP = reshape(strcmpi(force.type, 'pressure'), [], 1);

            cBC   = cellToBCMap*component;
            dc = c - cBC(isP);
        end
        
        cells = src.sourceCells;
        switch lower(cname)
            case {'concentration'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
                if model.diffusion == 1
                    qC = (isInj.*c + ~isInj.*component(cells)).*qW + 1e-3.*dc;
                end
            otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
        end
        eq(cells) = eq(cells) - qC;
        src.components{end+1} = qC;
    end
     %}   
        end
  
end