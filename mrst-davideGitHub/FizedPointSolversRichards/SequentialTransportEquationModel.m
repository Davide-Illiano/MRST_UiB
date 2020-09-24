classdef SequentialTransportEquationModel < ReservoirModel
    properties
        polymer
        Picard
        Newton 
        L_c
        L_p
        analyticalForce
        analyticalForce1
        analyticalForce2
    end
    
    methods
        function model = SequentialTransportEquationModel(G, rock, fluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            model.water = true;
            model.oil = false;
            model.gas = false;
            model.polymer = true;
            model.nonlinearTolerance = 1e-5;
            model = merge_options(model);
            model = merge_options(model,varargin{:});
        end
        
        function forces = getValidDrivingForces(model)

        forces = getValidDrivingForces@ReservoirModel(model);
        forces.analyticalForce = [];
        forces.analyticalForce1 = [];
        forces.analyticalForce2 = [];
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
           if model.Newton
            [problem, state] = equationsTransportSequential(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
           end
           
           if model.Picard
            [problem, state] = equationsTransportSequentialPicard(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
           end
           
           if model.L_c>0
                [problem, state] = equationsTransportSequentialLScheme(state0, state, model,...
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
       
        function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@ReservoirModel(model);
            if model.polymer
                names{end+1} = 'qWc';
            end
        end
       
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
        
            function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
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
     %}   
        end
  
end