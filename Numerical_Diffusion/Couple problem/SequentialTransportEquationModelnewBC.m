classdef SequentialTransportEquationModelnewBC < ReservoirModel
    properties
        polymer
        Picard
        Newton 
        L_c
        L_p
        forces
        D
    end
    
    methods
        function model = SequentialTransportEquationModelnewBC(G, rock, fluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            model.water = true;
            model.oil = false;
            model.gas = false;
            model.polymer = true;
            model = merge_options(model);
            model.nonlinearTolerance = 1e-5;
            model = merge_options(model,varargin{:});
        end
        
        function forces = getValidDrivingForces(model)

        forces = getValidDrivingForces@ReservoirModel(model);
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
           if model.Newton
            [problem, state] = equationsTransportSequential(state0, state, model,...
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
        
       function [eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 p, s, mob, rho, ...
                                                                 dissolved, components, ...
                                                                 forces)
        % Add in the boundary conditions and source terms to equations
        %
        % SYNOPSIS:
        %   [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
        %                                                        p, sat, mob, rho, ...
        %                                                        rs, components, ...
        %                                                        drivingForces);
        %
        % PARAMETERS:
        %   model  - Class instance.
        %   eqs    - Cell array of equations that are to be updated.
        %
        %   names  - The names of the equations to be updated. If
        %            phase-pseudocomponents are to be used, the names must
        %            correspond to some combination of "water", "oil", "gas"
        %            if no special component treatment is to be introduced.
        %
        %   types  - Cell array with the types of "eqs". Note that these
        %            types must be 'cell' where source terms is to be added.
        %
        %   src    - Struct containing all the different source terms that
        %            were computed and added to the equations.
        %
        %   p      - Cell array of phase pressures.
        %
        %   s      - Cell array of phase saturations.
        %
        %   mob    - Cell array of phase mobilities
        %
        %   rho    - Cell array of phase densities
        %
        %   dissolved - Cell array of dissolved components for black-oil
        %               style pseudocompositional models.
        %
        %   components - Cell array of equal length to the number of
        %                components. The exact representation may vary
        %                based on the model, but the respective
        %                sub-component is passed onto
        %                `addComponentConbctributions`.
        %
        %   forces - DrivingForces struct (see `getValidDrivingForces`)
        %            containing (possibily empty) `src` and `bc` fields.
        %
        % RETURNS:
        %   eqs   - Equations with corresponding source terms added.
        %   state - Reservoir state. Can be modified to store e.g. boundary
        %           fluxes due to boundary conditions.
        %   src   - Normalized struct containing the source terms used.
        %
        % NOTE:
        %  This function accomodates both the option of black-oil
        %  pseudocomponents (if the model equations are named "oil", "gas"
        %  or water) and true components existing in multiple phases.
        %  Mixing the two behaviors can lead to unexpected source terms.
        
        %[src_terms, bnd_cond] = computeSourcesAndBoundaryConditionsAD(model, p, s, mob, rho, dissolved, forces);
        c = state.c;
        [src_terms, bnd_cond, dC, ez_contribution] = computeSourcesAndBoundaryConditionsandConcAD(model, p, c, s, mob, rho, dissolved, forces);        
        [~, longNames] = getPhaseNames(model);
        rhoS = model.getSurfaceDensities();
        % We first consider pseudocomponents that correspond to phases,
        % e.g. the black-oil model, most immiscible models and other models
        % where the number of phases is approximately equal to the number
        % of components. If the equations "water", "oil" and "gas" exist,
        % these will get direct source terms added. Note that the
        % corresponding source terms will already have added the effect of
        % dissolution (rs/rv).
        %
        % For fully compositional problems, this branch will not execute.
        for i = 1:numel(s)
            sub = strcmpi(names, longNames{i});
            if any(sub)
                assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
                sc = src_terms.sourceCells;
                if ~isempty(sc)
                    eqs{sub}(sc) = eqs{sub}(sc) - src_terms.phaseMass{i}./rhoS(i);
                end
                
                if isfield(forces.bc, 'phaseMass')
                    bnd_cond.phaseMass{i} = forces.bc.phaseMass(:, i);
                end

                bc = bnd_cond.sourceCells;
                if ~isempty(bc)
                    if isempty(bnd_cond.mapping)
                        q = bnd_cond.phaseMass{i}./rhoS(i);
                    else
                        q = (bnd_cond.mapping*bnd_cond.phaseMass{i})./rhoS(i);
                    end
                    eqs{sub}(bc) = eqs{sub}(bc) - q;
                end
            end
        end
        % Get the fluxes and store them in the state.
        if nargout > 1 && model.outputFluxes
            act = model.getActivePhases();
            tmp = cell(numel(act), 1);
            tmp(act) = bnd_cond.phaseVolume;
            state = model.storeBoundaryFluxes(state, tmp{:}, forces);
        end
        % Finally deal with actual components that exist in the different
        % phases to varying degrees.
        cnames = model.getComponentNames();
        for i = 1:numel(cnames)
            % Iterate over individual components
            name = cnames{i};
            sub = strcmpi(name, names);
            if any(sub)
                eq = eqs{sub};
            else
                eq = zeros(model.G.cells.num, 1);
            end
            C = components{i};
            assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
            % Add BC component source terms
%             
            [eq, bnd_cond] = model.addComponentContributions(name, eq, C, bnd_cond, forces.bc, dC, ez_contribution, boundary_info);
            [eq, src_terms] = model.addComponentContributions(name, eq, C, src_terms, forces.src);
            if any(sub)
                eqs{sub} = eq;
            end
        end
        % If requested, provide the computed values for source and bc for
        % further manipulations outside this function.
        if nargout > 2
            src = struct('src', src_terms, 'bc', bnd_cond);
        end
    end
       
       
       
       
       function [eq, src, new_vw] = addComponentContributions(model, cname, eq, component, src, force, dC, boundary_info)
                
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
                
                %% Change adding actual vw not qw. Also have to add gravity contribuiton at boundary.
                % basically import K and ez now we only have vw.
                vW = boundary_info.vW;
                VW = boundary_info.vW;
                T = boundary_info.T; 
                T_all = boundary_info.T_all;
                p = boundary_info.p;
                bfaces  = find(any(model.G.faces.neighbors==0,2));
                X = model.G.cartDims(1);
                Y = model.G.cartDims(2);
                top_int = bfaces(4*X:5*X-1);    % TOP int
                bottom_int = bfaces(3*X:4*X-1); % BOTTOM int
%
N = model.G.faces.neighbors(force.face,:);                
BCcells = sum(N, 2);
nbc = numel(force.face);
cellToBCMap = sparse((1:nbc)', BCcells, 1, nbc, model.G.cells.num);

        ez =  model.G.cells.centroids(:,2);
        yTop = max(model.G.nodes.coords(:,2));
        ez_val = [yTop * ones((model.G.cartDims(1)),1)' zeros((model.G.cartDims(1)),1)' ]';
        dez = ez_val - cellToBCMap*ez;
%         T = [T(size(VW.val,1)-X+1:size(VW.val,1))' T(size(VW.val,1)/2-X-4+X:size(VW.val,1)/2-5+X)' ]'; %         
        T = [ T(size(T,1)/2-X-4+X:size(T,1)/2-5+X)' T((size(T,1))-X+1:size(T,1))']'; 
%         d_ez(1:X) = 0;
        d_ez = T.*dez;
        
        %faceFlux2cellVelocity(model.G, vW.val)
                
%                 vW = [vW.val(1:X)' vW.val((size(vW.val,1))-X+1:size(vW.val,1))']';
                vW = [ vW.val(size(vW.val,1)-X+1:size(vW.val,1))' vW.val(size(vW.val,1)/2-X-4+X:size(vW.val,1)/2-5+X)']';
                vW = vW;
%                 vW(1:X) = 0*ones(X,1);
%                 vW(X+1:size(vW,1)) = -0.0334 * ones(X,1);
%                 vW = -T.*dp;
                isInj = vW > 0;
%                 qC = (isInj.*c + ~isInj.*component(cells)).*vW;
%                 qC = c .* vW;
                %}
                
                %% Itroduce also contribution from dispersion
                qC = qC + 1e-3.*dC; 
                
                %update full flux
%                 boundary_info.vW( size(VW.val,1)-X+1:size(VW.val,1)  ) = vW( 1:(model.G.cartDims(1)) ); 
%                 boundary_info.vW( size(VW.val,1)/2-X-4+X:size(VW.val,1)/2-5+X ) = -0.033* ones(X,1 );
%                 boundary_info.vW( size(boundary_info.vW.val,1)-X+1:size(boundary_info.vW.val,1) ) = -0.033* ones(X,1 );

%                 
%{
%                  boundary_info.vW( 1:X  ) = vW( 1:(model.G.cartDims(1)) ); 
                boundary_info.vW( size(boundary_info.vW.val,1)/2-X-4+X:size(boundary_info.vW.val,1)/2-5+X ) = vW( (size(vW,1))-model.G.cartDims(1)+1:size(vW,1) );
                boundary_info.vW( size(boundary_info.vW.val,1)-X+1:size(boundary_info.vW.val,1) ) = vW( (size(vW,1))-model.G.cartDims(1)+1:size(vW,1) );
% %                  
%                 boundary_info.vW( size(boundary_info.vW.val,1)/2-X-4+X:size(boundary_info.vW.val,1)/2-5+X ) = mean(mean(boundary_info.vW.val)).*ones(X,1)/4.92;
%                 boundary_info.vW( size(boundary_info.vW.val,1)-X+1:size(boundary_info.vW.val,1) ) = mean(mean(boundary_info.vW.val)).*ones(X,1)/4.92;
% %                 
%                 VW = boundary_info.vW;
%                 VW = 0 .* VW;
%                 VW(size(VW.val,1)-X+1:size(VW.val,1)) = ones(X,1)/4.962;% top vy
%                 VW(size(VW.val,1)/2-X-4+X:size(VW.val,1)/2-5+X) = ones(X,1)/4.962; % bottom vy
%                 boundary_info.vW = boundary_info.vW.*0;
%                 boundary_info.vW(1:2*X) = ones(X,1);
%                 boundary_info.vW = vW;
              %}
                new_vw = boundary_info.vW; 
                
            otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
        end
        eq(cells) = eq(cells) - qC;
        src.components{end+1} = qC;
    end
     %}   
        end
  
end