function [eqs, state, src] = addBoundaryConditionsAndSourcesNEW(model, eqs, names, types, state, ...
                                                                 p, s, mob, rho, ...
                                                                 dissolved, components, ...
                                                                 forces,vW)
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
        %                `addComponentContributions`.
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
        
        % we need to introduce the diffusion effects od the transport
        % equation. We create a different computeSourcesAndBoundaryConditionsAD
        c = state.c; 
        [src_terms, bnd_cond, dC] = computeSourcesAndBoundaryConditionsandConcAD(model, p, c, s, mob, rho, dissolved, forces);
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

                bc = bnd_cond.sourceCells;
                if ~isempty(bc)
                    eqs{sub}(bc) = eqs{sub}(bc) - bnd_cond.phaseMass{i}./rhoS(i);
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
            [eq, bnd_cond, newvW] = model.addComponentContributions(name, eq, C, bnd_cond, forces.bc, dC, vW);
            [eq, src_terms] = model.addComponentContributions(name, eq, C, src_terms, forces.src);
            
            %state.NewvW = double(newvW);
%             state = model.storeFluxes(state, newvW, [], []);
               
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
