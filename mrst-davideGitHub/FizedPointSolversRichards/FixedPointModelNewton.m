%% Fixed point method based on Modified Picard
classdef FixedPointModelNewton < PhysicalModel
    properties
        parentModel
    end
    
    methods
        function model = FixedPointModelNewton(parentModel, varargin)
            model = model@PhysicalModel(parentModel.G);
            model = merge_options(model, varargin{:});
            
            model.parentModel = parentModel;
        end
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            %% 
            resOnly = false;
            [problem, state] = model.parentModel.getEquations(state0, state, dt, drivingForces, ...
                                       'ResOnly', resOnly, ...
                                       'iteration', iteration, ...
                                       varargin{:});

            problem.iterationNo = iteration;
            problem.drivingForces = drivingForces;
            [convergence, values, resnames] = model.parentModel.checkConvergence(problem);

            %% Assemble system
            problem = problem.assembleSystem();
                % Store number of primary variables         
                ix = find(cellfun(@(x) isa(x, 'ADI'), problem.equations), 1);
                state.numVars = cellfun(@(x) size(x, 2), problem.equations{ix}.jac)';
               
                J = -problem.A;
                F = problem.b;

            % The update of the system is obtained with the expression
            dx =  -J\F;
             
            % Store output
            numVars = state.numVars;
            cumVars = cumsum(numVars);
            
            ii = [[1;cumVars(1:end-1)+1], cumVars];
            % Store increments for each primary variable
            eqn = size(ii,1);
            deltax = cell(eqn,1);
            for i = 1:eqn
                deltax{i} = dx(ii(i,1):ii(i,2), :);
            end
            
            doneMinIts = iteration > nonlinsolve.minIterations;
            isConverged = (all(convergence) && doneMinIts) || model.parentModel.stepFunctionIsLinear;
            if model.parentModel.verbose
                printConvergenceReport(resnames, values, convergence, iteration);
            end
[state, updateReport] = model.parentModel.updateState(state, problem, deltax, drivingForces);
            
            report = model.parentModel.makeStepReport(...
                            'UpdateState',  updateReport, ...
                            'Failure',      false, ...
                            'FailureMsg',   '', ...
                            'Converged',    isConverged, ...
                            'Residuals',    values, ...
                            'ResidualsConverged', convergence);
            %figure(1); clf
            %plot(state.pressure)
            %drawnow
                     
        end
        
        function varargout = getActivePhases(model)
            varargout = cell(1, nargout);
            [varargout{:}] = model.parentModel.getActivePhases();
        end
        function forces = getValidDrivingForces(model)
            forces = model.parentModel.getValidDrivingForces();
        end
        function state = validateState(model, state)
            state = model.parentModel.validateState(state);
        end

        function [model, state] = updateForChangedControls(model, state, forces)
            [model.parentModel, state] = model.parentModel.updateForChangedControls(state, forces);
        end

        function model = validateModel(model, varargin)
            model.parentModel = model.parentModel.validateModel(varargin{:});
        end

        function [fn, index] = getVariableField(model, name)
            [fn, index] = model.parentModel.getVariableField(name);
        end
        
        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = model.parentModel.updateAfterConvergence(varargin{:});
        end
    end
end