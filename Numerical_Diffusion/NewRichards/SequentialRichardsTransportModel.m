classdef SequentialRichardsTransportModel < ReservoirModel
    % Sequential meta-model which solves Richards and transport equations using a fixed
    % flux splitting
    properties
        % Model for computing the pressure (Richards equation)
        RichardsEquationModel
        % Model for the transport subproblem after pressure is found
        TransportEquationModel
        
        % NonLinearSolver instance used for pressure updates
        RichardsNonLinearSolver
        % NonLinearSolver instance used for saturation/mole fraction updates
        TransportNonLinearSolver
        % Utility prop for setting pressure linear solver
        pressureLinearSolver
        % Utility prop for setting transport linear solver
        transportLinearSolver
        
        % Outer tolerance, which, if stepFunctionIsLinear is set to false,
        % is used to check if the pressure must be recomputed after
        % transport has been solved, in order to converge to the fully
        % implicit solution.
        outerTolerance
        % Indicates if we check well values when outer loop is enabled
        outerCheckWellConvergence
        % Maximum outer loops for a given step. When maxOuterIterations is
        % reached, the solver will act as if the step converged and
        % continue.
        maxOuterIterations
        % Update pressure based on new mobilities before proceeding to next
        % step
        reupdatePressure
        domainType
        forces
    end
    
    methods
        
        function forces = getValidDrivingForces(model)
        forces = getValidDrivingForces@ReservoirModel(model);
        forces.analyticalForce1 = [];
        forces.analyticalForce2 = [];
        end
        
        function model = SequentialRichardsTransportModel(RichardsEquationModel, TransportEquationModel, varargin)
            model = model@ReservoirModel([]);
            % Set up defaults
            model.RichardsEquationModel  = RichardsEquationModel;
            model.TransportEquationModel = TransportEquationModel;
            model.outerTolerance = 1e-6;  %-5
            model.outerCheckWellConvergence = false;
            model.maxOuterIterations = 20000;   % it was 2
            % Default: We do not use outer loop.
            model.stepFunctionIsLinear = false;   %if true we skip outer loop
            model.reupdatePressure = false; %true
            model = merge_options(model, varargin{:});
            
          %Transport model determines the active phases
           model.water = model.TransportEquationModel.water;   
           model.oil   = model.TransportEquationModel.oil; 
           model.gas   = model.TransportEquationModel.gas;
           
           
            if isempty(model.RichardsNonLinearSolver)
                model.RichardsNonLinearSolver = NonLinearSolver();    
            end
            model.RichardsNonLinearSolver.identifier = 'PRESSURE';
            
            if isempty(model.TransportNonLinearSolver)
                model.TransportNonLinearSolver = NonLinearSolver();    
            end
            model.TransportNonLinearSolver.identifier = 'TRANSPORT';

       
            model.RichardsNonLinearSolver.maxTimestepCuts = 10;   % maybe it was 0 
            model.RichardsNonLinearSolver.errorOnFailure = false;
     
            
            model.TransportNonLinearSolver.errorOnFailure = false;
            
        end
        
        
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            % Solve pressure and transport sequentially
            psolver = model.RichardsNonLinearSolver;
            csolver = model.TransportNonLinearSolver;
            
            % Get the forces used in the step
            forceArg = model.RichardsEquationModel.getDrivingForces(drivingForces);
            
            % First, solve the pressure using the pressure nonlinear
            % solver.
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.RichardsEquationModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_ok = pressureReport.Converged || psolver.continueOnFailure;
            
            if pressure_ok    
                if ~isempty(drivingForces.bc)
                    isDir = strcmpi(drivingForces.bc.type, 'pressure');
                    if any(isDir)
                        % Convert Dirichlet boundary conditions to flux
                        % boundary conditions for the transport
                        transportForces = drivingForces;
                        dirFace = transportForces.bc.face(isDir);
                        transportForces.bc.value(isDir) = sum(state.flux(dirFace, :), 2);
                        [transportForces.bc.type{isDir}] = deal('resflux');
                        forceArg = model.TransportEquationModel.getDrivingForces(transportForces);
                    end
                end
                state.timestep = dt;
                state.pressure_full = state.pressure;
                % If pressure converged, we proceed to solve the transport
                [state, transportReport] = ...
                    csolver.solveTimestep(state0, dt, model.TransportEquationModel,...
                                'initialGuess', state, ...
                                forceArg{:});
                transport_ok = transportReport.Converged;
            else
                transport_ok = false;
                transportReport = [];
            end
            

            converged = pressure_ok && transport_ok;
            if converged %&& ~model.stepFunctionIsLinear   %we want always the outer loop
                % Alternate mode: If outer loop is enabled, we will revisit
                % the pressue equation to verify that the equation remains
                % converged after the transport step. This check ensures
                % that the assumption of fixed total velocity is reasonable
                % up to some tolerance.
                %
                if isa(model.RichardsEquationModel, 'PressureNaturalVariablesModel')
                    values = max(abs(sum(state.s, 2) - 1));
                else
                    problem = model.RichardsEquationModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                    % Is the pressure still converged when accounting for the
                    % updated quantities after transport (mobility, density
                    % and so on?)
                    [~, values] = model.RichardsEquationModel.checkConvergence(problem);
                    values = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals(1);
                    values_p = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals(1); 
                    values_c = transportReport.StepReports{end}.NonlinearReport{end}.Residuals(1);
                    %values = [values_p, values_c];
                end
                
                if model.outerCheckWellConvergence
                    lv = max(values);
                else
                    values = values(1);
                    lv = values(1);
                end
               converged = all(values < model.outerTolerance);
               converged = converged || iteration > model.maxOuterIterations;
               
               %converged_p = all(values_p < model.outerTolerance);
               %converged_c = all(values_c < model.outerTolerance);
               %converged_p1 = all(values_p > 0);
               %converged_c1 = all(values_c > 0);
               converged = all(values_p < model.outerTolerance) && all(values_c < model.outerTolerance); % && converged_c1 && converged_p1;
               
               
                if model.verbose
                    if converged
                        s = 'Converged.';
                    else
                        s = 'Not converged.';
                    end
                    fprintf('OUTER LOOP step #%d with tolerance %1.4e: Largest value %1.4e -> %s \n', ...
                                                            iteration, model.outerTolerance, lv, s);
                end
            else 
                % Need to have some value
                values = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals(1);
            end
            
            failure = false;
            FailureMsg = '';
            if ~pressure_ok
                converged = converged && false;
            end
            report = model.makeStepReport(...
                                    'Failure',         failure, ...
                                    'Converged',       converged, ...
                                    'FailureMsg',      FailureMsg, ...
                                    'ResidualsConverged', converged, ...
                                    'Residuals',       values ...
                                    );
                                
            report.PressureSolver =  pressureReport;
            report.TransportSolver = transportReport;
            
            if model.reupdatePressure && converged
                state = ...
                    psolver.solveTimestep(state0, dt, model.RichardsEquationModel,...
                            'initialGuess', state, ...
                            forceArg{:});
                [~, state] = model.TransportEquationModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
            end
        end
        
        function varargout = getActivePhases(model)
            % Transport model solves for saturations, so that is where the
            % active phases are defined
            varargout = cell(1, nargout);
            [varargout{:}] = model.TransportEquationModel.getActivePhases();
        end
        
        function state = validateState(model, state)
            % Pressure comes first, so validate that.
            state = model.RichardsEquationModel.validateState(state);
        end

        function [model, state] = updateForChangedControls(model, state, forces)
            [model.RichardsEquationModel, state] = model.RichardsEquationModel.updateForChangedControls@PhysicalModel(state, forces);
        end

        function model = validateModel(model, varargin)
            model.RichardsEquationModel = model.RichardsEquationModel.validateModel(varargin{:});
            model.TransportEquationModel = model.TransportEquationModel.validateModel(varargin{:});
            return
        end

        function [fn, index] = getVariableField(model, name)
            [fn, index] = model.RichardsEquationModel.getVariableField(name);
        end
        
    end
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
