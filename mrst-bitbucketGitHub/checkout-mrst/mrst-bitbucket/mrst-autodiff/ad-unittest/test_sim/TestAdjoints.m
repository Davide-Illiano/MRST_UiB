classdef TestAdjoints < matlab.unittest.TestCase
    properties
    end
    
    methods
        function test = TestAdjoints(varargin)
            mrstModule add ad-unittest ad-blackoil ad-core optimization
        end
        

        function testGradient(test, state0, schedule, model)
            G = model.G;
            % Gradient scaling
            scalFacs = struct('pressure', 100*barsa,...
                              'rate',     100/day);
            [wellSols, states] = simulateScheduleAD(state0, model, schedule);
            
            % This part is currently hard-written to use NPV for OW system
            [obj_adj, obj_num] = test.getObjectiveNPVOW(G, wellSols, schedule);
            
            numgrad = computeGradientPerturbationAD(state0, model, schedule, obj_num, 'scaling', scalFacs);
            adjgrad = computeGradientAdjointAD(state0, states, model, schedule, obj_adj);
            
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.IsEqualTo;
            reltol = RelativeTolerance(1e-2);
            
            test.assertLength(numgrad, numel(adjgrad));
            for i = 1:numel(numgrad)
                n = numgrad{i};
                g = adjgrad{i};
                test.assertLength(n, numel(g));
                % Compare each element by itself, because the magnitude
                % varies greatly
                for j = 1:numel(n)
                    test.verifyThat(n(j), IsEqualTo(g(j), 'Within', reltol));
                end
            end
        end
    end
    
    methods (Static)
        function [obj_adj, obj_num] = getObjectiveNPVOW(G, wellSols, schedule)
            % Setup function handles for NPV for oil/water system
            obj_adj = @(tstep)NPVOW(G, wellSols, schedule, ...
                            'ComputePartials', true, 'tStep', tstep);
            obj_num = @(wellSols, states, schedule)NPVOW(G, wellSols, schedule);
        end
    end
    
    methods (Test)
        function testGradientOW_NPV(test)
            [state0, schedule, model] = setupSimpleOW();
            test.testGradient(state0, schedule, model);
        end
    end
    methods
    end
end

