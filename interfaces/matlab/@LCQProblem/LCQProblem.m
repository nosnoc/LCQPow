classdef LCQProblem < handle
    properties%(Access=private)
        self % ptr to C++ LCQProblem object. On a 64 bit x86
    end

    methods
        % Constuctor
        function obj = LCQProblem(nV, nC, nComp)
            obj.constructProblem(nV, nC, nComp);
        end


        %
        loadLCQP(obj, args);

        %
        runSolver();

        function delete(obj)
        % obj is always scalar, if it isn't we panic :)
            obj.destructProblem();
        end
    end

    methods(Access=private)
        constructProblem(obj, nV, nC, nComp);
        destructProblem(obj);
    end
end
