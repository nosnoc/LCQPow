classdef LCQProblem < handle
    properties(Access=private)
        self % ptr to C++ LCQProblem object
    end

    methods
        % Constuctor
        obj = LCQProblem(nV, nC, nComp);

        %
        loadLCQP(obj, args);

        %
        runSolver();
    end
end
