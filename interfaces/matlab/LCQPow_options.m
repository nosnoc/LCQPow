function [ options ] = LCQPow_options( )
% setup options struct with default values
	  options = struct(	'stationarityTolerance',    1.0e3*eps, ...
        'complementarityTolerance', 1.0e6*eps, ...
        'initialPenaltyParameter',  0.01, ...
        'penaltyUpdateFactor',      2.0, ...
        'solveZeroPenaltyFirst',    true, ...
        'perturbStep',              true, ...
        'maxIterations',            1000, ...
        'maxPenaltyParameter',      1.0e8, ...
        'printLevel',               1, ...
        'nDynamicPenalty',          3, ...
        'etaDynamicPenalty',        0.9);
end
