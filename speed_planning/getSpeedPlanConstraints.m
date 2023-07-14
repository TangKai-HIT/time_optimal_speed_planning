function Constraints = getSpeedPlanConstraints(Phi, Alpha, Mu, v_init, v_end)
% GETSPEEDPLANCONSTRAINTS return a Constraints struct for speed planning
%   Inputs:
%       Phi: 1 X dim -- speed constraints
%       Alpha: 1 X dim -- acceleration constraints
%       Mu: 1 X dim -- torque constraints
%       v_init: 1 X dim -- initial speed
%       v_end: 1 X dim -- end speed

Constraints.Phi = Phi;
Constraints.Alpha = Alpha;
Constraints.Mu = Mu;
Constraints.b_init = v_init.^2;
Constraints.b_end = v_end.^2;