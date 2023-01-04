function d_state=jointTrajTrackModel(motionModel, state, jointStateTraj, t, t_sample)
% JOINTTRAJTRACKMODEL joint trajectory tracking error dynamic model for
% numerical simulation, using spline for state interpolation
%   Inputs:
%       motionModel: robot joint state motion model
%       state: actual current state
%       jointStateTraj: desire state trajectory to be followed -- [q; dq] (2*dim X  N)
%       t: current simulation time step
%       t_sample

% Use cubic-spline interpolation
desire_state = spline(t_sample, jointStateTraj, t);
    
% Compute state derivative
d_state = derivative(motionModel, state, desire_state);
