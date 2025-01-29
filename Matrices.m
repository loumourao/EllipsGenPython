%% Matrices.m is a function that computes the quadratic matrix and the  
%% center for an ellipse defined by a parameter set.
%
% Given the parameter set E of the ellipse, this function calculates the 
% quadratic matrix Q and the center o of the ellipse. It also provides 
% additional outputs for the rotation matrix R and the scaling matrix D.
%
% -------------------------------------------------------------------------
% Inputs:
%   E - Parameter set of the ellipse [a, b, theta, o_x, o_y]
%     where:
%       a     - Semi-major axis
%       b     - Semi-minor axis
%       theta - Rotation angle of the ellipse
%       o_x   - x-coordinate of the center
%       o_y   - y-coordinate of the center
%
% Outputs:
%   Q - Quadratic matrix representing the ellipse
%   o - Center of the ellipse [o_x, o_y]
%   varargout:
%     {1} - Rotation matrix R for the ellipse
%     {2} - Scaling matrix D for the ellipse
%
% -------------------------------------------------------------------------

function [Q, o, varargout] = Matrices(E)

    % Extract parameters from the input
    a = E(1);      % Semi-major axis
    b = E(2);      % Semi-minor axis
    theta = E(3);  % Rotation angle
    o_x = E(4);    % x-coordinate of the center
    o_y = E(5);    % y-coordinate of the center
    
    % Define the rotation matrix for the given angle
    R = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    % Define the scaling matrix for the given semi-axes
    D = @(a, b) [1 / a^2, 0; 0, 1 / b^2];
    
    % Compute the quadratic matrix Q
    Q = R(theta) * D(a, b) * R(-theta);
    
    % Define the center of the ellipse
    o = [o_x; o_y];
    
    % Provide the rotation and scaling matrices as additional outputs
    varargout(1) = {R(theta)};
    varargout(2) = {D(a, b)};
end
