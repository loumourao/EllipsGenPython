%% Contact_Solve.m is a function that computes the contact points between 
%% two ellipses.
%
% Given two ellipses defined by their parameter sets E_i and E_j, this 
% function finds the contact points x_i and x_j that satisfy the geometric 
% conditions.
%
% -------------------------------------------------------------------------
% Inputs:
%   E_i - Parameter set of the first ellipse [a, b, theta, o_x, o_y]
%   E_j - Parameter set of the second ellipse [a, b, theta, o_x, o_y]
%
% Outputs:
%   x_i - Contact point on the first ellipse
%   x_j - Contact point on the second ellipse
%   varargout:
%     {1} - Residual error for the first contact equation (fval1)
%     {2} - Residual error for the second contact equation (fval2)
%     {3} - Exit flag for the first solver (exitflag1)
%     {4} - Exit flag for the second solver (exitflag2)
%     {5} - Solver output details for the first contact point (output1)
%     {6} - Solver output details for the second contact point (output2)
%
% -------------------------------------------------------------------------

function [x_i, x_j, varargout] = Contact_Solve(E_i, E_j)

    % Compute quadratic matrices and centers of the ellipses
    [Q_i, o_i] = Matrices(E_i);
    [Q_j, o_j] = Matrices(E_j);
    
    % Define quadratic potential functions for each ellipse
    ei = @(x) x'*Q_i*x - 1;  % Ellipse i equation
    ej = @(x) x'*Q_j*x - 1;  % Ellipse j equation
    
    % Define the parametric form w(u) of co-gradient locus H_ij
    w = @(u) ((1 - u)*Q_i + u*Q_j) \ ((1 - u)*Q_i*o_i + u*Q_j*o_j);
    
    % Define contact conditions for each ellipse
    contact_point_1 = @(u) ei(w(u) - o_i);
    contact_point_2 = @(u) ej(w(u) - o_j);
    
    % Solve for the first contact point
    [parameter_1, fval1, exitflag1, output1] = fzero(contact_point_1, ...
        [0 1]);
    x_i = w(parameter_1);
    
    % Solve for the second contact point
    [parameter_2, fval2, exitflag2, output2] = fzero(contact_point_2, ...
        [0 1]);
    x_j = w(parameter_2);
    
    % Store additional output information
    varargout{1} = fval1;
    varargout{2} = fval2;
    varargout{3} = exitflag1;
    varargout{4} = exitflag2;
    varargout{5} = output1;
    varargout{6} = output2;
end
