%% Configurations.m is a function that generates a pair of ellipses and 
%% computes their contact points
%
% -------------------------------------------------------------------------
% Inputs:
%   gamma_i, omega_i - Shape parameters of the first ellipse
%   theta_i - Orientation of the first ellipse
%   gamma_j - Shape parameter of the second ellipse
%   theta_j - Orientation of the second ellipse
%   o_j - Center of the second ellipse
%   phi - Contact angle parameter
%   epsilon_bar - Distance parameter, using the contact point on (i)
%
% Outputs:
%   E_i - Parameters of the first ellipse [a, b, theta, o_x, o_y]
%   E_j - Parameters of the second ellipse [a, b, theta, o_x, o_y]
%   x_i - Contact point on the first ellipse
%   x_j - Contact point on the second ellipse
%
% -------------------------------------------------------------------------

function [E_i, E_j, x_i, x_j, epsilon] = Configurations(gamma_i, ...
    omega_i, theta_i, gamma_j, theta_j, o_j, phi, epsilon_bar)

    % Rotation matrix function
    R = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
    % Diagonal matrix function for ellipse shape
    D = @(a, b) [1 / a^2, 0; 0, 1 / b^2];

    % Compute semi-axes of the ellipses
    a_i = sqrt(omega_i * gamma_i);
    b_i = sqrt(omega_i / gamma_i);
    a_j = sqrt(gamma_j);
    b_j = 1 / a_j;

    % Compute the transformed shape matrix of the first ellipse (i) in the 
    % second ellipse (j)'s frame
    Q_i_hat = D(1/sqrt(a_j), 1/sqrt(b_j)) * R(-theta_j) * R(theta_i) * ...
        D(a_i, b_i) * R(-theta_i) * R(theta_j) * D(1/sqrt(a_j), ...
        1/sqrt(b_j));
    [a_i_hat, b_i_hat, theta_i_hat] = parametrize(Q_i_hat);

    % Compute the initial contact point in the local frame of (i)
    x_i_bar = [a_i_hat * cos(phi); b_i_hat * sin(phi)];
    g_x_i_b = [b_i_hat * cos(phi); a_i_hat * sin(phi)];
    n_i_bar = g_x_i_b / norm(g_x_i_b); % Normal vector at the contact point
    o_j_bar = x_i_bar + (1 + epsilon_bar) * n_i_bar; % Second ellipse pos

    % Transformations to global coordinates
    T1 = @(x_bar) R(theta_i_hat) * (x_bar - o_j_bar);
    T2 = @(x_hat) R(theta_j) * D(1/sqrt(a_j), 1/sqrt(b_j)) * x_hat + o_j;

    % Compute center of the first ellipse in the transformed frame
    o_i_bar = [0; 0];
    o_i_hat = T1(o_i_bar);
    o_i = T2(o_i_hat);

    % Define ellipse parameters in global coordinates
    E_i = [a_i, b_i, theta_i, o_i(1), o_i(2)];
    E_j = [a_j, b_j, theta_j, o_j(1), o_j(2)];

    % Solve for contact points on the ellipses
    [x_i, x_j] = Contact_Solve(E_i, E_j);

    % Compute the separation distance
    epsilon = norm(x_j - x_i);
end
