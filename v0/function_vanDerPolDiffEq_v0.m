%Medical Physics Department at Bariloche Atomic Center, Argentina.
%Author: Osvaldo Velarde.
%Project: Differential equation for the Van der Pol Oscillator.
%Date: 10/06/2019.

% Van der Pol Equation.
% d^2 (x) - mu*(1-x^2)* d(x) + w^2 * x = 0
%
% The VdP equation is equivalent to the system
% dz/dt = F(z,mu,w)
% where z = (x,y) and F = (y, mu*(1-x^2)*y - w^2 * x)

% Inputs:
% t ->  Time (scalar).
% z ->  Variable (array 2 x 1).
% mu -> Nonlinearity parameter of the VdP equation (scalar).
% w ->  Angular frequency parameter of the VdP equation (scalar).

% Outputs:
% dzdt -> Derivative of vector "z" over time "t" (array 2 x 1).

%Reference:
%Velarde O, Urdapilleta E, Mato G, and Dellavale D (2019), Bifurcation
%structure determines different phase-amplitude coupling patterns in the
%activity of biologically plausible neural networks, NeuroImage, In Press,
%(DOI: ...)

function dzdt = function_vanDerPolDiffEq_v0(t,z,mu,w)
    dzdt = [z(2);...
            mu * (1 - z(1)^2) * z(2) - w^2 * z(1)];
end

