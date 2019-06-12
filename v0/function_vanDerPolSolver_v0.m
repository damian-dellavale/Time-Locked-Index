%Medical Physics Department at Bariloche Atomic Center, Argentina.
%Author: Osvaldo Velarde (osva.m.velarde@gmail.com).
%Project: Solver for the Van der Pol differential equation.
%Date: 10/06/2019.

% Solutions of Van der Pol equation.
% d^2 (x) - mu (1-x^2) d(x) + w^2 * x = 0

% Inputs:
% w ->      Angular frequency parameter of the VdP equation (scalar).
% mu ->     Nonlinearity parameter of the VdP equation (array 1 x Nmu).
% time_f -> Time lenght of the simulation (scalar).
% dt ->     Time resolution (scalar).

% Outputs:
% signal: Solution of VdP equation (matrix Ns x Nmu), where Ns = length(0:+dt:time_f).

%Reference:
%Velarde O, Urdapilleta E, Mato G, and Dellavale D (2019), Bifurcation
%structure determines different phase-amplitude coupling patterns in the
%activity of biologically plausible neural networks, NeuroImage, In Press,
%(DOI: ...)

function signal = function_vanDerPolSolver_v0(mu, w, time_f, dt)
    
    timeVector = 0:+dt:time_f;
    Nmu = length(mu);
       
    %Memory pre-allocation to speed up the loop.
    signal = zeros(length(timeVector),Nmu);
    for ii = 1:+1:Nmu,
        [t,z] = ode45(@(t,z) function_vanDerPolDiffEq_v0(t,z,mu(ii),w),timeVector,[2;1]);
        signal(:,ii) = z(:,1);
    end
    
end

