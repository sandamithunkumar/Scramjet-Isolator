clear; close all; clc;

gamma_iso = 1.37; R_iso = 287;  cp_iso = 1063;
gamma_comb = 1.31; R_comb = 297; cp_comb = 1255;

Cf = 0.002;
D = 0.06;
fst = 0.029;
phi = 0.71;
h_PR = 1.2e8;
eta_c_tot = 0.8;
vartheta = 5.0;

M2 = 2.65;
T2 = 650;
P2 = 50e3;
A2 = pi*D^2/4;

Ht2 = 1.59e6;
Tt2 = Ht2/cp_iso;
Pt2 = P2*(Tt2/T2)^(gamma_iso/(gamma_iso-1));

x2 = 0.0; x3 = 0.2; x4 = 0.5;

% Property switches as functions of x
gamma_fun = @(x) (x <= x3).*gamma_iso + (x > x3).*gamma_comb;
cp_fun    = @(x) (x <= x3).*cp_iso    + (x > x3).*cp_comb;
% R_fun not used below, keep if needed:
R_fun     = @(x) (x <= x3).*R_iso     + (x > x3).*R_comb;

% Geometry and combustion scheduling
dAdx  = @(x, A) (x <= x3)*0 + (x > x3)*(A2/(x4 - x3));
X_fun = @(x) (x - x3) / (x4 - x3);
detadx= @(x) (eta_c_tot * vartheta / (x4-x3)) ./ (1 + (vartheta-1) * X_fun(x)).^2;

% ODE right-hand sides with x-dependent gamma, cp
dTtdx = @(x, Tt) (x <= x3)*0 + (x > x3)*((h_PR * fst * phi ./ cp_fun(x)) .* detadx(x));
dPtdx = @(x, Tt, Pt, M) Pt .* (-gamma_fun(x) .* M.^2) .* ( dTtdx(x, Tt) ./ Tt + (4 * Cf) / D );
dMdx  = @(x, Tt, M, A) (M .* (1 + ((gamma_fun(x) - 1)/2) .* M.^2) ./ (1 - M.^2)) .* ...
    ( ((1 + gamma_fun(x) .* M.^2) / 2) .* (dTtdx(x, Tt) ./ Tt) ...
    + (gamma_fun(x) .* M.^2) .* (2 * Cf / D) - dAdx(x, A) ./ A );
dPdx  = @(x, Tt, M, A, P) (P .* (gamma_fun(x) .* M.^2) ./ (1 - M.^2)) .* ...
    ( dAdx(x, A) ./ A - (1 + (gamma_fun(x) - 1) .* M.^2) .* (2 * Cf / D) ...
    - (1 + ((gamma_fun(x) - 1) .* M.^2) / 2) .* (dTtdx(x, Tt) ./ Tt) );
dTdx  = @(x, Tt, M, A, T) (T ./ (1 - M.^2)) .* ...
    ( (gamma_fun(x) - 1) .* M.^2 .* ( dAdx(x, A) ./ A - gamma_fun(x) .* M.^2 .* (2 * Cf / D) ) ...
    + (1 - gamma_fun(x) .* M.^2) .* (1 + (gamma_fun(x) - 1) .* M.^2 / 2) .* (dTtdx(x, Tt) ./ Tt) );

odes = @(x, y) [ ...
    dAdx(x, y(1));                   % dA/dx
    dTtdx(x, y(2));                  % dTt/dx
    dPtdx(x, y(2), y(3), y(4));      % dPt/dx
    dMdx(x, y(2), y(4), y(1));       % dM/dx
    dPdx(x, y(2), y(4), y(1), y(5)); % dP/dx
    dTdx(x, y(2), y(4), y(1), y(6))  % dT/dx
    ];

xspan = linspace(x2, x4, 300);
y0 = [A2; Tt2; Pt2; M2; P2; T2];

% Optional: stop at sonic condition
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
% opts = odeset(opts,'Events',@(x,y) sonicEvent(x,y));

[x, y] = ode45(odes, xspan, y0, opts);

A  = y(:, 1);
Tt = y(:, 2);
Pt = y(:, 3);
M  = y(:, 4);
P  = y(:, 5);
T  = y(:, 6);

A_A2 = A / A2;
Tt_Tt2 = Tt / Tt2;
Pt_Pt2 = Pt / Pt2;
P_P2 = P / P2;
T_T2 = T / T2;

figure;
plot(x, A_A2, '-r', 'LineWidth', 2); hold on;
plot(x, Tt_Tt2, '-b', 'LineWidth', 2);
plot(x, Pt_Pt2, '-g', 'LineWidth', 2);
plot(x, M, 'k-', 'LineWidth', 2);
xlabel('x(m)'); ylabel('A/A2, Pt/Pt2, Tt/Tt2, M');
legend('A(x)/A_2', 'Tt(x)/Tt_2', 'Pt(x)/Pt_2', 'M', 'Location', 'northwest');
title('Normalized Area, Total Temperature, Total Pressure, and Mach Number');
grid on;

xlim([0 0.5]);
ylim([0 5]);
xticks(0:0.1:0.5);
yticks(0:1:5);

figure;
plot(x, A_A2, '-r', 'LineWidth', 2); hold on;
plot(x, M, 'k-', 'LineWidth', 2);
plot(x, P_P2, '-g', 'LineWidth', 2);
plot(x, T_T2, '-b', 'LineWidth', 2);
xlabel('x(m)'); ylabel('A/A2, P/P2, T/T2, M');
legend('A(x)/A_2', 'M', 'P(x)/P_2', 'T(x)/T_2', 'Location', 'northwest');
title('Normalized Area, Mach, Static Pressure, and Static Temperature');
grid on;

xlim([0 0.5]);
ylim([0 5]);
xticks(0:0.1:0.5);
yticks(0:1:5);

function [value,isterminal,direction] = sonicEvent(x,y)
    M = y(4);
    value = 1 - M.^2;
    isterminal = 1;
    direction = -1;
end
