% Define symbolic variables
syms x z T

% Given constants
g = 9.81;
theta = 5; % Convert degrees to radians
rho = 1.225;
m = 0.027;

% Define the equations
eq1 = g*sind(theta) - (1/m)*0.0016428*rho*x^2 == 0;
eq2 = -g*cosd(theta) + (1/m)*T == 0;

% Solve the system of equations
[sol_x, sol_T] = solve(eq1, eq2, x, T);

% Display the solutions
sol_x = vpa(simplify(sol_x))
sol_T = vpa(simplify(sol_T))
