%non linear model (1.4)

%pkg load symbolic

%syms phi theta psi real
%syms omega_x omega_y omega_z real
%syms velocity_x velocity_y velocity_z real

clc
clear all



%Lander dynamic system modeling
zI = [0;0;1];
C = [   eye(3)    , zeros(3)   , zeros(3) , zeros(3)
        zeros(1,3), zeros(1,3) , zI' , zeros(1,3)  ];
D = zeros(4);

%define constants
g = 9.81;  % Gravitational constant on Earth (m/s^2)
rho = 0.01225; %air density
% Az=0.03*0.02
% Ay=0.01*0.01
% bz = 0.5*Cp*rho*Az;
% by = 0.5*Cp*rho*Ay;      % Drag coefficient
m = 0.027;      % Mass of the Drone
Cp = 1; %drg coefficient

l = 0.056650;%Drone's length (m)
w = 0.056650;%Drone's width (m)
h = 0.029;%Drone's height (m)
a = 0.025; %Drone's arm length(m)
Al=h*l
Az=l*w

CM = [0;0;0];%Center of mass [x;y;z]
%Thruster Positions relative to CM
T1_pos = [l/2;-w/2;h/2];%distances according to the drone's dimensions 
T2_pos = [-l/2;-w/2;h/2];%distances according to the drone's dimensions 
T3_pos = [-l/2;w/2;h/2];%distances according to the drone's dimensions 
T4_pos = [l/2;w/2;h/2];%distances according to the drone's dimensions 

% Rigid-body nonlinear model
f_g = m * g;% Gravity force
J = (m/12) * diag([h^2+l^2,h^2+l^2,l^2+w^2]); %Moment of Inertia Matrix
% v_dot = -skew(omega)*v + f_g + (1/m_t)*(f_a + f_p);
% omega_dot = -inv(J)*skew(omega)*J + inv(J)*(n_a + n_p);
Area = [l*h;w*h;l*w];
Cd = 1;

Teq=m*g %for OP1

T1 = Teq/4;
T2 = Teq/4;
T3 = Teq/4;
T4 = Teq/4;

n1 = [0;0;-1];
n2 = [0;0;1];
n3 = [0;0;-1];
n4 = [0;0;1];

%1.3 - Net force and moment vectors
f1 = [0;0;1];
f2 = [0;0;1];
f3 = [0;0;1];
f4 = [0;0;1];

Mu = [f1, f2, f3, f4];
fp = Mu * [T1;T2;T3;T4];

%np = [cross(T1_pos,f1) cross(T2_pos,f2) cross(T3_pos,f3) cross(T4_pos,f4)];
Mn = [n1+skew(T1_pos)*Mu(:,1) n2+skew(T2_pos)*Mu(:,2) n3+skew(T3_pos)*Mu(:,3) n4+skew(T4_pos)*Mu(:,4)];
np = Mn * [T1;T2;T3;T4];

%1.4 - 4x4 matrix that maps the four individual forces of the rotors

Mf = [Mu(3,:);Mn];

% simulation parameters
nx = 12;
ny = 4;
Dt = 0.01;
t = 0:Dt:10;
%u_L = [0.2*m_t*g_M;0.01;0.01;0.01]*(t>=0);
Te = m*g+0.01; % equilibirum thrust (divided by the number of motors)
u_NL = [Te;0;0;0]*ones(size(t)); %+ u_L;

% simulate nonlinear system
Nsim = length(t);
x = zeros(nx,Nsim);
y = zeros(ny,Nsim);

for k = 1:Nsim
    % prepare variables:
    p   = x(1:3,k);
    v   = x(4:6,k);
    lambda = x(7:9,k);
    omega  = x(10:12,k);
    R = Euler2R(lambda);
    T = u_NL(1,k);
    np = u_NL(2:4,k);

    %1.3 - Mars gravity and air drag force vectors
    beta = 0.5*rho*Area*Cd;
    drag = beta .* v.^2.*eye(3);
    
    % compute state derivative:
    p_dot = R*v;
    v_dot = -skew(omega)*v - g*R'*zI + 1/m*((-R*drag*R'*v + T).*zI);
    lambda_dot = Euler2Q(lambda)*omega;
    omega_dot = -inv(J)*skew(omega)*J*omega - R*drag*R'*omega + inv(J)*np;
    x_dot = [p_dot;v_dot;lambda_dot;omega_dot];

    % integrate state
    x(:,k+1) = x(:,k) + Dt*x_dot;

    % compute current output:
    y(:,k) = C*x(:,k);
end

figure(20000)
plot(t, y(3,:))
title('Non-Linear Simulation')
xlabel('Time')
ylabel('Height')
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times
%% Linearization around the given equilibrium points
%(considering phi, theta, psi = 0 in equilibrium)
R_lambda_diff_eq = [0 1 0;-1 0 0;0 0 0];
%OP1 -> zero vertical velocity => T1=T2=T3=T4=m*g
ve = 0;
% A_OP2_a = [zeros(3,3), eye(3), zeros(3,3), zeros(3,3)
%            zeros(3,3), zeros(3,3), g_M.*[0 -1 0;1 0 0;0 0 0], zeros(3,3)
%            zeros(3,3), zeros(3,3), zeros(3,3), eye(3)
%            zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3)];


dfv_dv=(-2/m)*[0.00164285*1.225, 0, 0
               0, 0.00164285*1.225, 0
               0, 0, 0.0032092225*1.225]

A_OP1 = [zeros(3,3), eye(3), ve.*R_lambda_diff_eq, zeros(3,3)
           zeros(3,3), 0*dfv_dv, g*R_lambda_diff_eq, ve.*R_lambda_diff_eq
           zeros(3,3), zeros(3,3), zeros(3,3), eye(3)
           zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3)];


%tirei isto 1/m*2*beta.*ve.*diag(zI) onde está 
%Te = m*g;
B_OP1 = [zeros(3,4);(1/m)*Mu;zeros(3,4);Mn];

B = [
    zeros(3, 4);
    (1/m) * [
        [0; 0; 1], [0; 0; 1], [0; 0; 1], [0; 0; 1]
    ];
    zeros(3, 4);
    Mn%J(:,1).*ones(3,4)%Mn
];

C = eye(12);
D = zeros(12,4);
Thrust = m*g*0.5/4;
u_eq = [Thrust;0;0;0]*ones(size(t)); 
x_init = [0;0;1.5;0;0;0;0;0;0;0;0;0];
sys=ss(A_OP1,B,C,D);
y=lsim(sys,u_eq, t,x_init);%,x_eq);


figure(2)

plot(t, y(:,1),t, y(:,2),t, y(:,3));
legend(["x" "y" "z"])
xlabel('Time')
ylabel('Position')
title('Response of State Variables')
grid on
figure(3)
plot(t, y(:,4),t, y(:,5),t, y(:,6));
%axis([0 5])
legend(["vx" "vy" "vz"])
xlabel('Time')
ylabel('Position')
title('Response of State Variables')
grid on