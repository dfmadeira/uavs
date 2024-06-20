clc
clear all

zI = [0;0;1];
C = eye(12);
D = zeros(4);

% define constants
g = 9.81;  % Gravitational constant on Earth (m/s^2)
rho = 0.01225; %air density
m = 0.027; % Mass of the Drone
Cp = 1; %drag coefficient
l = 0.056650;%Drone's length (m)
w = 0.056650;%Drone's width (m)
h = 0.029;%Drone's height (m)
a = 0.025; %Drone's arm length(m)
Al=h*l;
Az=l*w;

%drag
Area = [l*h;w*h;l*w];
Cd = 1;
beta = diag([0.5*Cp*Al*rho, 0.5*Cp*Al*rho, 0.5*Cp*Az*rho]);

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

%Thrust
Teq=m*g; %for OP1

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
%fp = Mu * [T1;T2;T3;T4];

%np = [cross(T1_pos,f1) cross(T2_pos,f2) cross(T3_pos,f3) cross(T4_pos,f4)];
Mn = [n1+skew(T1_pos)*Mu(:,1) n2+skew(T2_pos)*Mu(:,2) n3+skew(T3_pos)*Mu(:,3) n4+skew(T4_pos)*Mu(:,4)];

%1.4 - 4x4 matrix that maps the four individual forces of the rotors

Mf = [Mu(3,:);Mn];

% simulation parameters
nx = 12;
ny = 4;
Dt = 0.01;
t = 0:Dt:10;
%u_L = [0.2*m_t*g_M;0.01;0.01;0.01]*(t>=0);


%introducing variables
syms x y z vx vy vz phi theta psi w1 w2 w3 T
p = [x;y;z];
v = [vx;vy;vz];
lbd = [phi;theta;psi];
omega = [w1;w2;w3];
R = Euler2R(lbd);

np = Mn * [T/4;T/4;T/4;T/4];
u = [T;phi;theta;psi];

%Drag = beta*v.^2

% expression

p_dot = R*v;
v_dot = -skew(omega)*v - g*R'*zI - R*beta*R'*v + 1/m*T*zI;
lambda_dot = Euler2Q(lbd)*omega;
omega_dot = -inv(J)*skew(omega)*J*omega - R*beta*R'*omega + inv(J)*np;

%first row
df1 = diff(p_dot,x);
df2 = diff(p_dot,y);
df3 = diff(p_dot,z);
df4 = diff(p_dot,vx);
df5 = diff(p_dot,vy);
df6 = diff(p_dot,vz);

%2nd row

df7 = diff(v_dot,x);
df8 = diff(v_dot,y);
df9 = diff(v_dot,z);
df10 = diff(v_dot,vx);
df11 = diff(v_dot,vy);
df12 = diff(v_dot,vz);

A = [df1 df2 df3 df4 df5 df6
    df7 df8 df9 df10 df11 df12];

B = [zeros(3,3);eye(3)];

C = eye(6); D=zeros(6,3);

%% LQR, linearization around hover T=0

%q0 = 1;
%Q = q0*eye(6);
Q = [1 0 0 0 0 0 %penalize px
     0 1 0 0 0 0 %penalize py
     0 0 1 0 0 0 %penalize pz
     0 0 0 1 0 0 %penalize vx
     0 0 0 0 1 0 %penalize vy
     0 0 0 0 0 1]; %penalize vz 
R = 1; %penalize control action

vars = [vx vy vz phi theta psi w1 w2 w3 T];
vals = zeros(1, length(vars));

A_hover = double(subs(A,vars,vals));
[K, P] = lqr(A_hover,B,Q,R);

% Closed-loop system dynamics
Ac = A_hover - B*K;
Bc = B;
Cc = C;
Dc = D;

% Create state-space system for closed-loop
sys_cl = ss(Ac, Bc, Cc, Dc);

% Define initial conditions and time vector
%x0 = zeros(6,1); % Example initial state
x0 = [4 5 3 0 0 0];
t = 0:0.01:100; % Time vector from 0 to 10 seconds with 0.01s time step

% Define the input signal (e.g., step input)
%u = [0;0;m*g]*ones(size(t)); % No external input for closed-loop simulation

% Simulate the system response
%y = lsim(sys_cl, u, t, x0);
[y,t,x] = initial(sys_cl, x0, t);

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');

% Plot the output trajectories
figure(1);
subplot(3,2,1)
plot(t, y(:,1));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$p_x$','Interpreter', 'latex');

subplot(3,2,3)
plot(t, y(:,2));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$p_y$', 'Interpreter', 'latex');

subplot(3,2,5)
plot(t, y(:,3));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$p_z$','Interpreter', 'latex');

subplot(3,2,2)
plot(t, y(:,4));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$v_x$', 'Interpreter', 'latex');

subplot(3,2,4)
plot(t, y(:,5));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$v_y$','Interpreter', 'latex');

subplot(3,2,6)
plot(t, y(:,6));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$v_z$','Interpreter', 'latex');

%% non linear 

% A e B devem ser as mesmas


%time parameters
Dt = 0.1;
t = 0:Dt:10;
Nsim = length(t);
%sim parameters
x0 = [0 0 3 0 0 0]';
x = [x0, zeros(6,Nsim-1)];
xe = zeros(6,Nsim);
y = zeros(6,Nsim);
u = zeros(3,Nsim);
ue = [0;0;m*g]*ones(size(t));


vals = [0     0     0     0     0     0     0     0     0     m*g+0]
21
for k = 1:Nsim
    
    if length(x(1,:)) ~= Nsim
        break
    end

    %define xdot
    p_dot = R*v;
    v_dot = -skew(omega)*v - g*R'*zI - R*beta*R'*v + 1/m*T*zI; %acho que este T está em sym??? nao percebo e esta tudo a zero
    p_dot = double(subs(p_dot,vars,vals));
    v_dot = double(subs(v_dot,vars,vals));
    x_dot = [p_dot;v_dot];


    % define control
    u(:, k) = K * (x(:, k) - xe(:, k)) + ue(:, k);  % correct indexing for u
    

    % integrate state
    if k < Nsim  % ensure we do not access out of bounds
        x(:, k+1) = x(:, k) + Dt * x_dot;
    end
    
    % current output
    y(:, k) = C * x(:, k);
end

figure(2);
subplot(3,2,1)
plot(t, y(1,:));
title('Output Trajectories', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Outputs', 'Interpreter', 'latex');
legend('$p_x$', 'Interpreter', 'latex');

subplot(3,2,3)
plot(t, y(2,:));
title('Output Trajectories', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Outputs', 'Interpreter', 'latex');
legend('$p_y$', 'Interpreter', 'latex');

subplot(3,2,5)
plot(t, y(3,:));
title('Output Trajectories', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Outputs', 'Interpreter', 'latex');
legend('$p_z$', 'Interpreter', 'latex');

subplot(3,2,2)
plot(t, y(4,:));
title('Output Trajectories', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Outputs', 'Interpreter', 'latex');
legend('$v_x$', 'Interpreter', 'latex');

subplot(3,2,4)
plot(t, y(5,:));
title('Output Trajectories', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Outputs', 'Interpreter', 'latex');
legend('$v_y$', 'Interpreter', 'latex');

subplot(3,2,6)
plot(t, y(6,:));
title('Output Trajectories', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Outputs', 'Interpreter', 'latex');
legend('$v_z$', 'Interpreter', 'latex');


%% define errors, os ganhos K, Q e R foram definidos anteriormente

x_ref = [0 0 0 0 0 0];
xref0 = [4 5 3 0 0 0]; % condições iniciais do sistema de referência
sys_ref = ss(A_hover, B, C, D)

t = 0:0.01:100; % Time vector from 0 to 10 seconds with 0.01s time step

% Define the input signal (e.g., step input)
u = [0;0;0]*ones(size(t)); % No external input for closed-loop simulation

% Simulate the system response
[y_ref, t, xr] = lsim(sys_cl, u, t, xref0);

x_error = x - xr
y_error = C*x_error'

figure(4)
%plot(t, y_error(3,:),'g',t,y_error(2,:),'b',t,y_error(3,:),'r');

% Plot the output trajectories
subplot(3,2,1)
plot(t, y_error(1,:));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$p_x$','Interpreter', 'latex');

subplot(3,2,3)
plot(t, y_error(2,:));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$p_y$', 'Interpreter', 'latex');

subplot(3,2,5)
plot(t, y_error(3,:));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$p_z$', 'Interpreter', 'latex');

subplot(3,2,2)
plot(t, y_error(4,:));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$v_x$','Interpreter', 'latex');

subplot(3,2,4)
plot(t, y_error(5,:));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$v_y$','Interpreter', 'latex');

subplot(3,2,6)
plot(t, y_error(6,:));
title('Output Trajectories');
xlabel('Time (s)');
ylabel('Outputs');
legend('$v_z$', 'Interpreter', 'latex');


%% Part 2: Nonlinear Control and Trials

%2.1 - Design a nonlinear controller with actuation in body accelerations, 
% asssuming a zero yaw, that is able to achieve asymptotic stability in the 
% Lyapunov sense

drag = beta; %ver dados iniciais
% zB = [theta;-phi;1];
% p_dot = v;
% v_dot = -g*zI-drag*v+T/m*zB;
% vd = [0;0;2];
% 
% for i=1:Nsim
%     f = T*zB;
%     vd_dot = f/m;
%     ev = v - vd;
% 
%     ep_dot = ev;
%     ev_dot = -g*zI - drag*v + 1/m*f - vd_dot;
% end

% initialize variables for all drones:
Tend = 70;
dTo = 0.1;
dTi = 0.05;
Nsim = round(Tend/dTi)+1;
p_ref_static = [0.5;0.5;1];

t = 0:dTi:Tend;
nt = length(t);
nx = 6;
nu = 4;

p0 = [0;0;0];
v0 = [0;0;0];
x = zeros(nx,Nsim);
T = zeros(1,Nsim);
x(:,1) = [p0;v0];

% Gains for initernal controller
kp=0.2; kv=0.23;

% step reference
% p_ref = [zeros(3,20/dTi+1), p_ref_static*ones(1,(Tend-20)/dTi)];
% v_ref = zeros(3,nt);
% a_ref = zeros(3,nt);

% circle reference
Rad = 0.5;      % radius of circle
omn = 2*pi/20;  % rotation frequency

p_ref = [Rad*cos(omn*t);Rad*sin(omn*t);ones(size(t))];
v_ref = [-Rad*omn*sin(omn*t);Rad*omn*cos(omn*t);0*ones(size(t))];
a_ref = [-Rad*omn^2*cos(omn*t);-Rad*omn^2*sin(omn*t);0*ones(size(t))];


% main time loop for simulation
for k = 1:Nsim

    % get state vector and plot it
    p = x(1:3,k);
    v = x(4:6,k);
    p_d = p_ref(:,k);
    v_d = v_ref(:,k);
    a_d = a_ref(:,k);

    % outer-loop controller
    e_p = p - p_d;
    e_v = v - v_d;

    % Mellinger Controller (up to attitude commands)
    f_des = -kp*e_p - kv*e_v + m*g*zI + m*a_d;

    % compute desired rotation matrix
    zB_des(:,k) = f_des/norm(f_des);

    % compute thrust
    %T(:,k) = f_des'*zB_des(:,k);
    T(:,k) = max(min(0.35,f_des'*zB_des(:,k)),0);
    
    
    % nonlinear drone model
    dot_p = v;
    dot_v = -g*zI - drag*v + T(:,k)/m*zB_des(:,k);

    % discretization 
    pp = p + dTi*dot_p;
    vp = v + dTi*dot_v;
    if k~=Nsim
        x(:,k+1) = [pp;vp];
    end
end

p = x(1:3,:);

% drone_show_data;
% show results plot
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstlighterblue  = [50,220,240]/255;
sstlightestblue = [60,230,255]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstlightergreen = [160,255,225]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;

dcolors = { sstgreen, sstblue, sstlightblue, sstlighterblue, sstlightestblue, sstlightgreen, sstlightergreen, sstlightgray };

%angle plots
figure(101);

subplot(311);
plot(t,T,'Color',dcolors{1});
hold on;
grid on;
ylabel('$$T(t)$$ [N]');
title('Control variables');

subplot(312);
plot(t,-zB_des(2,:),'Color',dcolors{1}); %plots phi
hold on;
grid on;
ylabel('$$\phi(t)$$ [rad]');

subplot(313);
plot(t,zB_des(1,:),'Color',dcolors{1}); %plots theta
hold on;
grid on;
ylabel('$$\theta(t)$$ [rad]');
xlabel('$$t$$ [s]');

%position plots
figure(102);

subplot(311);
plot(t,p_ref(1,:),'Color',sstgray);
hold on;
plot(t,x(1,:),'Color',dcolors{1});
hold off;
grid on;
ylabel('$$x(t)$$ [m]');
title('Drone position and reference');

subplot(312);
plot(t,p_ref(2,:),'Color',sstgray);
hold on;
plot(t,x(2,:),'Color',dcolors{1});
hold off;
grid on;
ylabel('$$y(t)$$ [m]');

subplot(313);
plot(t,p_ref(3,:),'Color',sstgray);
hold on;
plot(t,x(3,:),'Color',dcolors{1});
hold off;
grid on;
xlabel('$$t$$ [s]');
ylabel('$$z(t)$$ [m]');

%drone_animate(p,p_ref,[-zB_des(2,:);zB_des(1,:);zeros(1,length(zB_des))],t,dcolors);
lbd = [-zB_des(2,:);zB_des(1,:);zeros(1,length(zB_des))];

sstgray = [70,70,70]/255;
nt = length(t);

figure(104);

hini = plot3(p(1,1),p(2,1),p(3,1),'o','Color',dcolors{1},'MarkerSize',2);
hold on;
href = plot3(p_ref(1,:),p_ref(2,:),p_ref(3,:),'--','Color',sstgray);
hp = plot3(p(1,1:2),p(2,1:2),p(3,1:2),'-','Color',dcolors{1});
hd = drone_plot(p(1:3,1),lbd(:,1),[],dcolors{1});

hold off;
grid on;
axis equal;
axis([-1.2 1.2 -1.2 1.2 0 3]);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
legend('start','end','trajectory');
for k = 2:2:nt
    set(hp,'XData',p(1,1:k),'YData',p(2,1:k),'ZData',p(3,1:k));
    drone_plot(p(1:3,k),lbd(:,k),hd);

    axis equal;
    drawnow;
    %pause(dt/10);
end



