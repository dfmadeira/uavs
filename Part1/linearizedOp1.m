clear,clc
syms px py pz fi theta yaw vx vy vz wx wy wz T g_earth rho Cd  m cte Jx Jy Jz npx npy npz

load('A_matrix.mat')
load('B_matrix.mat')

vars = [px py pz fi theta yaw vx vy vz wx wy wz T g_earth rho Cd  m cte Jx Jy Jz npx npy npz]

zI = [0;0;1];
C = eye(12);
D = zeros(12,4);

%% equilibrium conditions (OP1)
px = 0;
py = 0;
pz = 1.5;
fi = 0;
theta = 0;
yaw = 0;
vx = 0;
vy = 0;
vz = 0;
wx = 0;
wy = 0;
wz = 0;
g_earth = 9.81;
rho = 1.225;
Cd = 0;
m = 0.027;
T = m*g_earth*0
cte = 0;
Jx = 9.113e-6;
Jy = 9.113e-6;
Jz = 1.444e-05;
npx = 0;
npy = 0;
npz = 0;
% 1. Define OP1 parameters
valsOP1 = [px py pz fi theta yaw vx vy vz wx wy wz T g_earth rho Cd  m cte Jx Jy Jz npx npy npz]

AOP1 = double(subs(A,vars,valsOP1))
BOP1 = double(subs(B,vars,valsOP1))
%hover(AOP1, BOP1, C, D, 0.027*9.81,1)
nx = 12;
ny = 4;
Dt = 0.01;
t = 0:Dt:100;
u_eq = [0.027*9.81*0;0;0;0]*ones(size(t)); %divided by 4 to consider all motors;

% simulate nonlinear system
Nsim = length(t);
xe = zeros(nx,Nsim);
xe(:,1) = [0;0;0;0;0;0;0;0;0;0;0;0];
x_eq = [0;0;1.5;0;0;0;0;0;0;0;0;0]*ones(size(t));
%y = zeros(ny,Nsim);

sys=ss(AOP1,BOP1,C,D);
y=lsim(sys, u_eq, t, xe(:,1));%,x_eq);

figure(2)
subplot(1,2,1)
plot(t, y(:,1),t, y(:,2),t, y(:,3));
legend(["x" "y" "z"])
xlabel('Time')
ylabel('Position')
title('Response of State Variables')
grid on
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times
subplot(1,2,2)
plot(t, y(:,4),t, y(:,5),t, y(:,6));
legend(["vx" "vy" "vz"])
xlabel('Time')
ylabel('Position')
title('Response of State Variables')
grid on
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times

% controlability, observability and stability

[Vj,Jor] = schur(AOP1),
[V,DL,W] = eig(AOP1),
if any(real(diag(DL)) >=0 ), disp('Linearized system OP1 is not stable.'); end
n1_unstable_modes = rank(ctrb(AOP1,BOP1))-12;
if n1_unstable_modes > 0, disp('Linearized system OP1 is not controlable.'); end
n1_unobservable_modes = rank(obsv(AOP1,C))-12;
if n1_unobservable_modes > 0, disp('Linearized system OP1 is not observable.'); end


%% 2. Define OP2 parameters
theta=5;
%vars = [px py pz fi theta yaw vx vy vz wx wy wz T g_earth rho Cd  m cte Jx Jy Jz npx npy npz]
valsOP2 = [0 0 0 0 theta*pi/180 0 3.387 0 0 0 0 0 0.2639-0.027*9.81 9.81 1.225 1 0.027 0 9.113e-6 9.113e-6 1.444e-05 0 0 0]

AOP2 = double(subs(A,vars,valsOP2)) 
BOP2 = double(subs(B,vars,valsOP2))
%u_eq = [T, npx, npy, npz]
u_eq = [0.027*9.81-0.2639*cosd(5);0;0;0]
hover(AOP2, BOP2, C, D, u_eq,2);

[Vj,Jor] = schur(AOP2),
[V,DL,W] = eig(AOP2),
if any(real(diag(DL)) >=0 ), disp('Linearized system OP1 is not stable.'); end
n1_unstable_modes = rank(ctrb(AOP2,BOP2))-12;
if n1_unstable_modes > 0, disp('Linearized system OP1 is not controlable.'); end
n1_unobservable_modes = rank(obsv(AOP2,C))-12;
if n1_unobservable_modes > 0, disp('Linearized system OP1 is not observable.'); end
