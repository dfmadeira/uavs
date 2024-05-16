clear, clc
%% L2Data1
% z está a 1 porque mede 1*g
df1 = csvread("L2Data1.csv",1,0);
stats('L2Data1.csv', 14);

cvargyx = std(deg2rad(df1(:,11)))^2;
cvargyy = std(deg2rad(df1(:,12)))^2;
cvarphi = std(deg2rad(df1(:,8)))^2;
cvartheta = std(deg2rad(df1(:,9)))^2;

%% L2Data2
stats('L2Data2.csv', 14);

%% Estimate inclination 
% clear, clc

syms phi theta psi axm aym azm;

% define constants
g = 9.8;

lbd = [phi; theta; psi];
R = Euler2R(lbd);
e3 = [0; 0; 1];

R = subs(R, psi, 0);

Rapr = [1 phi*theta theta;
    0 1 -phi;
    -theta phi 1];
%Re3 = [theta; -phi; 1];
Re3 = Rapr * e3;

am = [axm,; aym; azm];
a = [0; 0; 0]; %impossivel estimar

%% 1.4

syms amx nx bx amy ny by ntheta nphi bz nz;

thetam = (amx-nx-bx)/g - ntheta;

phim = (amy-ny-by)/-g - nphi; 

varstheta = [amx, nx, bx, ntheta];
varsphi= [amy, ny, by, nphi];

%estimation of L2Data3

df3 = csvread("L2Data4.csv",1,0);
% df3(1,10) %first value of x
% df3(1,11)
% df3(1,12)
% df3(2000,8)
%valstheta = [df3(1,10), 0, 0, 0];

%round(subs(thetam, varstheta, valstheta),10)

Dt=1;
t = 0:Dt:1000;
Nsim = length(t);

ytheta = zeros(Nsim,1);
ythetam = zeros(Nsim,1);
meantheta = ythetam;

yphi = zeros(Nsim,1);
yphi = zeros(Nsim,1);
meanphi = yphi;


for k = 1:Nsim

    valstheta = [df3(k,10), 0, 0, 1.1327];
    ytheta(k,:) = round(subs(thetam, varstheta, valstheta),10);
    ythetam(k,:) = df3(k,9);
    meantheta = ytheta(k,:) - ythetam(k,:);

    valsphi = [df3(k,11), 0, 0, 3.4613];
    yphi(k,:) = round(subs(phim, varsphi, valsphi),10);
    yphim(k,:) = df3(k,8);
    meanphi = yphi(k,:) - yphim(k,:);

    %round(k*100/Nsim);
end

figure(110)
plot(t, ytheta, t, ythetam)
title('Estimation of theta')
legend('Estimation','Measurement')
set(gca, 'FontName', 'Times New Roman')

figure(111)
plot(t, yphi, t, yphim)
title('Estimation of phi')
legend('Estimation','Measurement')
set(gca, 'FontName', 'Times New Roman')


%% Kalman Filter Design

%2.1 - Design a linear Kalman filter to estimate both pitch and roll based
%on previous measurements

df4 = csvread("L2Data6.csv",1,0);

num_states = 6;

A = [0     0     1     0     0     0
     0     0     0     1     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0];

B = zeros(6,0); G = eye(6);

C = [0     0     1     0     1     0
     0     0     0     1     0     1
     1     0     0     0     0     0
     0     1     0     0     0     0];

D = zeros(4,0); H = zeros(4,6);

Q = 0.01*eye(6); %diag([1,1,1,1,1,1]); %
R = 500*diag([cvarphi,cvartheta,cvargyx,cvargyy]);

%2.2 - Filter observability and stability analysis
if rank(ctrb(A,Q)) >= num_states; disp('The model is controllable');
elseif rank(ctrb(A,Q)) < num_states; disp('The model is not controllable');end

if rank(obsv(A,C)) >= num_states; disp('The model is observable');
elseif rank(obsv(A,C)) < num_states; disp('The model is not observable');end

%2.3 - Filter implementation and comparison with the attitude data
%Ts = 1;
%sys = ss(A,[B G],C,[D H]);
%[K_filter, L, P] = kalman(sys,Q,R); % L -> Kalman Gain ; P -> covariance matrix


Ts = 1;
A = expm(A*Ts); % discrete state matrix
G = [1;0;0;0;0;0];
H = C;
ny = size(H,1);
nx = size(H,2);

% kalman gains
% Q = diag([1,1,100,100,1,1]);
% R = diag([cvarphi,cvartheta,cvaracx,cvaracy]);

% initial values
P0 = eye(nx);
x0 = zeros(nx,1);

% simulate system and estimation
Tsim = length(df4(:,1));
t = 0:Ts:Tsim-1;
Nsim = length(t);
u_th = 0;
x = zeros(nx,Nsim);
%y = zeros(ny,Nsim);
y = deg2rad([df4(:,8),df4(:,9),df4(:,11),df4(:,12)])';
xe = zeros(nx,Nsim);
Pe = zeros(nx,nx,Nsim);
%x(:,1) = x0;
xe(:,1) = x0;
Pe(:,:,1) = P0;
for k = 1:Nsim

    %predict next estimate:
    Gu = G*u_th;
    [xem,Pem] = kalman_predict(xe(:,k),Pe(:,:,k),A,Gu,Q);

    % get measurements (simulated)
    % v_noise = 0*sqrtm(R)*randn(ny,1);
    % y(:,k) = H*x(:,k) + v_noise;

    % update estimate with measurement info
    [xe(:,k),Pe(:,:,k),K] = kalman_update(xem,Pem,y(:,k),H,R);

    % simulate system dynamics and process noise
    w_noise = 0*sqrtm(Q)*randn(nx,1); % process noise not added
    xp = A*x(:,k) + Gu + w_noise; % process noise not added
    if k < Nsim % update next state until last instant
        x(:,k+1) = xp;
    end

    % xe(:,k) = A*xe(:,k) + Gu + w_noise;
    % Pe(:,:,k) = A*Pe(:,:,k)*A' + Q;
    % 
    % K = Pe(:,:,k)*C'*(C*Pe(:,:,k)*C' + R)^-1;
    % xe(:,k+1) = xe(:,k) + K*(y(k)-C*xe(:,k));
    % Pe(:,:,k) = (eye(num_states)-K*C)*Pe(:,:,k);
    % 
    % xe(:,k+1) = A*xe(:,k)+Gu;
    % y(:,k+1) = C*xe(:,k+1);
end



% Show results
figure(1001);
subplot(3,1,1)
plot(df4(:,1), df4(:,11), df4(:,1), xe(1,:)*180/pi, 'g')% df4(:,1), x(1,:)*180/pi, 'r') %Gyro x
legend('data','estimated');
title('Gyro X')
subplot(3,1,2)
plot(df4(:,1), df4(:,12), df4(:,1), xe(2,:)*180/pi, 'g')% df4(:,1), x(2,:)*180/pi, 'r') %Gyro y
title('Gyro Y')
subplot(3,1,3)
plot(df4(:,1), df4(:,13)) %Gyro z
title('Gyro Z')




%% Extended Kalman Filter Design

%3.1 - Design an extended Kalman filter to estimate both pitch and roll based
%on the original accelerometers and rate gyro measurements

syms phi theta p q bp bq ax ay az amx amz amz nz bz bx nx;



%define matrices
x = [phi; theta; p; q; bp; bq];
bac = [bx; by; bz];
nac = [nx;ny;nz];
a=[ax,ay,az];

R = Euler2R(lbd);

am = R'*(a - g*e3) + bac + nac;

%first row elements
df1dphi = diff(am(1),phi);
df1dtheta = diff(am(1),theta);
df1dp = 0;
df1dq = 0;
df1bphi = 0;
df1btheta = 0;

%2nd row elements
df2dphi = diff(am(2),phi);
df2dtheta = diff(am(2),theta);
df2dp = 0;
df2dq = 0;
df2bphi = 0;
df2btheta = 0;

%3rd row elements
df3dphi = diff(am(3),phi);
df3dtheta = diff(am(3),theta);
df3dp = 0;
df3dq = 0;
df3bphi = 0;
df3btheta = 0;


% Aekf = [df1dphi, df1dtheta, df1dp, df1dq, df1bphi, df1btheta
%         df2dphi, df2dtheta, df2dp, df2dq, df2bphi, df2btheta
%         df3dphi, df3dtheta, df3dp, df3dq, df3bphi, df3btheta];

Aekf = A;

Bekf = B; Bw = eye(4);

phim = acos((amz-nz-bz)/(g*sin(theta)) - nphi);%asin((amy-ny-by)/-g - nphi);
thetam = acos((amx-nx-bx)/(g*sin(phi)) - ntheta);

Cekf = [0     0     1     0     1     0
        0     0     0     1     0     1
        diff(phim,phi)     diff(phim,theta)     diff(phim,p)     diff(phim,q)     diff(phim,bp)     diff(phim,bq)
        diff(thetam,phi)     diff(thetam,theta)     diff(thetam,p)     diff(thetam,q)     diff(thetam,bp)     diff(thetam,bq)];

Dekf = D; Dv = eye(4);

%3.2 - Filter observability by linearizing the system

% Cekf_lin = [0     0     1     0     1     0
%             0     0     0     1     0     1
%             0     0     0     0     0     0
%             0     0     0     0     0     0];

vars = [phi theta bx nx amx amy amz nz bz ntheta nphi];
vals = [df4(1,11) df4(1,13) 0 0 df4(1,14) df4(1,15) df4(1,16) 0 0 0 0];

Cekf_lin = double(subs(Cekf,vars,vals));

if rank(ctrb(Aekf,Q)) >= num_states; disp('The model is controllable');
elseif rank(ctrb(Aekf,Q)) < num_states; disp('The model is not controllable');end

if rank(obsv(Aekf,Cekf_lin)) >= num_states; disp('The model is observable');
elseif rank(obsv(Aekf,Cekf_lin)) < num_states; disp('The model is not observable');end

%3.3 - Filter implementation and comparison with the linear filter

Aekf = expm(Aekf*Ts); % discrete state matrix
Q = 0.0005*eye(6);
R = diag([cvarphi,cvartheta,cvargyx,cvargyy]);

x_ext = zeros(num_states,Nsim);
%y = zeros(ny,Nsim);
y = deg2rad([df4(:,8),df4(:,9),df4(:,11),df4(:,12)])';
xe_ext = zeros(num_states,Nsim);
Pe_ext = zeros(num_states,num_states,Nsim);
%x(:,1) = x0;
xe_ext(:,1) = x0;
Pe_ext(:,:,1) = P0;


for k = 1:Nsim
    vals = [df4(k,11) df4(k,13) 0 0 df4(k,14) df4(k,15) df4(k,16) 0 0 0 0];
    Cekf = double(subs(Cekf,vars,vals));

    %predict next estimate:
    Gu = G*u_th;
    [xem,Pem] = kalman_predict(xe_ext(:,k),Pe_ext(:,:,k),Aekf,Gu,Q);

    % update estimate with measurement info
    [xe_ext(:,k),Pe_ext(:,:,k),K] = kalman_update(xem,Pem,y(:,k),Cekf,R);

    % simulate system dynamics and process noise
    w_noise = 0*sqrtm(Q)*randn(num_states,1); % process noise not added
    xp = Aekf*x_ext(:,k) + Gu + w_noise; % process noise not added
    if k < Nsim % update next state until last instant
        x_ext(:,k+1) = xp;
    end
    k
end

figure(1002);
subplot(3,1,1)
plot(df4(:,1), df4(:,11), df4(:,1), xe_ext(1,:)*180/pi, 'r')% df4(:,1), x(1,:)*180/pi, 'r') %Gyro x
legend('data','estimated');
title('Gyro X')
subplot(3,1,2)
plot(df4(:,1), df4(:,13), df4(:,1), xe_ext(2,:)*180/pi, 'r')% df4(:,1), x(2,:)*180/pi, 'r') %Gyro y
title('Gyro Y')

