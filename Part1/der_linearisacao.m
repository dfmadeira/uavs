
syms px py pz fi theta yaw vx vy vz wx wy wz T g_earth rho Cd m cte Jx Jy Jz npx npy npz

l = 0.056650;%Drone's length (m)
w = 0.056650;%Drone's width (m)
h = 0.029;%Drone's height (m)
Al=h*l;
Az=l*w;

Area = diag([Al,Al,Az]);
v=[vx
    vy
    vz];

%modv = sqrt(vx^2+vy^2+vz^2)
%Drag_Body = 0.5*Cd*rho*Area*v*v'
%Drag_Body=diag([dragx,dragy,dragz]);
fa=-Cd*rho*Area*v;%*modv;

omg=[wx
    wy
    wz];

fp=[0
    0
    T];
zI=[0
    0
    1];


%pos1=[r1x; r1y; r1z];
%pos2=[r2x; r2y; r2z];
%pos3=[r3x; r3y; r3z];
%pos4=[r4x; r4y; r4z];


%np=moment_vector_func([T1;T2;T3;T4], zI*T1,zI*T2,zI*T3,zI*T4,pos1,pos2,pos3,pos4,cte)

np=[npx;npy;npz];

J=[Jx 0 0;
   0 Jy 0;
    0 0 Jz];
lbd=[fi
    theta
    yaw];

R = Euler2R(lbd);

f1 = R*v;

disp("Aqui começa as derivadas do f1")
disp("Em função da posição")
A1x=diff(f1,px);
A1y=diff(f1,py);
A1z=diff(f1,pz);
A1=[A1x, A1y, A1z];
disp("Em função da velocidade")
A2x=diff(f1,vx);
A2y=diff(f1,vy);
A2z=diff(f1,vz);
A2=[A2x, A2y, A2z];
disp("Em função daos ângulos de Euler")
A3x=diff(f1,fi);
A3y=diff(f1,theta);
A3z=diff(f1,yaw);
A3=[A3x,A3y,A3z];
disp("Em função da velocidade angular")
A4x=diff(f1,wx);
A4y=diff(f1,wy);
A4z=diff(f1,wz);
A4=[A4x, A4y, A4z];

f2 = -skew(omg)*v - g_earth*R'*zI + (1/m)*(fa+fp);

disp("Aqui começa as derivadas do f2")
disp("Em função da posição")
A5x=diff(f2,px);
A5y=diff(f2,py);
A5z=diff(f2,pz);
A5=[A5x, A5y, A5z];
disp("Em função da velocidade")
A6x=diff(f2,vx);
A6y=diff(f2,vy);
A6z=diff(f2,vz);
A6=[A6x, A6y, A6z];
disp("Em função dos ângulos de Euler")
A7x=diff(f2,fi);
A7y=diff(f2,theta);
A7z=diff(f2,yaw);
A7=[A7x, A7y, A7z];
disp("Em função da velocidade angular")
A8x=diff(f2,wx);
A8y=diff(f2,wy);
A8z=diff(f2,wz);
A8=[A8x, A8y, A8z];

f3=Euler2Q(lbd)*omg;

disp("Aqui começa as derivadas do f3")
disp("Em função da posição")
A9x=diff(f3,px);
A9y=diff(f3,py);
A9z=diff(f3,pz);
A9=[A9x, A9y, A9z];
disp("Em função da velocidade")
A10x=diff(f3,vx);
A10y=diff(f3,vy);
A10z=diff(f3,vz);
A10=[A10x, A10y, A10z];
disp("Em função daos ângulos de Euler")
A11x=diff(f3,fi);
A11y=diff(f3,theta);
A11z=diff(f3,yaw);
A11=[A11x, A11y, A11z];
disp("Em função da velocidade angular")
A12x=diff(f3,wx);
A12y=diff(f3,wy);
A12z=diff(f3,wz);
A12=[A12x, A12y, A12z];

f4=-inv(J)*skew(omg)*J*omg + inv(J)*np;

disp("Aqui começa as derivadas do f4")
disp("Em função da posição")
A13x=diff(f4,px);
A13y=diff(f4,py);
A13z=diff(f4,pz);
A13=[A13x, A13y, A13z];
disp("Em função da velocidade")
A14x=diff(f4,vx);
A14y=diff(f4,vy);
A14z=diff(f4,vz);
A14=[A14x, A14y, A14z];
disp("Em função daos ângulos de Euler")
A15x=diff(f4,fi);
A15y=diff(f4,theta);
A15z=diff(f4,yaw);
A15=[A15x, A15y, A15z];
disp("Em função da velocidade angular")
A16x=diff(f4,wx);
A16y=diff(f4,wy);
A16z=diff(f4,wz);
A16=[A16x, A16y, A16z];

A=[A1 A2 A3 A4
   A5 A6 A7 A8
   A9 A10 A11 A12
   A13 A14 A15 A16]
%-----------------------------------------
%Matrix B
disp("Aqui começa as derivadas do f1 em função dos inputs")
B1=diff(f1,T);
B2=diff(f1,npx);
B3=diff(f1,npy);
B4=diff(f1,npz);
disp("Aqui começa as derivadas do f1 em função dos inputs")
B5=diff(f2,T);
B6=diff(f2,npx);
B7=diff(f2,npy);
B8=diff(f2,npz);
disp("Aqui começa as derivadas do f1 em função dos inputs")
B9=diff(f3,T);
B10=diff(f3,npx);
B11=diff(f3,npz);
B12=diff(f3,npz);
disp("Aqui começa as derivadas do f1 em função dos inputs")
B13=diff(f4,T);
B14=diff(f4,npx);
B15=diff(f4,npy);
B16=diff(f4,npz);
% B13=[1/Jx, 1/Jy, 1/Jz];
% B14=B13;%diff(f4,npx);
% B15=B13;%diff(f4,npy);
% B16=B14;%diff(f4,npz);

B=[B1 B2 B3 B4
   B5 B6 B7 B8
   B9 B10 B11 B12
   B13 B14 B15 B16]

save('A_matrix.mat','A')
save('B_matrix.mat','B')