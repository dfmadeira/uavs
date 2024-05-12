% crazyflie load usd card log csv file

% read file
csvfilename = '2024-04-04_log07.csv';
T = readtable(csvfilename);

% get data from table
time = table2array(T(:,1))'*1e-3;
pos = table2array(T(:,2:4))';
vel = table2array(T(:,5:7))';
% acc = table2array(T(:,8:10))';
lbd = table2array(T(:,8:10))'*pi/180;
om = table2array(T(:,11:13))'*pi/180;
pos_ref = table2array(T(:,14:16))';
yaw_ref = table2array(T(:,17))';
motors = table2array(T(:,18:21))';

% convert date to print format
t = time - time(1);
x = [pos;vel;lbd;om];
x_ref = [pos_ref;0*vel;lbd*0;om*0];
x_ref(9,:) = yaw_ref;
u = motors/65535.0;

% plot data
initPlots;
vehicle3d_ref_show_data(t,x,u,x_ref);