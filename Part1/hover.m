function hover(A, B, C, D, Eq_Thrust, figureno)

    % simulation parameters
    nx = 12;
    ny = 4;
    Dt = 0.01;
    t = 0:Dt:10;
    %u_eq = Eq_Thrust*ones(size(t)); %divided by 4 to consider all motors;
    u_eq = [[0.027*9.81-0.2639*cosd(5);0;0;0]*ones(size(0:Dt:1)) zeros(4,900)]
    % simulate nonlinear system
    Nsim = length(t);
    x_init = [0;0;1.5;0;0;0;0;0;0;0;0;0];
    x_eq = [0;0;1.5;3.387;0;0;0;5*pi/180;0;0;0;0];%*ones(size(t));
    %y = zeros(ny,Nsim);
    
    sys=ss(A,B,C,D);
    y=lsim(sys,u_eq, t,x_eq);%,x_eq);
    
    
    figure(figureno)
    subplot(1, 2,1)
    plot(t, y(:,1),t, y(:,2),t, y(:,3));
    legend(["x" "y" "z"])
    xlabel('Time')
    ylabel('Position')
    title('Response of State Variables')
    grid on
    get(gca,'fontname')  % shows you what you are using.
    set(gca,'fontname','times')  % Set it to times

    subplot(1, 2,2)
    plot(t, y(:,4),t,y(:,5),t, y(:,6));
    %axis([0 5])
    legend(["vx" "vy" "vz"])
    xlabel('Time')
    ylabel('Position')
    title('Response of State Variables')
    grid on
    get(gca,'fontname')  % shows you what you are using.
    set(gca,'fontname','times')  % Set it to times

    % figure(figureno+2)
    % plot(t, u_eq(1,:))
end