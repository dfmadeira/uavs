function stats(filename, num)

a = csvread(filename,1,0);

figure(1);
subplot(3,1,1)
plot(a(:,1), a(:,num)) %Accelerometer x
title('Accelerometer X')
set(gca, 'FontName', 'Times New Roman')
subplot(3,1,2)
plot(a(:,1), a(:,num+1)) %Accelerometer y
title('Accelerometer Y')
set(gca, 'FontName', 'Times New Roman')
subplot(3,1,3)
plot(a(:,1), a(:,num+2)) %Accelerometer z
title('Accelerometer Z')
set(gca, 'FontName', 'Times New Roman')


figure(2);
subplot(3,1,1)
plot(a(:,1), a(:,num-3)) %Gyro x
title('Gyro X')
set(gca, 'FontName', 'Times New Roman')
subplot(3,1,2)
plot(a(:,1), a(:,num-2)) %Gyro y
title('Gyro Y')
set(gca, 'FontName', 'Times New Roman')
subplot(3,1,3)
plot(a(:,1), a(:,num-1)) %Gyro z
title('Gyro Z')
set(gca, 'FontName', 'Times New Roman')

% disp(Mean X acceloremeter)
% mean_acc_x = mean(a(:,14)); %mean
% disp(Standard Deviation of X acceloremeter)
% std_dev_acc_x = std(a(:,14)); %standard deviation
% disp(Variance of X acceloremeter)
% var_acc_x = var(a(:,14)); %variance
% 

dispstats(num,a);

end