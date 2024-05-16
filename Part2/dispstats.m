function dispstats(num,a)
%% acc
disp("Mean X acceloremeter")
disp(mean(a(:,num))); %mean
disp("Standard Deviation of X acceloremeter")
disp(std(a(:,num))); %standard deviation
disp("Variance of X acceloremeter")
disp(var(a(:,num))); %variance

num=num+1;
disp("Mean of Y acceloremeter")
disp(mean(a(:,num))); %mean
disp("Standard Deviation of Y acceloremeter")
disp(std(a(:,num))); %standard deviation
disp("Variance of Y acceloremeter")
disp(var(a(:,num))); %variance

num=num+1;
disp("Mean Z acceloremeter")
disp(mean(a(:,num))); %mean
disp("Standard Deviation of Z acceloremeter")
disp(std(a(:,num))); %standard deviation
disp("Variance of Y acceloremeter")
disp(var(a(:,num))); %variance

%% Gyro
num=num-5;
disp("Mean X Gyro")
disp(mean(a(:,num))); %mean
disp("Standard Deviation of X Gyro")
disp(std(a(:,num))); %standard deviation
disp("Variance of X Gyro")
disp(var(a(:,num))); %variance

num=num+1;
disp("Mean of Y Gyro")
disp(mean(a(:,num))); %mean
disp("Standard Deviation of Y Gyro")
disp(std(a(:,num))); %standard deviation
disp("Variance of Y Gyro")
disp(var(a(:,num))); %variance


num=num+1;
disp("Mean Z Gyro")
disp(mean(a(:,num))); %mean
disp("Standard Deviation of Z Gyro")
disp(std(a(:,num))); %standard deviation
disp("Variance of Y Gyro")
disp(var(a(:,num))); %variance

end