
X= csvread('dados_recebidos.csv',1,0);

figure(1)
axis([0 1 0 18]);
plot(X(:,1))
