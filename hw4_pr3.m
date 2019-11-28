%hw4 pr 3.3
clear;
clc;

x = [2 2.5 3 5 9];
y = [-4.2 -5 2 1 24.3];

%part a
figure(1);
scatter(x,y)
title('hw3 problem 3.3 part a');

%part b
p = [x' ones(5,1)];
R = [sum(p(:,1).*p(:,1)) sum(p(:,2).*p(:,1));
     sum(p(:,2).*p(:,1)) sum(p(:,2).*p(:,2))];
P = [sum(y'.*p(:,1));
     sum(y'.*p(:,2))];

c = inv(R)*P;

y_hat = c(1)*x+c(2);

error = y-y_hat

figure(2);
scatter(x,y)
hold on
plot(x,y_hat)
title('least squares regression')

%part c
A = p;
w = diag([10 1 1 1 10]);
c = inv(A'*w*A)*A'*w*y';

y_hat = c(1)*x+c(2);

figure(3);
scatter(x,y)
hold on
plot(x,y_hat)
title('weighted least squares regression')

