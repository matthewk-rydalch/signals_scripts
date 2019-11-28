%hw4 pr 3.12
clear;
clc;

%i
A = [1 1; 2 1; 3 2; 5 3; 8 5; 13 8];

R = A'*A

%ii
A = [1 0; 1 1; 2 1; 3 2; 5 3; 8 5; 13 8; 0 13];

R = A'*A

%part b & c
%covariance method
A = [1 1; 2 1; 3 2; 5 3; 8 5];
R = A'*A;

p1 = A(:,1);
p2 = A(:,2);

d = [2; 3; 5; 8; 13];

P = [sum(d.*p1); sum(d.*p2)];

c = inv(R)*P;

e = d - c(1)*p1-c(2)*p2;

%autocorrelation method
A = [1 0; 1 1; 2 1; 3 2; 5 3; 8 5; 13 8; 0 13];
R = A'*A;

p1 = A(:,1);
p2 = A(:,2);

d = [1; 2; 3; 5; 8; 13; 0; 0];

P = [sum(d.*p1); sum(d.*p2)];

c = inv(R)*P;

e = d - c(1)*p1-c(2)*p2;