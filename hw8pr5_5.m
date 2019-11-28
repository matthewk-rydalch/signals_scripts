%hw8 problem 5.5
clear;
clc;

A = [4 6 4; 6 25 18; 4 18 22];

L = chol(A, 'lower');

A_lo = L*L'

D = diag([L(1,1)^2, L(2,2)^2, L(3,3)^2]);
U = inv(sqrt(D))*L';

A_up = U'*D*U