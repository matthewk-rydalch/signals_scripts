%hw 7 problem 5.4
clear;
clc;

A = [2 4 -5; 6 12.001 1; 4 -8 -3];
b = [-5; 33.002; -21];

%no permutation
[U_no, L_no] = no_perm(A,b);
x_no = solve_x(U_no, L_no, b)

%with permutation
[U_with, L_with, P] = with_perm(A,b);
x_with = solve_x(U_with, L_with, P*b)

function [U, L] = no_perm(A,b)
    % get U
    temp = round(A,3);
    E1 = [1 0 0; -temp(2,1)/temp(1,1), 1, 0; -temp(3,1)/temp(1,1), 0, 1];
    temp = round(E1*A,3);
    E2 = [1 0 0; 0 1 0; 0 -temp(3,2)/temp(2,2) 1];
    temp = round(E2*E1*A,3);
    U = temp;

    %get L
    L = round(inv(E2*E1),3);
end

function [U, L, P] = with_perm(A,b)
    % get U
    temp = A;
    P1 = find_perm(temp, 1);
    temp = round(P1*A,3);
    E1 = [1 0 0; -temp(2,1)/temp(1,1), 1, 0; -temp(3,1)/temp(1,1), 0, 1];
    temp = round(E1*P1*A,3);
    P2 = find_perm(temp, 2);
    temp = round(P2*E1*P1*A,3);
    E2 = [1 0 0; 0 1 0; 0 -temp(3,2)/temp(2,2) 1];
    temp = round(E2*P2*E1*P1*A,3);
    U = temp;

    %get L
    P = P2*P1;
    L = round(P*inv(E2*P2*E1*P1),3);
end

function P = find_perm(temp, i)
    [M,p1] = max(abs(temp(i:3,i)));
    p = eye(3,3);
    P = eye(3,3);
    P(i,:) = p(p1+i-1,:);
    P(p1+i-1,:) = p(i,:);
end


function x = solve_x(U,L,b)
    %LUx = b
    %y = Ux
    for i = 1:1:3
        y(i) = b(i);
        a = i-1;
        while a > 0
            y(i) = round(y(i)-y(a)*L(i,a),3);
            a = a-1;
        end
    end
    for i =3:-1:1
       x(i) = round(y(i)/U(i,i),3);
       a = i+1;
       while a <= 3
           x(i) = round(x(i) - x(a)*U(i,a),3);
           a = a+1;
       end
    end

    x;
end
