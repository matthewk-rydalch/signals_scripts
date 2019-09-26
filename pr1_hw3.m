%ECEN 671 HW 3
clear;
clc;

part = 'b'

switch part
    case 'a'
        [q1, q2, q3, q4, q5, t, dt] = gram_schmidt();
        plotter_A(q1, q2, q3, q4, q5, t);
    case 'b'
        [q1, q2, q3, q4, q5, t, dt] = gram_schmidt();
        f = exp(-t);
        f_hat = regression(q1, q2, q3, q4, q5, f', dt);
        error = f_hat-f;
        er_norm = sum(error.*error.*dt)
        plotter_B(f, f_hat, t)
        
    case 'c'
        
    case 'd'
        
end

function [q1, q2, q3, q4, q5, t, dt] = gram_schmidt()
    t = linspace(-1, 1, 10000);
    p = zeros(10000, 5);
    dt = 2/10000;
    for kk = 1:5
     p(:, kk) = t.^(kk-1);
    end 

    %step 1
    q1 = p(:,1)/(sqrt(sum(p(:,1).*p(:,1)*dt)));
    %step 2
    e2 = p(:,2) - sum(p(:,2).*q1.*dt)*q1;
    q2 = e2/(sqrt(sum(e2.*e2*dt)));
    %step3
    e3 = p(:,3) - sum(p(:,3).*q1.*dt)*q1 - sum(p(:,3).*q2.*dt)*q2;
    q3 = e3/(sqrt(sum(e3.*e3.*dt)));
    %step 4
    e4 = p(:,4) - sum(p(:,4).*q1.*dt)*q1 - sum(p(:,4).*q2.*dt)*q2 - sum(p(:,4).*q3.*dt)*q3;
    q4 = e4/(sqrt(sum(e4.*e4.*dt)));
    %step 5
    e5 = p(:,5) - sum(p(:,5).*q1.*dt)*q1 - sum(p(:,5).*q2.*dt)*q2 - sum(p(:,5).*q3.*dt)*q3 - sum(p(:,5).*q4.*dt)*q4;
    q5 = e5/(sqrt(sum(e5.*e5.*dt)));
end

function[f_hat] = regression(q1, q2, q3, q4, q5, f, dt)
    R = [sum(q1.*q1.*dt), sum(q1.*q2.*dt), sum(q1.*q3.*dt), sum(q1.*q4.*dt), sum(q1.*q5.*dt); ...
        sum(q2.*q1.*dt), sum(q2.*q2.*dt), sum(q2.*q3.*dt), sum(q2.*q4.*dt), sum(q2.*q5.*dt); ...
        sum(q3.*q1.*dt), sum(q3.*q2.*dt), sum(q3.*q3.*dt), sum(q3.*q4.*dt), sum(q3.*q5.*dt); ...
        sum(q4.*q1.*dt), sum(q4.*q2.*dt), sum(q4.*q3.*dt), sum(q4.*q4.*dt), sum(q4.*q5.*dt); ...
        sum(q5.*q1.*dt), sum(q5.*q2.*dt), sum(q5.*q3.*dt), sum(q5.*q4.*dt), sum(q5.*q5.*dt)];

    P = [sum(f.*q1.*dt); sum(f.*q2.*dt); sum(f.*q3.*dt); sum(f.*q4.*dt); sum(f.*q5.*dt)];
    
    C = P'*inv(R);
    
    f_hat = C*[q1'; q2'; q3'; q4'; q5'];
end

% function[x] = in(a,b,dt)
%     x = sum(a.*b.*dt)
% end

function[] = plotter_A(q1, q2, q3, q4, q5, t)
    plot(t,q1);
    hold on;
    plot(t,q2);
    plot(t,q3);
    plot(t,q4);
    plot(t,q5);
    legend('q1','q2','q3','q4','q5');
end

function[] = plotter_B(f, f_hat, t)
    plot(t,f);
    hold on;
    plot(t,f_hat);
    legend('f','estimation');
end

