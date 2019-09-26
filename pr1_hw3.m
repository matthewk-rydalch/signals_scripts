%ECEN 671 HW 3
clear;
clc;

t = linspace(-1, 1, 10000);
p = zeros(10000, 5);
dt = 2/10000;
for kk = 1:5
 p(:, kk) = t.^(kk-1);
end 

%step 1
q1 = p(:,1)/(sqrt(sum(p(:,1).*p(:,1))));
%step 2
e2 = p(:,2) - sum(p(:,2).*q1.*dt)*q1;
q2 = e2/(sqrt(sum(p(:,2).*p(:,2))));
%step3
e3 = p(:,3) - sum(p(:,3).*q1.*dt)*q1 - sum(p(:,3).*q2.*dt)*q2;
q3 = e3/(sqrt(sum(p(:,3).*p(:,3))));
%step 4
e4 = p(:,4) - sum(p(:,4).*q1.*dt)*q1 - sum(p(:,4).*q2.*dt)*q2 - sum(p(:,4).*q3.*dt)*q3;
q4 = e4/(sqrt(sum(p(:,4).*p(:,4))));
%step 5
e5 = p(:,5) - sum(p(:,5).*q1.*dt)*q1 - sum(p(:,5).*q2.*dt)*q2 - sum(p(:,5).*q3.*dt)*q3 - sum(p(:,5).*q4.*dt)*q4;
q5 = e5/(sqrt(sum(p(:,5).*p(:,5))));

plot(t,q1);
hold on;
plot(t,q2);
plot(t,q3);
plot(t,q4);
plot(t,q5);
legend('q1','q2','q3','q4','q5');


%check
% Legrande1 = ones(size(t));
% Legrande2 = t;
% % Legrande3 = 1/2*(3*t.^2-1);
% % Legrande4 = 1/2*(5*t.^3-3*t);
% % Legrande5 = 1/8*(35*t.^4-30*t.^2+3);
% Legrande3 = t.^2-1/3;
% Legrande4 = t.^3-3/5*t;
% Legrande5 = t.^4-5/7*t.^2+3;
% dif1 = mean(q1-Legrande1)
% dif2 = mean(q2-Legrande2)
% dif3 = mean(q3-Legrande3)
% dif4 = mean(q4-Legrande4)
% dif5 = mean(q5-Legrande5)