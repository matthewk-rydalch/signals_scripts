%ILS_test_driver
clear;
clc;

samples = 30;

%Values for low covariance case
Nhat = [14.6, 103.33]'; %real value least squares estimate of N.  Arbitraryily selected for demonstration
Q_Nhat = [.5 0; %covariance.  Symetric square positive semi-def.  Arbitrarily selected for demonstration
      0 .5];
X = 15; %chi for search region

% changes in N for each test
n = 100; %number of values to test
d_N = .05;

%changes in covariance for "with covariance case
d_var = 0.3;
d_cov = 0.3;
cov = Q_Nhat+[d_var d_cov; d_cov d_var];

for i=1:n
    Nhat_test(:,i) = Nhat+(i-1)*[-d_N, d_N]'; %float estimate
    rounding(:,i) = round(Nhat_test(:,i)); %rounding estimate
    N_normal(:,i) = ILS(Nhat_test(:,i), Q_Nhat, X); %low covariance estimate
    N_cov(:,i) = ILS(Nhat_test(:,i), cov, X); %with covariance estimate
    dv(i) = i; %just storing the iteration step for plotting
    i
end

%plotting
figure(1)
hold on
plot(dv, Nhat_test(2,:))
plot(dv, rounding(2,:), '*')
plot(dv, N_normal(2,:))
plot(dv, N_cov(2,:))
legend('Float Estimate', 'Rounded Estimate', 'ILS low covariance', 'ILS with covariance')
xlabel('Time step')
ylabel('Estimate')
