clear;
clc;

%init
n = 5;
steps = 100;
impulse = [1 2 3 4 5]';
del = .1;
Pp = 1/del*eye(n);
hp = zeros(n,1);
qt = zeros(n,1);
qd = zeros(5,1);

for t=1:1:steps
    f(t) = sin(t);
end

for t = 1:1:steps
    
    qd(5) = qd(4);
    qd(4) = qd(3);
    qd(3) = qd(2);
    qd(2) = qd(1);
    qd(1) = f(t);
    
    qt = qd(1:n);
    
    dt = qd'*impulse;
    kt = Pp*qt/(1+qt'*Pp*qt);
    eps = dt-qt'*hp;
    ht = hp+kt*eps;
    Pt = Pp-kt*qt'*Pp;
    
    Pp = Pt;
    hp = ht;
    
end
    

