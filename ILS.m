%671 project.  ILS using LAMBDA method
clear;
clc;

ahat = [401.6, 358.33]'; %real value least squares estimate of N
Q_ahat = [1.6 0.75; %covariance.  Symetric square positive semi-def
      0.75 1.9];

n = size(ahat);
n = n(1);
p=1; %Number of solutions.  Best and second best solutions
X2 = 100; %elipse space to search squared

%z transform to get ldl factorization, make L nearly diaganol and order D in decending order 
%this decorrelates Q/covariance of the variables.
[Z, L, D, ahat] = reduction(Q_ahat, ahat, n);
% Optis = search2(L,D,Q_ahat,X2,ahat,ahat,n);
zhat = Z'*ahat;
% z = msearch(L,D,zhat,p,n) %MLAMBDA search method.
zs = mysearch(L,D,zhat,p,n);
z = optimal_z(zs, zhat, Z, Q_ahat);
N = inv(Z')*z
% [Lz,Uz,Pz] = lu(Z');
% y = Lz\(Pz*z);
% N = Uz\y

function newer_z = mysearch(L,D,zhat,p,n)
    X = 20; %chi for search region
    for i=n:-1:1
        level = i;
        if i==n
            z = zeros(2,1);
            z_(n) = zhat(n);
            [z,candidates] = get_level(D,z_,z,n,X,p,level);
        else
            for j=1:p
                z_ = z(:,p)-L'*(z(:,p)-zhat);
                [new_z,candidates] = get_level(D,z_,z(:,j),n,X,p,level);
                new_z = z(:,j)+new_z;
                newer_z{j} = new_z;
            end
        end
        p = p*candidates;
    end
end

function [z,p] = get_level(D,z_,z,n,X,p,i)
    [lower, upper] = bounds(D, z_,z,n,X,i);
    k = ceil(lower); %integer to test
    if k > upper
        z(i,p) = 1000;
    else
        m=1;
        while k <= upper
           z(i,m) = k;
           m=m+1;
           k = k+1;
        end
        p = m-1; %set p to be the amount of vectors found on the first level
    end
end
function [lower, upper] = bounds(D, z_,z,n,X,i)
    S = 0; %summation used in bounds
    if i==n
        lower = z_(n)-D(n,n)^(-1/2)*X;
        upper = z_(n)+D(n,n)^(-1/2)*X;
    else
        for j=i+1:n
            S = S+(z(j)-z_(j))^2/D(j,j);
        end
        lower = z_(i)-D(i,i)^(1/2)*(X^2-S)^(1/2);
        upper = z_(i)+D(i,i)^(1/2)*(X^2-S)^(1/2);
    end
end

function z = optimal_z(zs, zhat, Z, Q_ahat)
    Q_zhat = Z'*Q_ahat*Z;
    S = size(zs);
    k = 1;
    for i = 1:S(2)
        len = size(zs{i});
        for j = 1:len(2)
            solution(k) = (zs{i}(:,j)-zhat)'*inv(Q_zhat)*(zs{i}(:,j)-zhat);
            z_index(:,k) = zs{i}(:,j);
            k = k+1;
        end
    end
    [z_min, index] = min(solution);
    z = z_index(:,index);
            
end
function optis = search2(L,D,Qa,X2,a,ahat,n)
    dist = a-ahat;
    right(n+1) = X2; %right bound
    left(n+1)=0; %left bound
    ende = false; %end search
    ncan = 0; %number of candidates
    for i=1:n-1 %calculate change in D for each element on diagnol
        dq(i) = D(i+1,i+1)/D(i,i);
    end
    dq(n) = 1/D(i,i);
    i = n+1;
    while ende == false %search until ende is true
        i = i-1;
        lef(i) = dot(L(i+1:n,i),dist(i+1:n));
        right(i) = (right(i+1)-left(i+1))*dq(i);
        reach = sqrt(right(i));
        del = a(i)-reach-lef(i);
        dist(i) = ceil(del)-a(i);
        if dist(i) > reach - lef(i)
            %what is en?
            [dist, left, ende] = BACKTS(n, i, en, dist, lef, left, ende); %backtrack
        else
            en(i) = reach-lef(i)-1;
            left(i)=(dist(i)+lef(i))^2;
        end
        if i == 1
            cands = COLLECT(n, 2, D, lef, left, right, dist, en, ncan, disall, cands, tmax, imax);
            [dist, left, ende] = BACKTS(n, i, en, dist, lef, left, ende);
        end
    end
    
end

function [dist, left, ende] = BACKTS(n, i, en, dist, lef, left, ende)
    j = i+1;
    for i = j:n
        if dist(i) <= en(i)
            dist(i) = dist(i)+1;
            left(i)=(dist(i)+lef(i))^2;
        else
            ende = true;
        end
    end
end
            
function cands = COLLECT(n, maxcan, D, lef, left, right, X2, dist, en, ncan, disall, cancds, tmax, imax)
    t = X2-(right(1)-left(1))*D(1,1);
    en(1) = en(1)+1;
    while dist(1) <= en(1)
        ncan = ncan+1;
        if ncan<maxcan
            cands = STOREs(ncan,ncan,imax,t,tmax,dist,cands,disall)
        else
            if t < tmax
                cands = STOREs(maxcan, imax, imax, t, tmax, dist, cands, disall)
            end
        end
        t=t+(2*(dist(1)+lef(1))+1)*D(1,1);
        dist(1) = dist(1)+1;
    end
end

function cands = STOREs(ican, ipos, imax, t, tmax, dist, cands, disall)
    for i = 1:n
        cands(i,ipos) = dist(i);
    end
    disall(ipos) = t;
    tmax = t;
    imax = ipos;
    for i = 1:ican
        if disall(i) > tmax
            imax=i;
            tmax=disall(i);
        end
    end
end

function [Z, L, D, ahat] = reduction(Q_ahat, ahat, n)
    [L, D] = ldl(Q_ahat); %ldl' factorization
%     L = L'
    Z = eye(n,n);
    k = n-1;
    k1 = n-1;
    while k>0
        if k<=k1
            for i = k+1:n
                [L,ahat,Z] = gauss(L,i,k,ahat,Z,n); %adjusts L, ahat, and Z to have L nearly diagnol
            end
        end
        D_(k+1,k+1)=D(k,k)+L(k+1,k)^2*D(k+1,k+1);
        if D_(k+1,k+1)<D(k+1,k+1) %determine if D_ needs to be permutated
            [L,D,ahat,Z] = Permute(L,D,k,D_(k+1,k+1),ahat,Z,n); %swap rows of D to get them in decending order.  Adjust L, ahat, and Z to match
            k1 = k;
            k = n-1;
        else
            k = k-1;
        end
    end
end

%used to make off diagnols of L_ small
function [L, ahat, Z] = gauss(L, i, j, ahat, Z, n)
    mu = round(L(i,j));
    if mu ~= 0
        L(i:n,j) = L(i:n,j)-mu*L(i:n,i);
        Z(1:n,j)=Z(1:n,j)-mu*Z(1:n,i);
        ahat(j) = ahat(j) - mu*ahat(i);
    end
end

% %used to make D_ in decending order
function [L,D,ahat,Z]=Permute(L,D,k,del,ahat,Z,n)
    eta = D(k,k)/del; %see eq 8
    lambda = D(k+1,k+1)*L(k+1,k)/del; %see eq 8
    D(k,k) = eta*D(k+1,k+1); %see eq 8
    D(k+1,k+1)=del;
    L(k:k+1,1:k-1)=[-L(k+1,k) 1;
                    eta lambda]*L(k:k+1,1:k-1); %see eq 10
    L(k+1,k) = lambda;
    %swap columns and entries, see eq 11
    column1 = L(k+2:n,k);
    column2 = L(k+1:n,k+1);
    L(k+2:n,k) = column1;
    L(k+1:n,k+1) = column2;
    column3 = Z(1:n,k);
    column4 = Z(1:n,k+1);
    Z(1:n,k) = column4;
    Z(1:n,k+1) = column3;
    entry1 = ahat(k);
    entry2 = ahat(k+1);
    ahat(k) = entry2;
    ahat(k+1) = entry1;
end

function Optis = msearch(L,D,zhat,p,n)
    maxDist = inf; %current X^2, starts around inf
    k = n;
    dist(k) = 0; 
    endSearch = false;
    count = 0; %count: the number of candidates
    S = zeros(n,n); %will be used for computing z_
    z_(n) = zhat(n); %see eq 17
    z(n) = z_(n);
    y = z_(n)-z(n);
    step(n) = sign(y);
    imax = p;
    while endSearch == false
        newDist=dist(k)+y^2/D(k,k);
        if newDist < maxDist
            if k~=1 %move down
                k = k-1;
                dist(k) = newDist;
                S(k,1:k) = S(k+1,1:k)+(z(k+1)-z_(k+1))*L(k+1,1:k);
                z_(k)=zhat(k)+S(k,k); %eq 17
                z(k)=z_(k);
                y=z_(k)-z(k);
                step(k)=sign(y);
            else %store the candidate and try next valid integer
                if count < p-1
                    count = count+1;
                    Optis(:,count)=z(1:n);
                    fun(count) = newDist; %store f(z)
                else
                    Optis(:,imax) = z(1:n);
                    fun(imax)=newDist;
                    [mvalue,imax] = max(fun);
                    maxDist = fun(imax);
                end
                z(1) = z(1)+step(1);
                %next valid integer
                y=z_(1)-z(1);
                step(1) = -step(1)-sign(step(1));
            end
        else
            if k == n
                endSearch = true;
            else
                k=k+1; %move up
                z(k) = z(k)+step(k);
                %next valid integer
                y = z_(k)-z(k);
                step(k) = -step(k)-sign(step(k));
            end
        end
    end
end

