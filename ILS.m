%671 project
clear;
clc;

%reduction process
ahat = [5.45, 3.10, 2.97]'; %real value least squares estimate of a
Q_ahat = [1.5 2.4 -0.9;
      -1.2 3.1 -2.1;
      0.2 1.1 -0.4];
Z = [1 1 1;
     1 1 2;
     2 1 3]; %an integer unimodular matrix (abs(det(z)) = 1

n = size(ahat);
n = n(1);
p=2;
X2 = 100; 

[Z, L, D, ahat] = reduction(Q_ahat, ahat, n);
% Optis = mysearch(L,D,Qa,X2,ahat,n);
% Optis = search2(L,D,Q_ahat,X2,ahat,ahat,n);
zhat = Z'*ahat;
Optis = msearch(L,D,zhat,p,n)
[Lz,Uz,Pz] = lu(Z')
y = Lz\(Pz*Optis);
a = Uz\y

function optis = mysearch(L,D,Qa,X2,ahat,n)
    upper_bound = ahat(i)-sqrt(right(i))-sum
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
    Z = eye(n,n);
    k = n-1;
    k1 = n-1;
    while k>0
        if k<=k1
            for i = k+1:n
                [L,ahat,Z] = gauss(L,i,k,ahat,Z,n);
            end
        end
        D_(k+1,k+1)=D(k,k)+L(k+1,k)^2*D(k+1,k+1);
        if D_(k+1,k+1)<D(k+1,k+1)
            [L,D,ahat,Z] = Permute(L,D,k,D_(k+1,k+1),ahat,Z,n);
            k1 = k;
            k = n-1;
        else
            k = k-1;
        end
    end
end

%used to make off diagnols of L_ small
function [L, ahat, Z] = gauss(L, i, j, ahat, Z, n)
    mu = L(i,j);
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

