%671 project.  ILS using LAMBDA method

function N = ILS(Nhat, Q_Nhat, X)
    n = size(Nhat);
    n = n(1);
    p=1; %Number of solutions on the previous level.  Initialize as 1.
    
    %reduction is a z transform to get ldl factorization, make L nearly diaganol and order D in decending order 
    %this decorrelates Q/covariance of the variables.
    [Z, L, D, Nhat, swapped] = reduction(Q_Nhat, Nhat, n);
    zhat = Z'*Nhat;
    zs = search(L,D,zhat,p,X,n); %search for valid integers.  Algorithm is only set for 2 levels
    z = optimal_z(zs, zhat, Z, Q_Nhat); %select the optimal integer z
    N = inv(Z')*z; %Fixed ambiguity for ILS
%     N
%     swapped
%     if swapped == 1
%         N_up = N(2);
%         N_down = N(1);
%         N = [N_up, N_down]';
%     end
end

function [Z, L, D, Nhat, swapped] = reduction(Q_Nhat, Nhat, n)
    [L, D] = ldl(Q_Nhat); %ldl' factorization
    Z = eye(n,n);
    k = n-1;
    k1 = n-1;
    while k>0
        if k<=k1
            for i = k+1:n
                [L,Nhat,Z] = gauss(L,i,k,Nhat,Z,n); %adjusts L, Nhat, and Z to have L nearly diagnol
            end
        end
        D_(k+1,k+1)=D(k,k)+L(k+1,k)^2*D(k+1,k+1);
        if D_(k+1,k+1)<D(k+1,k+1) %determine if D_ needs to be permutated
            [L,D,Nhat,Z] = Permute(L,D,k,D_(k+1,k+1),Nhat,Z,n); %swap rows of D to get them in decending order.  Adjust L, Nhat, and Z to match
            k1 = k;
            k = n-1;
            swapped = 1; %to tell if the result needs to be swapped back
        else
            k = k-1;
            swapped = 0;
        end
    end
end

%used to make off diagnols of L_ small
function [L, Nhat, Z] = gauss(L, i, j, Nhat, Z, n)
    mu = round(L(i,j));
    if mu ~= 0
        L(i:n,j) = L(i:n,j)-mu*L(i:n,i);
        Z(1:n,j)=Z(1:n,j)-mu*Z(1:n,i);
        Nhat(j) = Nhat(j) - mu*Nhat(i);
    end
end

% %used to make D_ in decending order
function [L,D,Nhat,Z]=Permute(L,D,k,del,Nhat,Z,n)
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
    entry1 = Nhat(k);
    entry2 = Nhat(k+1);
    Nhat(k) = entry2;
    Nhat(k+1) = entry1;
end

function newer_z = search(L,D,zhat,p,X,n)
    for i=n:-1:1
        level = i; %level (index of integers)
        if i==n
            z = zeros(2,1);
            z_(n) = zhat(n);
            [z,candidates] = get_level(D,z_,z,n,X,p,level); %get the z for the first level
        else
            for j=1:p %get a set of z's for each z in level 1
                z_ = z(:,p)-L'*(z(:,p)-zhat);
                [new_z,candidates] = get_level(D,z_,z(:,j),n,X,p,level);
                new_z = z(:,j)+new_z; %add the z elements in level 2 to those in level 1
                newer_z{j} = new_z;
            end
        end
        p = p*candidates; %increase the amount of valid solution vectors for the next level by the curent levels candidates
    end
end

function [z,p] = get_level(D,z_,z,n,X,p,i)
    [lower, upper] = bounds(D, z_,z,n,X,i); %get bounds for level
    k = ceil(lower); %first integer to test
    if k > upper
        z(i,p) = inf; %if there are no integers in bounds set vector value too high to be selected
    else
        m=1;
        while k <= upper
           z(i,m) = k; %valid z values for level
           m=m+1;
           k = k+1;
        end
        p = m-1; %set p to be the amount of vectors found on the level
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

function z = optimal_z(zs, zhat, Z, Q_Nhat)
    Q_zhat = Z'*Q_Nhat*Z; %z tranform of covariance
    S = size(zs);
    k = 1;
    for i = 1:S(2)
        len = size(zs{i});
        for j = 1:len(2) %put N's and z's in two vectors
            N_vector(k) = (zs{i}(:,j)-zhat)'*inv(Q_zhat)*(zs{i}(:,j)-zhat); %
            z_vector(:,k) = zs{i}(:,j); %possible fixed integer
            k = k+1;
        end
    end
    [N_min, z_index] = min(N_vector);
    z = z_vector(:,z_index); %z transform of fixed ambiguity
            
end

