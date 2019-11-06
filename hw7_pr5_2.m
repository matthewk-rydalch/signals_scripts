%hw 7 problem 5.2
clear;
clc;

A = [1 2 0; 0 1 1; -1 2 0];
b = [5; 1; -1];
[lu,indx] = newlu(A);
[L, U, P] = getvals(lu,indx);

%LUx = Pb
%y = Ux
Pb = P*b;
for i = 1:1:length(indx)
    y(i) = Pb(i);
    a = i-1;
    while a > 0
        y(i) = y(i)-y(a)*L(i,a);
        a = a-1;
    end
end
for i =length(indx):-1:1
   x(i) = y(i)/U(i,i);
   a = i+1;
   while a <= length(indx)
       x(i) = x(i) - x(a)*U(i,a);
       a = a+1;
   end
end

x
check = A*x' - b


function [lu,indx] = newlu(A)
% Compute the lu factorization of A
% Copyright 1999 by Todd K. Moon

    [n,m] = size(A);
    if(n ~= m)
       error('Matrix should be square');
    end
    indx = 1:n;
    for k=1:n-1
      % For pivoting, determine largest element in this column (m = index of largest)
      [p,m] = max(abs(A(k:n,k)));
      % The previous index was out of k:n; adjust so it is indexed on 1:n
      m = m+k-1;
      % interchange the kth and mth rows
      dum = A(k,1:n);  A(k,1:n) = A(m,1:n);  A(m,1:n) = dum;
      % record which row the kth row was swapped with
      dum1 = indx(k); indx(k) = indx(m); indx(m) = dum1;
      if(A(k,k) == 0)
        error('linearly dependent columns to machine precision')
      else
        for j=k+1:n
          mult = A(j,k)/A(k,k);
          % do the row operation
          A(j,k:n) = A(j,k:n) - mult*A(k,k:n);
          % store the multiplier element in the lower triangle
          A(j,k) = mult;
        end
      end
    end
    lu= A;
end
function [L, U, P] = getvals(lu, indx)
    L = tril(lu)-diag(diag(lu))+eye(size(lu));
    U = triu(lu);
    P = zeros(size(lu));
    for i = 1:1:length(lu)
        P(i,indx(i)) = 1;
    end
end