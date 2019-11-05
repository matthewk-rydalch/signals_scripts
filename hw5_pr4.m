clear;
clc;

A1 = [4 3; 3 6];
A2 = [1 2; 3 0];
A3 = [1 2; 0 1];

%solving l2 norms
A1_norm = l2_norm(A1)
A2_norm = l2_norm(A2)
A3_norm = l2_norm(A3)


function A_norm = l2_norm(A)
    lambda = eig(A'*A);
    ro = max(lambda);
    A_norm = sqrt(ro);
end