function [out, theta] = Givens(A)
%
% INPUT -> symmetric matrix A 
% OUTPUT -> tridiagonal matrix
%

% we start with position (1, 3) -> (ihat, j)
n = length(A);
start = n - 1;
stop = n - 2;
for ihat = 1 : stop
    for j = start : n
        % we want to create a zero in position ihat, j
        i = ihat + 1;
        theta = -atan(((A(ihat, j)) / (A(ihat, i))));
        P = eye(n);
        P(i, i) = cos(theta);
        P(i, j) = -sin(theta);
        P(j, i) = sin(theta);
        P(j, j) = cos(theta);
        A = P * A * transpose(P);
    end
    start = start + 1;
end
out = A;
end