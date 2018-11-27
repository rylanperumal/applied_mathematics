function out = Householders(A)
%
% INPUT -> Takes in a symmetric diagonal matrix A
% OUTPUT -> Produces a tridiagonal matrix Abar
%

n = length(A);
stop = n - 2; % we only want a tridiagonal matrix
start = 2;
for r = 1 : stop
    zero = zeros(r, n - r);
    Ip = eye(r, r);
    x = transpose(A(r, start:n));
    alpha = norm(x);
    e = zeros(length(x), 1);
    e(1, 1) = 1;
    if(x(1) > 0)
        alpha = -alpha;
    end
    u = x - alpha * e;
    w = (1 / norm(u)) * u;
    H = eye(length(x)) - (2 * w * transpose(w));
    P = [Ip, zero; transpose(zero), H];
    A = P * A * P;
    start = start + 1;
end
out = A;
end
