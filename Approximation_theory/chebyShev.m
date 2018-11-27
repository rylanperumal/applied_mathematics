function chebyShev(f, x, n)

w = @(x) 1./sqrt(1 - x.^2);
xs = @(x) x;
g = zeros(1, length(x));

for i = 1:n
    integrand = @(xs) w(xs).*f(xs).*genPhins(xs, i);
%     integrand
    g = g + integral(integrand, -1, 1) * genPhins(x, i);
end
plot(x, f(x), x, g);
end

function phin = genPhins(x, n)
% x is a row vector 
% n is the order

if n == 0
    phin = 1 / sqrt(pi) * ones(1, length(x));
    return
end
phin = sqrt(2 / pi) * genTns(x, n);
end
function Tn = genTns(x, n)
% x is a row vector
% n is the order
if n == 0
    Tn = ones(1, length(x));
    return
end
if n == 1
    Tn = x;
    return
end
Tn = (2 * x.*genTns(x, n - 1)) - genTns(x, n - 2);
end