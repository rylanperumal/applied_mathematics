function chebyShev = genChebyShev(f, x, n)
% f is the function it is trying to approximate
% xs is the number of x points
% n is the order

w = @(x) 1./sqrt(1 - x.^2); % weight function, will be specified
xs = @(x) x; % note: integeral function can only take in a function handle
g = zeros(1, length(x)); % vector of points defining the approximated function
for i = 1:n % we want to sum up to order n
    integrand = @(xs) w(xs).*f(xs).*genPhis(i, xs); 
    g = g + integral(integrand, -1, 1)*genPhis(i, x);
end
chebyShev = g;
end