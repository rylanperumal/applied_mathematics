% Fitzhugh Nagumo equation with Dirichlet bounday conditions using
% the Crank Nicolson semi-implicit scheme

N = 100; % this is equal to jmax
x = linspace(0, 1, N); % domain vector
beta = 0.5;
D = 1;
dx = x(2) - x(1); % space between each point in the domain
dt = (dx ^ 2) / (2 * D); % space between timesteps 
timeSteps = (0:dt:1); 
T = length(timeSteps);
% Parts of the A matrix
r = (dt / dx.^2);
a = -beta * r;
b = 1 + r;
c = -beta *  r;

jmax = N - 2; % 2 less the total length
un = zeros(1, N);
% initial conditions:
for i = 1:N
    if (x(i) >= 0 && x(i) <= 0.5)
        un(i) = 0;
    end
    if (x(i) > 0.5 && x(i) <= 1)
        un(i) = 1;
    end
end
unp1 = zeros(1, jmax);
A = diag( b * ones(1, N-2)) + diag(a * ones(1, N-3), -1) + diag(c * ones(1,N-3), 1); %Tridiagonal
PlotMat = zeros(T, N); % matrix for 3D plot
PlotMat(1, :) = un;
d = zeros(1, jmax);
for j = 1:T
    % This loop is used to build the d matrix, from Au = d
    for i = 2:jmax    
        if i == 1
             d(i) = un(2) - r*un(2) + beta*r*un(1) + beta*r*un(3) + beta*r*0; 
        elseif i == jmax
            d(jmax) = un(end-1) - r*un(end-1) + beta*r*un(end-2) + beta*r*un(end) + beta*r*1;
        else        
        d(i) = un(i+1) - r*un(i+1) + beta*r*un(i) + beta*r*un(i+2); % Row vector
        end
    end
    unp1 = d * inv(A) + (dt*(unp1.*(1- unp1).*(unp1 - 0.3)));
    un(1, 2:jmax + 1) = unp1;
    PlotMat(j, :) = un;
    figure(1);
    plot(x, un);
    title('Crank-Nicolson scheme using Dirichlet boundary conditions')
    ylabel('u(x, t)') 
    xlabel('x')
    axis([0 1 0 1]);
    pause(0.01);
end
% 3D plot
figure()
[xplot, timeplot] = meshgrid(x, timeSteps);
surf(xplot, timeplot, PlotMat);
axis([0 1 0 1 0 1])
title('Crank-Nicolson scheme using Dirichlet boundary conditions')
zlabel('u(x, t)') 
xlabel('x')
ylabel('t');
shading interp

