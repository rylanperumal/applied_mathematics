% Fitzhugh Nagumo equation with Neumann bounday conditions using
% the Crank Nicolson semi-implicit scheme
N = 10; % Number of spacial points
x = linspace(0, 1, N); % domain vector
D = 1;
dx = x(2) - x(1); % sapce between each point in the domain
dt = (dx.^2) / (2*D); % "time" between each step in time
timeSteps = 0:dt:1;
T = length(timeSteps);
 
x = [x(1) - dx, x, x(N) + dx];
ghostN = N + 2;

beta = 0.5;
s = (dt/dx.^2);

un = zeros(1, ghostN);
% initial condition
for i = 1:ghostN
    if (x(i) >= 0 && x(i) <= 0.5)
        un(i) = 0;
    end
    if (x(i) > 0.5 && x(i) <= 1)
        un(i) = 1;
    end
end
unp1 = un(1, 2:N+1);
% Building A matrix
A = diag( (1 + 2*beta*s) * ones(1, N));
A = A + diag((-beta*s)*ones(1, N-1), 1) + diag((-beta*s)*ones(1, N-1), -1);
A(1, 1) = 1 + 2*s;
A(end, end) = 1 + 2*s;
A(1, 2) = -2*s;
A(end, end-1) = -2*s;

% Matrix for 3D plot
PlotMat = zeros(T, N);
PlotMat(1, :) = unp1;
% Neuman Boundary values
cn1 = 0; % alpha
cn2 = 0; % beta

d = zeros(1, N);
for i = 2:T
    % Creating our d matrix
    for j = 1:N
        if j == 1
            d(1) = un(2) - 2*s*dx*cn1;
        elseif j == N
            d(N) = un(ghostN-1) + 2*s*dx*cn2;
        else
            d(j) = un(j+1) - 2*(1-beta)*s*un(j+1) + (1-beta)*s*un(j) + (1-beta)*s*un(j+2);
        end
        
    end
    nonlinU = un(1, 2:N+1);
    % update rule
    unp1 = inv(A) * d' + (dt*(nonlinU.*(1 - nonlinU).*(nonlinU - 0.3)))';
    un(1, 2:N+1) = unp1';
    PlotMat(i, :) = unp1;
    % simulating the equation
    figure(1);
    plot(x, un);
    title('Crank-Nicolson scheme using Neumann boundary conditions')
    ylabel('u(x, t)') 
    xlabel('x')
    axis([0 1 0 1]);
    pause(0.1);    
end
% 3D plot of solution
x = x(1, 2:ghostN-1);
figure()
[xs, timeplot] = meshgrid(x, timeSteps);
surf(xs, timeplot, PlotMat);
axis([0 1 0 1 0 1])
title('Crank-Nicolson scheme using Neumann boundary conditions')
zlabel('u(x, t)') 
xlabel('x')
ylabel('t');
shading interp

    