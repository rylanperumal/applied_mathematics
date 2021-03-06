% Fitzhugh Nagumo equation with Neumann bounday conditions using
% the explicit finite difference scheme

% N = 15 for simulation
N = 30;% number of spacial points, set low so the simulation in quick can be higher
x = linspace(0, 1, N); % x number line
D = 1; % diffusivity constant
dx = x(2) - x(1); % the distance between the spacial points
dt = (dx^2) / (4 * D);
timeSteps = (0:dt:1);
T = length(timeSteps);
un = zeros(1, N);
% setting the initial condition
for i = 1:N
    if (x(i) > 0.5 && x(i) <= 1)
        un(i) = 1;
    end
    if (x(i) >= 0 && x(i) <= 0.5)
        un(i) = 0;
    end
end

unp1 = un;
PlotMat = zeros(T, N); % matrix for surface plot
PlotMat(1, :) = unp1;
lambda1 = (D * dt) / (dx^2);
lambda2 = 1 - 2*lambda1;
A = diag(lambda2 * ones(1, N)) + diag(lambda1 * ones(1, N-1), -1) + diag(lambda1 * ones(1, N-1), 1);
A(1, 2) = 2 * lambda1;
A(end, end-1) = 2 * lambda1;
b = zeros(N, 1);
alpha = 0;
beta = 0;
b(1) = -2*lambda1*dx*alpha;
b(N) = 2*lambda1*dx*beta;

for i = 2 : T
    unp1 = A * un' + b + (dt*(un.*(1-un).*(un-0.3)))';
    un = unp1';
    PlotMat(i, :) = un;
%     figure(1);
%     plot(x, unp1);
%     title('Fitzhugh-Nagumo equation using Neumann boundary conditions')
%     ylabel('u(x, t)') 
%     xlabel('x')
%     axis([0 1 0 1]);
%     pause(0.01);
end
% 3D plot
figure() 
[xplot, timeplot] = meshgrid(x, timeSteps);
surf(xplot, timeplot, PlotMat);
axis([0 1 0 1 0 1])
title('Fitzhugh-Nagumo equation using Neumann boundary conditions')
zlabel('u(x, t)') 
xlabel('x')
ylabel('t');
shading interp

    
    
    

