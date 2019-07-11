clear 

%% Define relevant parameters.   

m = -1; % minimum value of f
M = 1; % maximum value of f
R_1 = 0.5;
R_2 = 8*R_1;
M_0 = 8*R_2*M;

plotW = 4*R_2;
plotH = max(M*plotW*2, 1.5*M_0);

N = 1000;
x = linspace(-plotW, plotW, N);

%% Approximate the index for -R_1, R_1, -R_2, R_2
dx = 2*plotW/N;
n1 = round((plotW-R_1)/dx); % index of -R_1
n2 = round((plotW - R_2)/dx); % index of - R_2
N1 = round((plotW+R_1)/dx); % index of R_1
N2= round((plotW + R_2)/dx); % index of R_2

%% Plot the lines y = mx and y = Mx in regions (-R_1, R_1), and x<R-R_2, x>R_2

figure(1)
plot(x(1:n2), M*x(1:n2), 'g');
hold on
plot(x(n1:N1), m*x(n1:N1), 'g');
plot(x(N2:N), M*x(N2:N), 'g');
plot(x, 0*x, 'k');
axis([-plotW plotW -plotH plotH]);

%% plot vertical lines for boundaries
line(x(n2)*ones(1, 100), linspace(M*x(n2), -plotH), 'Color', 'g');
line(x(n1)*ones(1, 100), linspace(m*x(n1), plotH), 'Color', 'g');
line(x(N2)*ones(1, 100), linspace(M*x(N2), -plotH), 'Color', 'g');
line(x(N1)*ones(1, 100), linspace(m*x(N1), plotH), 'Color', 'g');

%% Fill in the regions

for i=1:floor(n2/4)
  line(x(4*i)*ones(1, 100), linspace(M*x(4*i), -plotH), 'Color', 'g');
end

for i=1:floor((N1-n1)/4)
  line(x(n1+4*i - 1)*ones(1, 100), linspace(m*x(n1+4*i - 1), plotH), 'Color', 'g');
end

for i=1:floor((N-N2)/4)
  line(x(N2+4*i - 1)*ones(1, 100), linspace(M*x(N2+4*i - 1), -plotH), 'Color', 'g');
end

%% Plot the bounds (-M_0, M_0) on g between -R_2 and R_2

plot(x(n2:N2), M_0*ones(1, N2-n2+1), 'y');
plot(x(n2:N2), -M_0*ones(1, N2-n2+1), 'y');
line(x(n2)*ones(1, 100), linspace(M_0, plotH), 'Color', 'y');
line(x(n2)*ones(1, 100), linspace(-M_0, -plotH), 'Color', 'y');
line(x(N2)*ones(1, 100), linspace(M_0, plotH), 'Color', 'y');
line(x(N2)*ones(1, 100), linspace(-M_0, -plotH), 'Color', 'y');

%% Fill in the new regions

for i=1:floor((N2-n2)/4)
  line(x(n2+4*i)*ones(1, 100), linspace(M_0, plotH), 'Color', 'y');
end

for i=1:floor((N2-n2)/4)
  line(x(n2+4*i)*ones(1, 100), linspace(-M_0, -plotH), 'Color', 'y');
end

%% Plot some possible functions g

g1 = x.^2 - M_0/4;
plot(x, g1, 'b');
g2 = 0.125*x.^4 - M_0/1.5;
plot(x, g2, 'r');
g3 = 2.5*sqrt(x.^2 + 1) - M_0/8;
plot(x, g3, 'k');

hold off
