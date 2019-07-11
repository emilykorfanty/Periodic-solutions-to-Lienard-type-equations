clear

%% Understanding conditions for Theorem 3.1
% This is a program that, given a continuous function f defined on [0,T], 
% shows the possible regions where g(x) can be defined to satisfy (A_1) to
% (A_4); the approach first specifies R_1 and R_2 along with the min and
% max of f, and looks at the conditions to determine what g can be.

%% Setting up parameters
T = 1;
R_1 = 1/4;
R_2 = 2;

plotW = 100;
plotH = 100;
xf = plotW;
xi = -plotW;
N = 1000;
x=linspace(xi, xf, N);

%Define f  as a sinusoidal here, given M and m:
M =  -0.55;
m = -0.85;
y0 = (M+m)/2;
A = abs(M-y0);  
t = linspace(0,T, 10*N);
f = A*sin(2*pi*t/T) + y0;

% Calculate ave f
ave_f = mean(f);

maxk = sqrt(12)/(T^2) - M;
k_1 = maxk - 0.001;
k_2 = 1;

%% Understanding condition (A_2)
% We know that this condition is satisfied if: whenever |x| is outside of
% the region R_1, R_2, g(x) \neq cx for all c \in [m, M]. So, we plot the
% regions where g is NOT allowed to pass though.

% First, let's plot the cone defined by y = mx and y=Mx

dx = (xf-xi)/N;  %spacing of linspace
N1 = floor((R_1-xi)/dx);
N2 = ceil((R_2-xi)/dx);
N3 = floor((-R_2-xi)/dx);
N4 = ceil((-R_1-xi)/dx);

Mline = M*x;
mline = m*x;

plot(x(1:N3), mline(1:N3), 'g');
axis([-plotW plotW -plotH plotH])
hold on
plot(x(N4:N1), mline(N4:N1), 'g');
plot(x(N2:N), mline(N2:N), 'g');
plot(x(1:N3), Mline(1:N3), 'g');
plot(x(N4:N1), Mline(N4:N1), 'g');
plot(x(N2:N), Mline(N2:N), 'g');
plot(x, 0*x, 'k');

y1 = linspace(mline(N1), Mline(N1));
line(R_1*ones(1, 100),y1, 'Color','green');
y2 = linspace(mline(N2), Mline(N2));
line(R_2*ones(1, 100),y2, 'Color','green');
y3 = linspace(mline(N3), Mline(N3));
line(-R_2*ones(1, 100),y3, 'Color','green');
y4 = linspace(mline(N4), Mline(N4));
line(-R_1*ones(1, 100),y4, 'Color','green');

% Let's fill in the regions where g(x) can't pass through! Using a bunch of
% vertical lines.

for i=1:floor(N3/4)
  line(x(4*i)*ones(1, 100), linspace(mline(4*i), Mline(4*i)), 'Color', 'g');
end

for i=1:floor((N1-N4)/4)
  line(x(N4+4*i)*ones(1, 100), linspace(mline(N4+4*i), Mline(N4+4*i)), 'Color', 'g');
end

for i=1:floor((N-N2)/4)-1
  line(x(N2+4*i)*ones(1, 100), linspace(mline(N2+4*i), Mline(N2+4*i)), 'Color', 'g');
end

%% Adding condition (A_1): g must be sublinear

K1 = zeros(1, N);
for i=1:N
    K1(i) = k_1*abs(x(i)) + k_2;
end
K2 = -K1;

plot(x, K1, 'c');
plot(x, K2, 'c');

% Fill in the not-allowed region!

xk1 = -(plotH - k_2)/k_1;
xk2 = (plotH - k_2)/k_1;

Nk1 = floor((xk1-xi)/dx);
Nk2 = ceil((xk2-xi)/dx);

for i=1:floor(Nk2/8)
  line(x(Nk1+4*i)*ones(1, 100), linspace(K1(Nk1+4*i), plotH), 'Color', 'c');
end

for i=1:floor(Nk2/8)
  line(x(Nk1+4*i)*ones(1, 100), linspace(-plotH, K2(Nk1+4*i)), 'Color', 'c');
end

%% Investigating condition (A_4)

% Let's plot the bound on the end behaviour of g; we want g< ave f *x + b
% for some constant b.

AVE = ave_f*x;
plot(x, AVE, 'y');

R = x(7*N/8);

line(R*ones(1, 100), linspace(AVE(7*N/8), plotH, 100), 'Color', 'y');

for i = 1 : floor(N/64)
   line(x(7*N/8  + 8*i)*ones(1, 100), linspace(AVE(7*N/8  + 8*i), plotH, 100), 'Color', 'y');
end

line(-R*ones(1, 100), linspace(AVE(N/8), plotH, 100), 'Color', 'y');

for i = 1 : floor(N/64)
   line(x(8*i)*ones(1, 100), linspace(AVE(8*i), plotH, 100), 'Color', 'y');
end


%% g(0) \neq 0
plot(0,0, 'mo');

%% Let's plot some possible g!


%a = 1;
%b = (M+k_2)/2;
%c = 0.05;
%g1 = -b*sqrt((x.^2/a^2) + 1) + c;
%gplot1 = plot(x, g1, 'b');

%a = 1;
%b = 2;
%c = 1.75;
%g2 = -b*sqrt((x.^2/a^2) + 1) + c;
%gplot2 = plot(x, g2, 'k');

%a = 1;
%b = 2.75;
%c = 3;
%g3 = -b*sqrt((x.^2/a^2) + 1) + c;
%gplot3 = plot(x, g3, 'r');

% choose legend data
%hleglines = [gplot1(1) gplot2(1) gplot3(1)];
% create the legend
%legend(hleglines,'b = 0.5, c = 0.05','b = 2, c = 1.75','b = 2.75, c = 3');
%title('max(f) = 0.5, min(f) = -0.25, k_1 = 2.9631, k_2 = 0.5, R_1 = 0.25, R_2 = 0.75, R = 1');
xlabel('x');
%ylabel('g(x) = -b(x^2 +1)^{1/2} + c');

g = zeros(1, N);
for j=1:N
    g(j) = 2*(exp(-1/x(j)^2)+0.1)*(1-x(j))*exp(4*x(j)+1)*exp(0.25*sin(4*x(j)+1))/(exp(4*x(j)+1)+1)-x(j)/2;
end

plot(x, g, 'b');

hold off
