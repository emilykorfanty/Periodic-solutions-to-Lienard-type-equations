clear

ParameterTest;

%% Setting up plot parameters
plotW = 100;
plotH = k_1*plotW + k_2;
xf = plotW;
xi = 0;
N = 1000;
x=linspace(xi, xf, N);

%Define f  as a sinusoidal here, given M and m:
y0 = (M-m)/2;
A = M-y0;  
t = linspace(0,T, 10*N);
f = A*sin(2*pi*t/T) + y0;

% Calculate ave f
ave_f = mean(f);

%% Understanding condition (B_2)
% We know that this condition is satisfied if: whenever |x| is outside of
% the region R_1, R_2, g(x) \neq cx for all c \in [m, M]. So, we plot the
% regions where g is NOT allowed to pass though.

% First, let's plot the cone defined by y = mx and y=Mx

dx = (xf-xi)/N;  %spacing of linspace
N1 = floor((R_1-xi)/dx);
N2 = ceil((R_2-xi)/dx);

Mline = M*x;
mline = m*x;

plot(x, 0*x, 'k');
hold on
axis([0 plotW -plotH plotH])
plot(x(1:N1), mline(1:N1), 'g');
plot(x(N2:N), mline(N2:N), 'g');
plot(x(1:N1), Mline(1:N1), 'g');
plot(x(N2:N), Mline(N2:N), 'g');


% Let's fill in the regions where g(x) can't pass through! Using a bunch of
% vertical lines.

line(x(N1)*ones(1, 100), linspace(mline(N1), Mline(N1)), 'Color', 'g');
line(x(N2)*ones(1, 100), linspace(mline(N2), Mline(N2)), 'Color', 'g');

for i=1:floor(N1/4)
  line(x(4*i)*ones(1, 100), linspace(mline(4*i), Mline(4*i)), 'Color', 'g');
end

for i=1:floor((N-N2)/4)
  line(x(N2+4*i - 1)*ones(1, 100), linspace(mline(N2+4*i - 1), Mline(N2+4*i - 1)), 'Color', 'g');
end

%% Adding condition (B_1): g must be sublinear

K = zeros(1, N);
for i=1:N
    K(i) = k_1*x(i) + k_2;
end

plot(x, K, 'c');
plot(x, -K, 'c');

% Fill in the not-allowed region!
i=1;

while K(4*i)< plotH && 4*i<N
    line(x(4*i)*ones(1, 100), linspace(K(4*i), plotH), 'Color', 'c');
     line(x(4*i)*ones(1, 100), linspace(-K(4*i), -plotH), 'Color', 'c');
    i = i+1;
end

%% Plot R

%plot(x(1:N1), R*ones(1, N1), 'r'); %% want blue > red

%% How close are we?

% Let's see how much room g has between 0 and R_1; g must be greater than
% g_min.
%plot(x(1:N1), g_min*ones(1, N1), 'b');

%for i=1:floor(N1/8)
 % line(x(8*i)*ones(1, 100), linspace(-plotH, g_min), 'Color', 'b');
%end

%% Understanding condition (B_5) 

% First, we need g(x)> ave f x in (0, R_1)

plot(x(1:N1), ave_f*x(1:N1), 'y');
line(x(N1)*ones(1, 100), linspace(-plotH, ave_f*x(N1)), 'Color', 'y');
for i=1:floor(N1/8);
  line(x(8*i)*ones(1, 100), linspace(-plotH, ave_f*x(8*i)), 'Color', 'y');
end

% Lastly, we need g(x) < ave f x for x>R_2

plot(x(N2:N), ave_f*x(N2:N), 'y');
line(x(N2)*ones(1, 100), linspace(ave_f*x(N2), plotH), 'Color', 'y');

for i=1:floor((N-N2)/8);
  line(x(N2 + 8*i)*ones(1, 100), linspace(ave_f*x(N2 + 8*i), plotH), 'Color', 'y');
end

%% Can we actually find an example??

V = (R_2 + R_1)/2;

g1 = zeros(1,N);
g2 = zeros(1, N);
g3 = zeros(1, N);
g4 = zeros(1, N);

for i = 1:N
    g1(i) = (-8*(x(i) - V)/((x(i)-V)^2 + 4))*((((x(i)-V)/40)^2)+1);
    g2(i) = (-6*(x(i) - V)/((x(i)-V)^2 + 2))*((((x(i)-V)/20)^2)+1);
    g3(i) = (-3*(x(i) - V)/((x(i)-V)^2 + 0.5))*((((x(i)-V)/10)^2)+1);
end


plot(x, g1, 'k');
plot(x, g2, 'b');
plot(x, g3, 'r');

gplot1 = plot(x, g1, 'k');
gplot2 = plot(x, g2, 'b');
gplot3 = plot(x, g3, 'r');

xlabel('x');
ylabel('g(x)');
hleglines = [gplot1(1) gplot2(1) gplot3(1)];
legend(hleglines,'a = 8, b = 4, c = 40','a = 6, b = 2, c = 20','a = 3, b = 0.5, c = 10');

%% Plot R/R_1 for reference

plot(x(1:N1), (R/R_1)*ones(1,N1), 'm')

hold off
