%% Part A: Understanding the Hill equation
%% 1.
%% 1(a)
fig1 = figure(1);
S = 0:0.1:100; % mM
Vmax = 5; % mM/s
K = 20.0; % mM

h1 = 1;
v1 = (Vmax.*(S.^h1))./((K.^h1) + (S.^h1));
plot(S, v1)
hold on

h2 = 2;
v2 = (Vmax.*(S.^h2))./((K.^h2) + (S.^h2));
plot(S, v2)
hold on

h3 = 10;
v3 = (Vmax.*(S.^h3))./((K.^h3) + (S.^h3));
plot(S, v3)

xlabel('Substrate Concentration (mM)')
ylabel('Reaction Velocity (mM/s)')
legend('h=1', 'h=2', 'h=10', 'Location','southeast')
hold off

%% 1(b)
fig2 = figure(2);
S = 0:0.1:100; % mM
Vmax = 5; % mM/s
h = 2;

K1 = 10.0; % mM
v1 = (Vmax.*(S.^h))./((K1.^h) + (S.^h));
plot(S, v1)
hold on

K2 = 20.0;
v2 = (Vmax.*(S.^h))./((K2.^h) + (S.^h));
plot(S, v2)
hold on

K3 = 40.0;
v3 = (Vmax.*(S.^h))./((K3.^h) + (S.^h));
plot(S, v3)

xlabel('Substrate Concentration (mM)')
ylabel('Reaction Velocity (mM/s)')
legend('K_{1/2}= 10mM', 'K_{1/2}= 20mM', 'K_{1/2}= 40mM', 'Location','southeast')
hold off

%% 1(c)
fig3 = figure(3);
S = 0:0.1:100; % mM
K = 20.0; % mM
h = 2;

Vmax1 = 2; % mM/s
v1 = (Vmax1.*(S.^h))./((K.^h) + (S.^h));
plot(S, v1)
hold on

Vmax2 = 5; % mM/s
v2 = (Vmax2.*(S.^h))./((K.^h) + (S.^h));
plot(S, v2)
hold on

Vmax3 = 10; % mM/s
v3 = (Vmax3.*(S.^h))./((K.^h) + (S.^h));
plot(S, v3)

xlabel('Substrate Concentration (mM)')
ylabel('Reaction Velocity (mM/s)')
legend('V_{max}= 2mM/s', 'V_{max}= 5mM/s', 'V_{max}= 10mM/s', 'Location','northwest')
hold off

%% 2.
%% 2(b)
clear
fig4 = figure(4);
S = 0:0.1:100; % mM
Vmax = 20; % mM/s
K = 50; % mM
h = 4;
v = (Vmax.*(S.^h))./((K.^h) + (S.^h));
plot(S, v)
xlabel('Substrate Concentration (mM)')
ylabel('Reaction Velocity (mM/s)')

%% 3.
%% 3(a)
clear
fig5 = figure(5);
Vmax = 20; % mM/s
K = 50; % mM
h = 4;

deltaT = 0.01; % s
maxT = 10; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;
S0 = 100; % initial concentration (mM)
S = zeros(size(t));
S(1) = S0;
for i = 1:numiterations
    dSdt = -(Vmax*(S(i)^h))/((K^h) + (S(i)^h));
    S(i+1) = S(i) + dSdt*deltaT;
end
plot(t, S)
xlabel('Time (s)')
ylabel('Substrate Concentration (mM)')
