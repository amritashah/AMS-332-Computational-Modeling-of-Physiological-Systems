%% Part B: An autoregulatory gene
%% 1.
clear

u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM

deltaT = 0.01; % s
maxT = 20; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

%% 1(a)
fig6 = figure(6);
(Xprot(i)^2) / ((K^2)+(Xprot(i)^2)) - (Xr*Xrna(i));
dXprotdt = (w*Xrna(i)) - (Xp*Xprot(i));
Xrna_0 = 0.5; % initial concentration (mM)
Xprot_0 = 0.5; % initial concentration (mM)
Xrna = zeros(size(t));
Xprot = zeros(size(t));
Xrna(1) = Xrna_0;
Xprot(1) = Xprot_0;

for i = 1:numiterations
    dXrnadt = (u*a(i) + dXrnadt*deltaT;
    Xprot(i+1) = Xprot(i) + dXprotdt*deltaT;
end
plot(t, Xrna)
hold on
plot(t, Xprot)
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend('[X_{RNA}]', '[X_{protein}]', 'Location', 'southeast')
hold off

%% 1(b)
fig7 = figure(7);
Xrna_0 = 0; % initial concentration (mM)
Xprot_0 = 0.2; % initial concentration (mM)
Xrna = zeros(size(t));
Xprot = zeros(size(t));
Xrna(1) = Xrna_0;
Xprot(1) = Xprot_0;

for i = 1:numiterations
    dXrnadt = ((u*(Xprot(i)^2)) / ((K^2)+(Xprot(i)^2))) - (Xr*Xrna(i));
    dXprotdt = (w*Xrna(i)) - (Xp*Xprot(i));
    Xrna(i+1) = Xrna(i) + dXrnadt*deltaT;
    Xprot(i+1) = Xprot(i) + dXprotdt*deltaT;
end
plot(t, Xrna)
hold on
plot(t, Xprot)
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend('[X_{RNA}]', '[X_{protein}]')
hold off

%% 1(c)
fig8 = figure(8);
Xrna_0 = 0; % initial concentration (mM)
Xprot_0 = 0.5; % initial concentration (mM)
Xrna = zeros(size(t));
Xprot = zeros(size(t));
Xrna(1) = Xrna_0;
Xprot(1) = Xprot_0;

for i = 1:numiterations
    dXrnadt = ((u*(Xprot(i)^2)) / ((K^2)+(Xprot(i)^2))) - (Xr*Xrna(i));
    dXprotdt = (w*Xrna(i)) - (Xp*Xprot(i));
    Xrna(i+1) = Xrna(i) + dXrnadt*deltaT;
    Xprot(i+1) = Xprot(i) + dXprotdt*deltaT;
end
plot(t, Xrna)
hold on
plot(t, Xprot)
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend('[X_{RNA}]', '[X_{protein}]', 'Location', 'southeast')
hold off

%% 1(d)
fig9 = figure(9);
Xrna_0 = 0.2; % initial concentration (mM)
Xprot_0 = 0; % initial concentration (mM)
Xrna = zeros(size(t));
Xprot = zeros(size(t));
Xrna(1) = Xrna_0;
Xprot(1) = Xprot_0;

for i = 1:numiterations
    dXrnadt = ((u*(Xprot(i)^2)) / ((K^2)+(Xprot(i)^2))) - (Xr*Xrna(i));
    dXprotdt = (w*Xrna(i)) - (Xp*Xprot(i));
    Xrna(i+1) = Xrna(i) + dXrnadt*deltaT;
    Xprot(i+1) = Xprot(i) + dXprotdt*deltaT;
end
plot(t, Xrna)
hold on
plot(t, Xprot)
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend('[X_{RNA}]', '[X_{protein}]')
hold off

%% 1(e)
fig10 = figure(10);
Xrna_0 = 0.5; % initial concentration (mM)
Xprot_0 = 0; % initial concentration (mM)
Xrna = zeros(size(t));
Xprot = zeros(size(t));
Xrna(1) = Xrna_0;
Xprot(1) = Xprot_0;

for i = 1:numiterations
    dXrnadt = ((u*(Xprot(i)^2)) / ((K^2)+(Xprot(i)^2))) - (Xr*Xrna(i));
    dXprotdt = (w*Xrna(i)) - (Xp*Xprot(i));
    Xrna(i+1) = Xrna(i) + dXrnadt*deltaT;
    Xprot(i+1) = Xprot(i) + dXprotdt*deltaT;
end
plot(t, Xrna)
hold on
plot(t, Xprot)
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend('[X_{RNA}]', '[X_{protein}]', 'Location', 'southeast')
hold off

%% 2.
%% 2(a-b)
clear
close all
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
deltaT = 0.01; % s
maxT = 20; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

fig11 = figure(11);
for Xrna_0 = 0:0.2:1.4 % initial concentration (mM)
    for Xprot_0 = 0:0.2:1.4 % initial concentration (mM)
        Xrna = zeros(size(t));
        Xprot = zeros(size(t));
        Xrna(1) = Xrna_0;
        Xprot(1) = Xprot_0;
        for i = 1:numiterations
            dXrnadt = ((u*(Xprot(i)^2)) / ((K^2)+(Xprot(i)^2))) - (Xr*Xrna(i));
            dXprotdt = (w*Xrna(i)) - (Xp*Xprot(i));
            Xrna(i+1) = Xrna(i) + dXrnadt*deltaT;
            Xprot(i+1) = Xprot(i) + dXprotdt*deltaT;
        end
        plot(Xprot, Xrna)
        hold on
    end
end
xlabel('[X_{protein}] (mM)')
ylabel('[X_{RNA}] (mM)')
hold off
