%% Part A: Hodgkin-Huxley model
%% 1. Main simulation
clear
close all

GNa = 400; % nS
GK = 200; % nS
GL = 2; % nS
ENa = 99; % mV
EK = -85; % mV
VL = -65; % mV
C = 2; % pF

maxT = 200; % ms
deltaT = 0.01; % ms
Tswap = 40;
numiterations = maxT/deltaT;
time = 0:deltaT:maxT;

V = zeros(size(time));
m = zeros(size(time));
h = zeros(size(time));
n = zeros(size(time));

V(1) = VL;

Ie = zeros(size(time));
Ie((Tswap/deltaT)+1:end) = 200;

Im = zeros(size(time));
Na_cond = zeros(size(time));
K_cond = zeros(size(time));

Im(1) = GL*(V(1)-VL) + GNa*(m(1)^3)*h(1)*(V(1)-ENa) + GK*(n(1)^4)*(V(1)-EK);
Na_cond(1) = GNa*(m(1)^3)*h(1);
K_cond(1) = GK*(n(1)^4);

% transition rates for V(1)
a_m = (0.1*(V(1)+40)) / (1- exp(-0.1*(V(1)+40)));
a_h = 0.07 * exp(-0.05*(V(1)+65));
a_n = (0.01*(V(1)+55)) / (1- exp(-0.1*(V(1)+55)));
b_m = 4 * exp(-0.0556*(V(1)+65));
b_h = 1 / (1+ exp(-0.1*(V(1)+35)));
b_n = 0.125 * exp(-0.0125*(V(1)+65));

m_inf = a_m / (a_m + b_m);
h_inf = a_h / (a_h + b_h);
n_inf = a_n / (a_n + b_n);

m(1) = m_inf;
h(1) = h_inf;
n(1) = n_inf;

for i = 1:numiterations

    % transition rates for this V
    a_m = (0.1*(V(i)+40)) / (1- exp(-0.1*(V(i)+40)));
    a_h = 0.07 * exp(-0.05*(V(i)+65));
    a_n = (0.01*(V(i)+55)) / (1- exp(-0.1*(V(i)+55)));
    b_m = 4 * exp(-0.0556*(V(i)+65));
    b_h = 1 / (1+ exp(-0.1*(V(i)+35)));
    b_n = 0.125 * exp(-0.0125*(V(i)+65));

    dm_dt = a_m*(1-m(i)) - b_m*(m(i));
    dh_dt = a_h*(1-h(i)) - b_h*(h(i));
    dn_dt = a_n*(1-n(i)) - b_n*(n(i));

    m(i+1) = m(i) + dm_dt*deltaT;
    h(i+1) = h(i) + dh_dt*deltaT;
    n(i+1) = n(i) + dn_dt*deltaT;

    dV_dt = ( -GL*(V(i)-VL) - GNa*(m(i)^3)*h(i)*(V(i)-ENa) - GK*(n(i)^4)*(V(i)-EK) + Ie(i)) / C;
    V(i+1) = V(i) + dV_dt*deltaT;
    
    Im(i+1) = GL*(V(i+1)-VL) + GNa*(m(i+1)^3)*h(i+1)*(V(i+1)-ENa) + GK*(n(i+1)^4)*(V(i+1)-EK);
    Na_cond(i+1) = GNa*(m(i+1)^3)*h(i+1);
    K_cond(i+1) = GK*(n(i+1)^4);

end

fig1 = figure(1);
subplot(5,1,1)
plot(time, V)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Membrane Potential vs. Time')

subplot(5,1,2)
plot(time, Im)
xlabel('Time (ms)')
ylabel('Total Membrane Current (pA)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Total Membrane Current vs. Time')

subplot(5,1,3)
plot(time,Na_cond)
xlabel('Time (ms)')
ylabel('Na+ Conductance (nS)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Sodium Conductance vs. Time')

subplot(5,1,4)
plot(time,K_cond)
xlabel('Time (ms)')
ylabel('K+ Conductance (nS)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Potassium Conductance vs. Time')

subplot(5,1,5)
plot(time, Ie)
xlabel('Time (ms)')
ylabel('Input Current (pA)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Current I_{e} vs. Time')

%% 2. Spike detection
clear
close all

GNa = 400; % nS
GK = 200; % nS
GL = 2; % nS
ENa = 99; % mV
EK = -85; % mV
VL = -65; % mV
C = 2; % pF
Vspk = 0; % mV

maxT = 200; % ms
deltaT = 0.01; % ms
Tswap = 40;
numiterations = maxT/deltaT;
time = 0:deltaT:maxT;

spikes = 0;
Tspikes = [];

V = zeros(size(time));
m = zeros(size(time));
h = zeros(size(time));
n = zeros(size(time));

V(1) = VL;

Ie = zeros(size(time));
Ie((Tswap/deltaT)+1:end) = 200;

Im = zeros(size(time));
Na_cond = zeros(size(time));
K_cond = zeros(size(time));

Im(1) = GL*(V(1)-VL) + GNa*(m(1)^3)*h(1)*(V(1)-ENa) + GK*(n(1)^4)*(V(1)-EK);
Na_cond(1) = GNa*(m(1)^3)*h(1);
K_cond(1) = GK*(n(1)^4);

% transition rates for V(1)
a_m = (0.1*(V(1)+40)) / (1- exp(-0.1*(V(1)+40)));
a_h = 0.07 * exp(-0.05*(V(1)+65));
a_n = (0.01*(V(1)+55)) / (1- exp(-0.1*(V(1)+55)));
b_m = 4 * exp(-0.0556*(V(1)+65));
b_h = 1 / (1+ exp(-0.1*(V(1)+35)));
b_n = 0.125 * exp(-0.0125*(V(1)+65));

m_inf = a_m / (a_m + b_m);
h_inf = a_h / (a_h + b_h);
n_inf = a_n / (a_n + b_n);

m(1) = m_inf;
h(1) = h_inf;
n(1) = n_inf;

for i = 1:numiterations

    % transition rates for this V
    a_m = (0.1*(V(i)+40)) / (1- exp(-0.1*(V(i)+40)));
    a_h = 0.07 * exp(-0.05*(V(i)+65));
    a_n = (0.01*(V(i)+55)) / (1- exp(-0.1*(V(i)+55)));
    b_m = 4 * exp(-0.0556*(V(i)+65));
    b_h = 1 / (1+ exp(-0.1*(V(i)+35)));
    b_n = 0.125 * exp(-0.0125*(V(i)+65));

    dm_dt = a_m*(1-m(i)) - b_m*(m(i));
    dh_dt = a_h*(1-h(i)) - b_h*(h(i));
    dn_dt = a_n*(1-n(i)) - b_n*(n(i));

    m(i+1) = m(i) + dm_dt*deltaT;
    h(i+1) = h(i) + dh_dt*deltaT;
    n(i+1) = n(i) + dn_dt*deltaT;

    dV_dt = ( -GL*(V(i)-VL) - GNa*(m(i)^3)*h(i)*(V(i)-ENa) - GK*(n(i)^4)*(V(i)-EK) + Ie(i)) / C;
    V(i+1) = V(i) + dV_dt*deltaT;
    if V(i+1) > Vspk && V(i) < Vspk
        spikes = spikes+1;
        Tspikes(end+1) = i/100;
    end

    Im(i+1) = GL*(V(i+1)-VL) + GNa*(m(i+1)^3)*h(i+1)*(V(i+1)-ENa) + GK*(n(i+1)^4)*(V(i+1)-EK);
    Na_cond(i+1) = GNa*(m(i+1)^3)*h(i+1);
    K_cond(i+1) = GK*(n(i+1)^4);

end

fig2 = figure(2);
sgtitle('Second Spike')

subplot(5,1,1)
plot(time, V)
xlim([47.7 57.7])
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Membrane Potential vs. Time')
hold on
plot(time, zeros(size(V)), ':r')
hold off

subplot(5,1,2)
plot(time, Im)
xlim([47.7 57.7])
xlabel('Time (ms)')
ylabel('Total Membrane Current (pA)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Total Membrane Current vs. Time')

subplot(5,1,3)
plot(time,Na_cond)
xlim([47.7 57.7])
xlabel('Time (ms)')
ylabel('Na+ Conductance (nS)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Sodium Conductance vs. Time')

subplot(5,1,4)
plot(time,K_cond)
xlim([47.7 57.7])
xlabel('Time (ms)')
ylabel('K+ Conductance (nS)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Potassium Conductance vs. Time')

subplot(5,1,5)
plot(time, Ie)
xlim([47.7 57.7])
xlabel('Time (ms)')
ylabel('Input Current (pA)', 'Rotation', 0, 'HorizontalAlignment','right')
title('Current I_{e} vs. Time')

spikes
Tspikes

%% 3. Action potential threshold
clear
close all

GNa = 400; % nS
GK = 200; % nS
GL = 2; % nS
ENa = 99; % mV
EK = -85; % mV
VL = -65; % mV
C = 2; % pF

maxT = 200; % ms
deltaT = 0.01; % ms
Tswap = 40;
numiterations = maxT/deltaT;
time = 0:deltaT:maxT;

V = zeros(size(time));
m = zeros(size(time));
h = zeros(size(time));
n = zeros(size(time));

V(1) = VL;

Ie = zeros(size(time));
Ie((Tswap/deltaT)+1:end) = 18.43;

% transition rates for V(1)
a_m = (0.1*(V(1)+40)) / (1- exp(-0.1*(V(1)+40)));
a_h = 0.07 * exp(-0.05*(V(1)+65));
a_n = (0.01*(V(1)+55)) / (1- exp(-0.1*(V(1)+55)));
b_m = 4 * exp(-0.0556*(V(1)+65));
b_h = 1 / (1+ exp(-0.1*(V(1)+35)));
b_n = 0.125 * exp(-0.0125*(V(1)+65));

m_inf = a_m / (a_m + b_m);
h_inf = a_h / (a_h + b_h);
n_inf = a_n / (a_n + b_n);

m(1) = m_inf;
h(1) = h_inf;
n(1) = n_inf;

for i = 1:numiterations

    % transition rates for this V
    a_m = (0.1*(V(i)+40)) / (1- exp(-0.1*(V(i)+40)));
    a_h = 0.07 * exp(-0.05*(V(i)+65));
    a_n = (0.01*(V(i)+55)) / (1- exp(-0.1*(V(i)+55)));
    b_m = 4 * exp(-0.0556*(V(i)+65));
    b_h = 1 / (1+ exp(-0.1*(V(i)+35)));
    b_n = 0.125 * exp(-0.0125*(V(i)+65));

    dm_dt = a_m*(1-m(i)) - b_m*(m(i));
    dh_dt = a_h*(1-h(i)) - b_h*(h(i));
    dn_dt = a_n*(1-n(i)) - b_n*(n(i));

    m(i+1) = m(i) + dm_dt*deltaT;
    h(i+1) = h(i) + dh_dt*deltaT;
    n(i+1) = n(i) + dn_dt*deltaT;

    dV_dt = ( -GL*(V(i)-VL) - GNa*(m(i)^3)*h(i)*(V(i)-ENa) - GK*(n(i)^4)*(V(i)-EK) + Ie(i)) / C;
    V(i+1) = V(i) + dV_dt*deltaT;
    
end

fig3 = figure(3);
plot(time, V)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
title('Membrane Potential vs. Time')

%% 4. Rheobase current
clear
close all

GNa = 400; % nS
GK = 200; % nS
GL = 2; % nS
ENa = 99; % mV
EK = -85; % mV
VL = -65; % mV
C = 2; % pF

maxT = 1000; % ms
deltaT = 0.01; % ms
Tswap = 40;
numiterations = maxT/deltaT;
time = 0:deltaT:maxT;

V = zeros(size(time));
m = zeros(size(time));
h = zeros(size(time));
n = zeros(size(time));

V(1) = VL;

Ie = zeros(size(time));
Ie((Tswap/deltaT)+1:end) = 108.62;

% transition rates for V(1)
a_m = (0.1*(V(1)+40)) / (1- exp(-0.1*(V(1)+40)));
a_h = 0.07 * exp(-0.05*(V(1)+65));
a_n = (0.01*(V(1)+55)) / (1- exp(-0.1*(V(1)+55)));
b_m = 4 * exp(-0.0556*(V(1)+65));
b_h = 1 / (1+ exp(-0.1*(V(1)+35)));
b_n = 0.125 * exp(-0.0125*(V(1)+65));

m_inf = a_m / (a_m + b_m);
h_inf = a_h / (a_h + b_h);
n_inf = a_n / (a_n + b_n);

m(1) = m_inf;
h(1) = h_inf;
n(1) = n_inf;

for i = 1:numiterations

    % transition rates for this V
    a_m = (0.1*(V(i)+40)) / (1- exp(-0.1*(V(i)+40)));
    a_h = 0.07 * exp(-0.05*(V(i)+65));
    a_n = (0.01*(V(i)+55)) / (1- exp(-0.1*(V(i)+55)));
    b_m = 4 * exp(-0.0556*(V(i)+65));
    b_h = 1 / (1+ exp(-0.1*(V(i)+35)));
    b_n = 0.125 * exp(-0.0125*(V(i)+65));

    dm_dt = a_m*(1-m(i)) - b_m*(m(i));
    dh_dt = a_h*(1-h(i)) - b_h*(h(i));
    dn_dt = a_n*(1-n(i)) - b_n*(n(i));

    m(i+1) = m(i) + dm_dt*deltaT;
    h(i+1) = h(i) + dh_dt*deltaT;
    n(i+1) = n(i) + dn_dt*deltaT;

    dV_dt = ( -GL*(V(i)-VL) - GNa*(m(i)^3)*h(i)*(V(i)-ENa) - GK*(n(i)^4)*(V(i)-EK) + Ie(i)) / C;
    V(i+1) = V(i) + dV_dt*deltaT;
    
end

fig4 = figure(4);
plot(time, V)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
title('Membrane Potential vs. Time')

%% 5. f-I curve
clear
close all

runs = 100:0.5:125;
fr = zeros(1, length(runs));

for r = 1:length(runs)

    GNa = 400; % nS
    GK = 200; % nS
    GL = 2; % nS
    ENa = 99; % mV
    EK = -85; % mV
    VL = -65; % mV
    C = 2; % pF
    Vspk = 0; % mV
    
    maxT = 1000; % ms
    deltaT = 0.01; % ms
    Tswap = 40;
    numiterations = maxT/deltaT;
    time = 0:deltaT:maxT;
    
    spikes = 0;
    Tspikes = [];
    
    V = zeros(size(time));
    m = zeros(size(time));
    h = zeros(size(time));
    n = zeros(size(time));
    
    V(1) = VL;
    
    Ie = zeros(size(time));
    Ie((Tswap/deltaT)+1:end) = runs(r);
    
    % transition rates for V(1)
    a_m = (0.1*(V(1)+40)) / (1- exp(-0.1*(V(1)+40)));
    a_h = 0.07 * exp(-0.05*(V(1)+65));
    a_n = (0.01*(V(1)+55)) / (1- exp(-0.1*(V(1)+55)));
    b_m = 4 * exp(-0.0556*(V(1)+65));
    b_h = 1 / (1+ exp(-0.1*(V(1)+35)));
    b_n = 0.125 * exp(-0.0125*(V(1)+65));
    
    m_inf = a_m / (a_m + b_m);
    h_inf = a_h / (a_h + b_h);
    n_inf = a_n / (a_n + b_n);
    
    m(1) = m_inf;
    h(1) = h_inf;
    n(1) = n_inf;
    
    for i = 1:numiterations
    
        % transition rates for this V
        a_m = (0.1*(V(i)+40)) / (1- exp(-0.1*(V(i)+40)));
        a_h = 0.07 * exp(-0.05*(V(i)+65));
        a_n = (0.01*(V(i)+55)) / (1- exp(-0.1*(V(i)+55)));
        b_m = 4 * exp(-0.0556*(V(i)+65));
        b_h = 1 / (1+ exp(-0.1*(V(i)+35)));
        b_n = 0.125 * exp(-0.0125*(V(i)+65));
    
        dm_dt = a_m*(1-m(i)) - b_m*(m(i));
        dh_dt = a_h*(1-h(i)) - b_h*(h(i));
        dn_dt = a_n*(1-n(i)) - b_n*(n(i));
    
        m(i+1) = m(i) + dm_dt*deltaT;
        h(i+1) = h(i) + dh_dt*deltaT;
        n(i+1) = n(i) + dn_dt*deltaT;
    
        dV_dt = ( -GL*(V(i)-VL) - GNa*(m(i)^3)*h(i)*(V(i)-ENa) - GK*(n(i)^4)*(V(i)-EK) + Ie(i)) / C;
        V(i+1) = V(i) + dV_dt*deltaT;
        if V(i+1) > Vspk && V(i) < Vspk
            spikes = spikes+1;
            Tspikes(end+1) = i/100;
        end 
    end
    
    if spikes < 15
        fr(r) = 0;
    else     
        fr(r) = 1/(mean(diff(Tspikes)));
    end
end
    
fig5 = figure(5);
plot(runs, fr*1000)
xlabel('Input Current (pA)')
ylabel('Firing Rate (spks/s)')
title('f-I Curve')
