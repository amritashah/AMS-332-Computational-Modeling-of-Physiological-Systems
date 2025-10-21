%% Part B: Leaky Integrate-and-Fire neuron model
%% 1. Main simulation
clear
close all

GL = 50/1000; % uS
VL = -65; % mV
C = 1; % nF
Vspk = -45; % mV
Vr = -65; % mV
Tau = 2; % ms

maxT = 200; % ms
deltaT = 0.1; % ms
numiterations = maxT/deltaT;
time = 0:deltaT:maxT;

spikes = 0;
Tspikes = [];
V = zeros(size(time));
V(1) = VL;

Ie = 1.1*ones(size(time));

for i = 1:numiterations
    if V(i+1) == Vr
        continue
    end
    dV_dt = ((-GL*(V(i)-VL)) + Ie(i)) / C;
    V(i+1) = V(i) + dV_dt*deltaT;
    if V(i+1) > Vspk && V(i) < Vspk
        spikes = spikes+1;
        Tspikes(end+1) = i/10;
        V((i+1):(i+1+(Tau/deltaT))) = Vr;
    end 
end

fig1 = figure(1);
subplot(2,1,1)
plot(time, V(1:numiterations+1))
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
title('Membrane Potential vs. Time')

subplot(2,1,2)
plot(time, Ie)
xlabel('Time (ms)')
ylabel('Input Current (nA)')
title('Current I_{e} vs. Time')

%% 2. f-I curve
clear
close all

runs = 0.5:0.05:3;
fr = zeros(1, length(runs));

for r = 1:length(runs)
    GL = 50/1000; % uS
    VL = -65; % mV
    C = 1; % nF
    Vspk = -45; % mV
    Vr = -65; % mV
    Tau = 2; % ms
    
    maxT = 1000; % ms
    deltaT = 0.1; % ms
    numiterations = maxT/deltaT;
    time = 0:deltaT:maxT;
    
    spikes = 0;
    Tspikes = [];
    V = zeros(size(time));
    V(1) = VL;
    
    Ie = runs(r);
    
    for i = 1:numiterations
        if V(i+1) == Vr
            continue
        end
        dV_dt = ((-GL*(V(i)-VL)) + Ie) / C;
        V(i+1) = V(i) + dV_dt*deltaT;
        if V(i+1) > Vspk && V(i) < Vspk
            spikes = spikes+1;
            Tspikes(end+1) = i/10;
            V((i+1):(i+1+(Tau/deltaT))) = Vr;
        end 
    end
    spikes
    if spikes < 1
        fr(r) = 0;
    else     
        fr(r) = 1/(mean(diff(Tspikes)));
    end
end

fig5 = figure(5);
plot(runs, fr*1000)
xlabel('Input Current (nA)')
ylabel('Firing Rate (spks/s)')
title('f-I Curve')
