%% 5. f-I curve
clear
close all

runs = 98.5:0.5:123;
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
    
    sum = 0;
    for j = 1:(spikes-1)
        ISI = Tspikes(j+1) - Tspikes(j);
        sum = sum+ISI;
    end
        
    if spikes < 37
        fr(r) = 0;
    else     
        fr(r) = 1/(sum/spikes);
    end
    
end

fig5 = figure(5);
plot(runs, fr*1000)
xlabel('Input Current (pA)')
ylabel('Firing Rate (spks/s)')
title('f-I Curve')
