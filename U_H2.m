%%
clear
w = 6;
u = 4;
Xp = 3;
Xr = 2;
K = 2;

Xprot = 0:0.1:3;
Xrna1 = (Xp/w).*Xprot;
Xrna2 = (u/Xr).*((Xprot.^2) ./ ((K^2) + (Xprot.^2)));

plot(Xprot, Xrna1)
hold on
plot(Xprot, Xrna2)

xlabel('[X_{prot}]')
ylabel('[X_{rna}]')
hold off

%%
clear
w = 4;
u = 2;
Xp = 2;
Xr = 1;
K = 2;

deltaT = 0.01; % s
maxT = 20; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

Xrna_0 = u / (2*Xr); % initial concentration (mM)
Xprot_0 = K; % initial concentration (mM)
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

%% (all initial concentrations 0)
clear

XcI_rna = 0.8; % 1/s
XcI_prot = 0.8; % 1/s
Xcro_rna = 1.2; % 1/s
Xcro_prot = 1.2; % 1/s
wcI = 50; % 1/s
ucI = 50; % 1/s
wcro = 40; % 1/s
ucro = 40; % 1/s
KcI = 10; % mlcl/cell*s
Kcro = 80; % mlcl/cell*s

deltaT = 0.01; % s
maxT = 20; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

fig1 = figure(1);
cIprot = zeros(size(t));
cIrna = zeros(size(t));
croprot = zeros(size(t));
crorna = zeros(size(t));
cIprot(1) = 0; % initial concentration
cIrna(1) = 0; % initial concentration
croprot(1) = 0; % initial concentration
crorna(1) = 0; % initial concentration

for i = 1:numiterations
    dcIprotdt = (wcI*cIrna(i)) - (XcI_prot*cIprot(i));
    dcIrnadt = ucI*(1 - ((croprot(i)^2) / ((Kcro^2) + (croprot(i)^2)))) - XcI_rna*cIrna(i);
    dcroprotdt = (wcro*crorna(i)) - (Xcro_prot*croprot(i));
    dcrornadt = ucro*(1 - ((cIprot(i)^2) / ((KcI^2) + (cIprot(i)^2)))) - Xcro_rna*crorna(i);

    cIprot(i+1) = cIprot(i) + dcIprotdt*deltaT;
    cIrna(i+1) = cIrna(i) + dcIrnadt*deltaT;
    croprot(i+1) = croprot(i) + dcroprotdt*deltaT;
    crorna(i+1) = crorna(i) + dcrornadt*deltaT;
end

plot(t, cIprot)
hold on
plot(t, cIrna)
plot(t, croprot)
plot(t, crorna)
xlabel('Time (s)')
ylabel('Concentration (mlcls/cell)')
title('All Initial Concentrations = 0 mlcls/cell')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location', 'northwest')
hold off

%% One simulation
clear
close all

XcI_rna = 0.8; % 1/s
XcI_prot = 0.8; % 1/s
Xcro_rna = 1.2; % 1/s
Xcro_prot = 1.2; % 1/s
wcI = 50; % 1/s
ucI = 50; % 1/s
wcro = 40; % 1/s
ucro = 40; % 1/s
KcI = 10; % mlcl/cell*s
Kcro = 80; % mlcl/cell*s

steps = 50000;
time_points(1) = 0;

cIprot = zeros(1,steps);
cIrna = zeros(1,steps);
croprot = zeros(1,steps);
crorna = zeros(1,steps);
cIprot(1) = 0; % initial concentration
cIrna(1) = 0; % initial concentration
croprot(1) = 0; % initial concentration
crorna(1) = 0; % initial concentration

v1 = zeros(1,steps);
v2 = zeros(1,steps);
v3 = zeros(1,steps);
v4 = zeros(1,steps);
v5 = zeros(1,steps);
v6 = zeros(1,steps);
v7 = zeros(1,steps);
v8 = zeros(1,steps);

for i = 1:(steps-1)
    v1 = wcI*cIrna(i);
    v2 = XcI_prot*cIprot(i);
    v3 = ucI*(1- ( ((croprot(i))^2) / ((Kcro^2)+((croprot(i))^2)) ));
    v4 = XcI_rna*cIrna(i);
    v5 = wcro*crorna(i);
    v6 = Xcro_prot*croprot(i);
    v7 = ucro*(1- ( ((cIprot(i))^2) / ((KcI^2)+((cIprot(i))^2)) ));
    v8 = Xcro_rna*crorna(i);
    vtotal = v1+v2+v3+v4+v5+v6+v7+v8; %calc total rxn rate

    y = rand();
    tau = -log(y)/(vtotal); %calc time to next rxn
    time_points(i+1) = time_points(i) + tau;

    y_mod = vtotal*rand();
    if y_mod <= v1
        cIprot(i+1) = cIprot(i) + 1; %cIprot syn
        cIrna(i+1) = cIrna(i);
        croprot(i+1) = croprot(i);
        crorna(i+1) = crorna(i);
    elseif y_mod <= (v1 + v2)
        cIprot(i+1) = cIprot(i) - 1; %cIprot deg
        cIrna(i+1) = cIrna(i);
        croprot(i+1) = croprot(i);
        crorna(i+1) = crorna(i);
    elseif y_mod <= (v1 + v2 + v3)
        cIprot(i+1) = cIprot(i);
        cIrna(i+1) = cIrna(i) + 1; %cIrna syn
        croprot(i+1) = croprot(i);
        crorna(i+1) = crorna(i);
    elseif y_mod <= (v1 + v2 + v3 + v4)
        cIprot(i+1) = cIprot(i);
        cIrna(i+1) = cIrna(i) - 1; %cIrna deg
        croprot(i+1) = croprot(i);
        crorna(i+1) = crorna(i);
    elseif y_mod <= (v1 + v2 + v3 + v4 + v5)
        cIprot(i+1) = cIprot(i);
        cIrna(i+1) = cIrna(i);
        croprot(i+1) = croprot(i) + 1; %croprot syn
        crorna(i+1) = crorna(i);
    elseif y_mod <= (v1 + v2 + v3 + v4 + v5 + v6)
        cIprot(i+1) = cIprot(i);
        cIrna(i+1) = cIrna(i);
        croprot(i+1) = croprot(i) - 1; %croprot deg
        crorna(i+1) = crorna(i);
    elseif y_mod <= (v1 + v2 + v3 + v4 + v5 + v6 + v7)
        cIprot(i+1) = cIprot(i);
        cIrna(i+1) = cIrna(i);
        croprot(i+1) = croprot(i);
        crorna(i+1) = crorna(i) + 1; %crorna syn
    elseif y_mod <= (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8)
        cIprot(i+1) = cIprot(i);
        cIrna(i+1) = cIrna(i);
        croprot(i+1) = croprot(i);
        crorna(i+1) = crorna(i) - 1; %crorna deg
    else
        Disp('error');

    end
end

fig6 = figure(6);
hold on
plot(time_points, cIprot)
plot(time_points, cIrna)
plot(time_points, croprot)
plot(time_points, crorna)
xlabel('Time (s)')
ylabel('Concentration (mlcls/cell)')
title('1 Simulation of Stochastic Model')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location','northwest')
hold off

%% project 2, part a, 4. plotting null clines - doesn't work
clear
close all

XcI_rna = 1.2; % 1/s
XcI_prot = 1.2; % 1/s
Xcro_rna = 0.8; % 1/s
Xcro_prot = 0.8; % 1/s
wcI = 50; % 1/s
ucI = 50; % 1/s
wcro = 50; % 1/s
ucro = 50; % 1/s
KcI = 10; % mlcl/cell*s
Kcro = 10; % mlcl/cell*s

cIprot_scale = 0:.1:10000;
croprot_scale = 0:.1:10000;

cIprot_nc = ( (ucI*wcI)/(XcI_rna*XcI_prot) ) .* (1- ( (croprot_scale.^2) ./ ((100)+(croprot_scale.^2)) ));
croprot_nc = ( (ucro*wcro)/(Xcro_rna*Xcro_prot) ) .* (1- ( (cIprot_scale.^2) ./ ((100)+(cIprot_scale.^2)) ));
hold on
plot(cIprot_nc, croprot_scale)
plot(cIprot_scale, croprot_nc)
xlabel('[cI_{protein}] (mlcls/cell)')
ylabel('[cro_{protein}] (mlcls/cell)')
xlim([0 200])
ylim([0 5000])
title('Null Clines')
hold off