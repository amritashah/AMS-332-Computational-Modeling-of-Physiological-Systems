%% Part B: Stochastic model of the cro-cI genetic network
%% 1. One simulation
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

%% 2. 20 simulations
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

fig7 = figure(7);
fig8 = figure(8);
for s = 1:20
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

    set(0, 'CurrentFigure', fig7)
    hold on
    plot(cIprot, croprot)
    
    set(0, 'CurrentFigure', fig8)
    hold on
    plot(time_points, cIprot, 'r')
    plot(time_points, cIrna, 'm')
    plot(time_points, croprot, 'b')
    plot(time_points, crorna, 'c')
end

set(0, 'CurrentFigure', fig7)
xlabel('[cI_{protein}] (mlcls/cell)')
ylabel('[cro_{protein}] (mlcls/cell)')
title('20 Simulations of Stochastic Model (Phase Plane)')
xlim([-50 1500])
ylim([-50 1500])
hold off

set(0, 'CurrentFigure', fig8)
xlabel('Time (s)')
ylabel('Concentration (mlcls/cell)')
title('20 Simulations of Stochastic Model (Over Time)')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location','northwest')
hold off

%% 3. crorna0 = 20
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

steps = 50000;
time_points(1) = 0;

cIprot = zeros(1,steps);
cIrna = zeros(1,steps);
croprot = zeros(1,steps);
crorna = zeros(1,steps);
cIprot(1) = 0; % initial concentration
cIrna(1) = 0; % initial concentration
croprot(1) = 0; % initial concentration
crorna(1) = 20; % initial concentration

v1 = zeros(1,steps);
v2 = zeros(1,steps);
v3 = zeros(1,steps);
v4 = zeros(1,steps);
v5 = zeros(1,steps);
v6 = zeros(1,steps);
v7 = zeros(1,steps);
v8 = zeros(1,steps);

fig9 = figure(9);
fig10 = figure(10);
for s = 1:20
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

    set(0, 'CurrentFigure', fig9)
    hold on
    plot(cIprot, croprot)
    
    set(0, 'CurrentFigure', fig10)
    hold on
    plot(time_points, cIprot, 'r')
    plot(time_points, cIrna, 'm')
    plot(time_points, croprot, 'b')
    plot(time_points, crorna, 'c')
end

set(0, 'CurrentFigure', fig9)
xlabel('[cI_{protein}] (mlcls/cell)')
ylabel('[cro_{protein}] (mlcls/cell)')
title('20 Simulations of Stochastic Model (Phase Plane)')
xlim([-50 1500])
ylim([-50 1500])
hold off

set(0, 'CurrentFigure', fig10)
xlabel('Time (s)')
ylabel('Concentration (mlcls/cell)')
title('20 Simulations of Stochastic Model (Over Time)')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location','northwest')
hold off

%% 3. cIrna0 = 20
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

steps = 100000;
time_points(1) = 0;

cIprot = zeros(1,steps);
cIrna = zeros(1,steps);
croprot = zeros(1,steps);
crorna = zeros(1,steps);
cIprot(1) = 0; % initial concentration
cIrna(1) = 20; % initial concentration
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

fig11 = figure(11);
fig12 = figure(12);
for s = 1:20
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

    set(0, 'CurrentFigure', fig11)
    hold on
    plot(cIprot, croprot)
    
    set(0, 'CurrentFigure', fig12)
    hold on
    plot(time_points, cIprot, 'r')
    plot(time_points, cIrna, 'm')
    plot(time_points, croprot, 'b')
    plot(time_points, crorna, 'c')
end

set(0, 'CurrentFigure', fig11)
xlabel('[cI_{protein}] (mlcls/cell)')
ylabel('[cro_{protein}] (mlcls/cell)')
title('20 Simulations of Stochastic Model (Phase Plane)')
xlim([-50 1500])
ylim([-50 1500])
hold off

set(0, 'CurrentFigure', fig12)
xlabel('Time (s)')
ylabel('Concentration (mlcls/cell)')
title('20 Simulations of Stochastic Model (Over Time)')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location','northwest')
hold off

%% 4. 
clear
close all

XcI_rna = 1.2; % 1/s
XcI_prot = 6.02; % 1/s
Xcro_rna = 0.8; % 1/s
Xcro_prot = 0.8; % 1/s
wcI = 50; % 1/s
ucI = 50; % 1/s
wcro = 50; % 1/s
ucro = 50; % 1/s
KcI = 10; % mlcl/cell*s
Kcro = 10; % mlcl/cell*s

steps = 50000;
time_points(1) = 0;

cIprot = zeros(1,steps);
cIrna = zeros(1,steps);
croprot = zeros(1,steps);
crorna = zeros(1,steps);
cIprot(1) = 1735.8; % initial concentration
cIrna(1) = 41.6597; % initial concentration
croprot(1) = 0.1296; % initial concentration
crorna(1) = 0.0021; % initial concentration

v1 = zeros(1,steps);
v2 = zeros(1,steps);
v3 = zeros(1,steps);
v4 = zeros(1,steps);
v5 = zeros(1,steps);
v6 = zeros(1,steps);
v7 = zeros(1,steps);
v8 = zeros(1,steps);

fig14 = figure(14);
fig15 = figure(15);
for s = 1:20
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

    set(0, 'CurrentFigure', fig14)
    hold on
    plot(cIprot, croprot)
    
    set(0, 'CurrentFigure', fig15)
    hold on
    plot(time_points, cIprot, 'r')
    plot(time_points, cIrna, 'm')
    plot(time_points, croprot, 'b')
    plot(time_points, crorna, 'c')
end

set(0, 'CurrentFigure', fig14)
xlabel('[cI_{protein}] (mlcls/cell)')
ylabel('[cro_{protein}] (mlcls/cell)')
title('20 Simulations of Stochastic Model (Phase Plane)')
xlim([0 1500])
ylim([0 1500])
hold off

set(0, 'CurrentFigure', fig15)
xlabel('Time (s)')
ylabel('Concentration (mlcls/cell)')
title('20 Simulations of Stochastic Model (Over Time)')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location','northwest')
hold off
