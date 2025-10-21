%% Part A: Deterministic model of the cro-cI genetic network
clear
close all
%% 1.(a) (all initial concentrations 0)
clear

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
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location', 'east')
hold off

%% 1.(b) (initial crorna 20)
clear

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

deltaT = 0.01; % s
maxT = 20; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

fig2 = figure(2);
cIprot = zeros(size(t));
cIrna = zeros(size(t));
croprot = zeros(size(t));
crorna = zeros(size(t));
cIprot(1) = 0; % initial concentration
cIrna(1) = 0; % initial concentration
croprot(1) = 0; % initial concentration
crorna(1) = 20; % initial concentration

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
title('Initial [cro_{RNA}] = 20 mlcls/cell')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location', 'east')
hold off

%% 1.(c) (initial cI rna 50)
clear

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

deltaT = 0.01; % s
maxT = 20; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

fig3 = figure(3);
cIprot = zeros(size(t));
cIrna = zeros(size(t));
croprot = zeros(size(t));
crorna = zeros(size(t));
cIprot(1) = 0; % initial concentration
cIrna(1) = 50; % initial concentration
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
title('Initial [cI_{RNA}] = 50 mlcls/cell')
legend('[cI_{protein}]', '[cI_{RNA}]', '[cro_{protein}]', '[cro_{RNA}]', 'Location', 'east')
hold off

%% 2. (a) (initial cIrna & crorna 0 to 20)
clear

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

deltaT = 0.01; % s
maxT = 50; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

fig4 = figure(4);
for cIrna_0 = 0:20 % initial concentration (mlcls/cell)
    for crorna_0 = 0:20 % initial concentration (mlcls/cell)
        cIprot = zeros(size(t));
        cIrna = zeros(size(t));
        croprot = zeros(size(t));
        crorna = zeros(size(t));
        cIprot(1) = 0; % initial concentration
        cIrna(1) = cIrna_0; % initial concentration
        croprot(1) = 0; % initial concentration
        crorna(1) = crorna_0; % initial concentration

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
        plot(cIprot, croprot)
        hold on
    end
end
xlabel('[cI_{protein}] (mlcls/cell)')
ylabel('[cro_{protein}] (mlcls/cell)')
title('Varying Initial [cI_{RNA}] and [cro_{RNA}] from 0 to 20 mlcl/cell')
xlim([0 500])
ylim([0 1000])
hold off

%% 2. (b) (initial cIrna & crorna 0 to 2000)
clear

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

deltaT = 0.01; % s
maxT = 100; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

fig5 = figure(5);
for cIrna_0 = 0:500:2000 % initial concentration (mlcls/cell)
    for crorna_0 = 0:500:2000 % initial concentration (mlcls/cell)
        cIprot = zeros(size(t));
        cIrna = zeros(size(t));
        croprot = zeros(size(t));
        crorna = zeros(size(t));
        cIprot(1) = 0; % initial concentration
        cIrna(1) = cIrna_0; % initial concentration
        croprot(1) = 0; % initial concentration
        crorna(1) = crorna_0; % initial concentration

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
        plot(cIprot, croprot)
        hold on
    end
end
xlabel('[cI_{protein}] (mlcls/cell)')
ylabel('[cro_{protein}] (mlcls/cell)')
title('Varying Initial [cI_{RNA}] and [cro_{RNA}] from 0 to 2000 mlcl/cell')
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

deltaT = 0.01; % s
maxT = 300; % s
numiterations = maxT/deltaT;
t = 0:deltaT:maxT;

fig1 = figure(1);
cIprot = zeros(size(t));
cIrna = zeros(size(t));
croprot = zeros(size(t));
crorna = zeros(size(t));
cIprot(1) = 1735.8; % initial concentration
cIrna(1) = 41.6597; % initial concentration
croprot(1) = 0.1296; % initial concentration
crorna(1) = 0.0021; % initial concentration

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

