%% Part A: A diffusing autoregulatory gene
%% 1. initial RNA & protein = 0.5mM everywhere
clear
close all

% constants
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
Dp = 1E-4; % um/s
Dr = 1E-4; % um/s

% system parameters
totalT = 30; % s
deltaT = 0.01; % s
deltaX = 0.02; % um
Xmin = 0; % um
Xmax = 3.0; % um

% numbers of points
spacesteps = ceil( (Xmax-Xmin) / deltaX) + 1;
timesteps = ceil(totalT/deltaT);

% set arrays & initial conditions
Xprot = 0.5.*(ones(spacesteps, timesteps));
Xrna = 0.5.*(ones(spacesteps, timesteps));

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) ); % rna synthesis
        v2 = Xr*Xrna(x,t); % rna degradation
        v3 = w*Xrna(x,t); % protein synthesis
        v4 = Xp*Xprot(x,t); % protein degradation
        
        % RNA Laplacian
        del2rna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2rna = del2rna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2rna = del2rna + 0 + Xrna(x-1, t);
        else
            del2rna = del2rna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2rna = del2rna/(deltaX^2); % Laplacian

        % RNA derivative
        drna_dt = v1 - v2 + Dr*del2rna;

        % RNA update
        Xrna(x,t+1) = Xrna(x,t) + (drna_dt*deltaT);

        % protein Laplacian
        del2prot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2prot = del2prot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2prot = del2prot + 0 + Xprot(x-1, t);
        else
            del2prot = del2prot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2prot = del2prot/(deltaX^2); % Laplacian

        % protein derivative
        dprot_dt = v3 - v4 + Dp*del2prot;

        % protein update
        Xprot(x,t+1) = Xprot(x,t) + (dprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)))
plot(time, (Xprot(1,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location', 'southeast')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna(2,:)))
plot(time, (Xprot(2,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Second Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location', 'southeast')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(3,:)))
plot(time, (Xprot(3,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Third Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location', 'southeast')
hold off

fig2 = figure(2);
subplot(1,2,1)
hold on
plot(space, (Xrna(:,1)))
plot(space, (Xprot(:,1)), '--')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', 'Location', 'southeast')
hold off

subplot(1,2,2)
hold on
plot(space, (Xrna(:,timesteps)))
plot(space, (Xprot(:,timesteps)))
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', 'Location', 'southeast')
hold off

%% 2. initial RNA & protein = 0.5mM at 1st position, 0 elsewhere
clear
close all

% constants
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
Dp = 1E-4; % um/s
Dr = 1E-4; % um/s

% system parameters
totalT = 30; % s
deltaT = 0.01; % s
deltaX = 0.02; % um
Xmin = 0; % um
Xmax = 3.0; % um

% numbers of points
spacesteps = ceil( (Xmax-Xmin) / deltaX) + 1;
timesteps = ceil(totalT/deltaT);

% set arrays & initial conditions
Xprot = (zeros(spacesteps, timesteps));
Xrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 0.5;
Xrna(1,1) = 0.5;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) ); % rna synthesis
        v2 = Xr*Xrna(x,t); % rna degradation
        v3 = w*Xrna(x,t); % protein synthesis
        v4 = Xp*Xprot(x,t); % protein degradation
        
        % RNA Laplacian
        del2rna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2rna = del2rna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2rna = del2rna + 0 + Xrna(x-1, t);
        else
            del2rna = del2rna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2rna = del2rna/(deltaX^2); % Laplacian

        % RNA derivative
        drna_dt = v1 - v2 + Dr*del2rna;

        % RNA update
        Xrna(x,t+1) = Xrna(x,t) + (drna_dt*deltaT);

        % protein Laplacian
        del2prot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2prot = del2prot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2prot = del2prot + 0 + Xprot(x-1, t);
        else
            del2prot = del2prot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2prot = del2prot/(deltaX^2); % Laplacian

        % protein derivative
        dprot_dt = v3 - v4 + Dp*del2prot;

        % protein update
        Xprot(x,t+1) = Xprot(x,t) + (dprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)))
plot(time, (Xprot(1,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna(2,:)))
plot(time, (Xprot(2,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Second Position in Space')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(3,:)))
plot(time, (Xprot(3,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Third Position in Space')
legend('[X_{rna}]', '[X_{protein}]')
hold off

fig2 = figure(2);
subplot(1,2,1)
hold on
plot(space, (Xrna(:,1)))
plot(space, (Xprot(:,1)), '--')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,2,2)
hold on
plot(space, (Xrna(:,timesteps)))
plot(space, (Xprot(:,timesteps)))
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

%% 3. initial RNA & protein = 1.0mM at 1st position, 0 elsewhere
clear
close all

% constants
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
Dp = 1E-4; % um/s
Dr = 1E-4; % um/s

% system parameters
totalT = 30; % s
deltaT = 0.01; % s
deltaX = 0.02; % um
Xmin = 0; % um
Xmax = 3.0; % um

% numbers of points
spacesteps = ceil( (Xmax-Xmin) / deltaX) + 1;
timesteps = ceil(totalT/deltaT);

% set arrays & initial conditions
Xprot = (zeros(spacesteps, timesteps));
Xrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) ); % rna synthesis
        v2 = Xr*Xrna(x,t); % rna degradation
        v3 = w*Xrna(x,t); % protein synthesis
        v4 = Xp*Xprot(x,t); % protein degradation
        
        % RNA Laplacian
        del2rna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2rna = del2rna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2rna = del2rna + 0 + Xrna(x-1, t);
        else
            del2rna = del2rna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2rna = del2rna/(deltaX^2); % Laplacian

        % RNA derivative
        drna_dt = v1 - v2 + Dr*del2rna;

        % RNA update
        Xrna(x,t+1) = Xrna(x,t) + (drna_dt*deltaT);

        % protein Laplacian
        del2prot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2prot = del2prot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2prot = del2prot + 0 + Xprot(x-1, t);
        else
            del2prot = del2prot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2prot = del2prot/(deltaX^2); % Laplacian

        % protein derivative
        dprot_dt = v3 - v4 + Dp*del2prot;

        % protein update
        Xprot(x,t+1) = Xprot(x,t) + (dprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)))
plot(time, (Xprot(1,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna(2,:)))
plot(time, (Xprot(2,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Second Position in Space')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(3,:)))
plot(time, (Xprot(3,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Third Position in Space')
legend('[X_{rna}]', '[X_{protein}]')
hold off

fig2 = figure(2);
subplot(1,2,1)
hold on
plot(space, (Xrna(:,1)))
plot(space, (Xprot(:,1)), '--')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,2,2)
hold on
plot(space, (Xrna(:,timesteps)))
plot(space, (Xprot(:,timesteps)))
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

%% 4. (a) total time = 50 s
clear
close all

% constants
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
Dp = 1E-4; % um/s
Dr = 1E-4; % um/s

% system parameters
totalT = 50; % s
deltaT = 0.01; % s
deltaX = 0.02; % um
Xmin = 0; % um
Xmax = 3.0; % um

% numbers of points
spacesteps = ceil( (Xmax-Xmin) / deltaX) + 1;
timesteps = ceil(totalT/deltaT);

% set arrays & initial conditions
Xprot = (zeros(spacesteps, timesteps));
Xrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) ); % rna synthesis
        v2 = Xr*Xrna(x,t); % rna degradation
        v3 = w*Xrna(x,t); % protein synthesis
        v4 = Xp*Xprot(x,t); % protein degradation
        
        % RNA Laplacian
        del2rna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2rna = del2rna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2rna = del2rna + 0 + Xrna(x-1, t);
        else
            del2rna = del2rna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2rna = del2rna/(deltaX^2); % Laplacian

        % RNA derivative
        drna_dt = v1 - v2 + Dr*del2rna;

        % RNA update
        Xrna(x,t+1) = Xrna(x,t) + (drna_dt*deltaT);

        % protein Laplacian
        del2prot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2prot = del2prot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2prot = del2prot + 0 + Xprot(x-1, t);
        else
            del2prot = del2prot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2prot = del2prot/(deltaX^2); % Laplacian

        % protein derivative
        dprot_dt = v3 - v4 + Dp*del2prot;

        % protein update
        Xprot(x,t+1) = Xprot(x,t) + (dprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)))
plot(time, (Xprot(1,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna(2,:)))
plot(time, (Xprot(2,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Second Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(3,:)))
plot(time, (Xprot(3,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Third Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

fig2 = figure(2);
subplot(1,2,1)
hold on
plot(space, (Xrna(:,1)))
plot(space, (Xprot(:,1)), '--')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,2,2)
hold on
plot(space, (Xrna(:,timesteps)))
plot(space, (Xprot(:,timesteps)))
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

%% 4. (b) total time = 100 s
clear
close all

% constants
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
Dp = 1E-4; % um/s
Dr = 1E-4; % um/s

% system parameters
totalT = 100; % s
deltaT = 0.01; % s
deltaX = 0.02; % um
Xmin = 0; % um
Xmax = 3.0; % um

% numbers of points
spacesteps = ceil( (Xmax-Xmin) / deltaX) + 1;
timesteps = ceil(totalT/deltaT);

% set arrays & initial conditions
Xprot = (zeros(spacesteps, timesteps));
Xrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) ); % rna synthesis
        v2 = Xr*Xrna(x,t); % rna degradation
        v3 = w*Xrna(x,t); % protein synthesis
        v4 = Xp*Xprot(x,t); % protein degradation
        
        % RNA Laplacian
        del2rna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2rna = del2rna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2rna = del2rna + 0 + Xrna(x-1, t);
        else
            del2rna = del2rna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2rna = del2rna/(deltaX^2); % Laplacian

        % RNA derivative
        drna_dt = v1 - v2 + Dr*del2rna;

        % RNA update
        Xrna(x,t+1) = Xrna(x,t) + (drna_dt*deltaT);

        % protein Laplacian
        del2prot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2prot = del2prot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2prot = del2prot + 0 + Xprot(x-1, t);
        else
            del2prot = del2prot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2prot = del2prot/(deltaX^2); % Laplacian

        % protein derivative
        dprot_dt = v3 - v4 + Dp*del2prot;

        % protein update
        Xprot(x,t+1) = Xprot(x,t) + (dprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)))
plot(time, (Xprot(1,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna(2,:)))
plot(time, (Xprot(2,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Second Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(3,:)))
plot(time, (Xprot(3,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Third Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

fig2 = figure(2);
subplot(1,2,1)
hold on
plot(space, (Xrna(:,1)))
plot(space, (Xprot(:,1)), '--')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,2,2)
hold on
plot(space, (Xrna(:,timesteps)))
plot(space, (Xprot(:,timesteps)))
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

%% 4. (c) total time = 200 s
clear
close all

% constants
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
Dp = 1E-4; % um/s
Dr = 1E-4; % um/s

% system parameters
totalT = 200; % s
deltaT = 0.01; % s
deltaX = 0.02; % um
Xmin = 0; % um
Xmax = 3.0; % um

% numbers of points
spacesteps = ceil( (Xmax-Xmin) / deltaX) + 1;
timesteps = ceil(totalT/deltaT);

% set arrays & initial conditions
Xprot = (zeros(spacesteps, timesteps));
Xrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) ); % rna synthesis
        v2 = Xr*Xrna(x,t); % rna degradation
        v3 = w*Xrna(x,t); % protein synthesis
        v4 = Xp*Xprot(x,t); % protein degradation
        
        % RNA Laplacian
        del2rna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2rna = del2rna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2rna = del2rna + 0 + Xrna(x-1, t);
        else
            del2rna = del2rna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2rna = del2rna/(deltaX^2); % Laplacian

        % RNA derivative
        drna_dt = v1 - v2 + Dr*del2rna;

        % RNA update
        Xrna(x,t+1) = Xrna(x,t) + (drna_dt*deltaT);

        % protein Laplacian
        del2prot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2prot = del2prot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2prot = del2prot + 0 + Xprot(x-1, t);
        else
            del2prot = del2prot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2prot = del2prot/(deltaX^2); % Laplacian

        % protein derivative
        dprot_dt = v3 - v4 + Dp*del2prot;

        % protein update
        Xprot(x,t+1) = Xprot(x,t) + (dprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)))
plot(time, (Xprot(1,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna(2,:)))
plot(time, (Xprot(2,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Second Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(3,:)))
plot(time, (Xprot(3,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Third Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

fig2 = figure(2);
subplot(1,2,1)
hold on
plot(space, (Xrna(:,1)))
plot(space, (Xprot(:,1)), '--')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

subplot(1,2,2)
hold on
plot(space, (Xrna(:,timesteps)))
plot(space, (Xprot(:,timesteps)))
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]')
hold off

%% 4. (d) total time = 400 s
clear
close all

% constants
u = 1; % 1/s
w = 1; % 1/s
Xp = 1; % 1/s
Xr = 1; % 1/s
K = 0.33; % mM
Dp = 1E-4; % um/s
Dr = 1E-4; % um/s

% system parameters
totalT = 400; % s
deltaT = 0.01; % s
deltaX = 0.02; % um
Xmin = 0; % um
Xmax = 3.0; % um

% numbers of points
spacesteps = ceil( (Xmax-Xmin) / deltaX) + 1;
timesteps = ceil(totalT/deltaT);

% set arrays & initial conditions
Xprot = (zeros(spacesteps, timesteps));
Xrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) ); % rna synthesis
        v2 = Xr*Xrna(x,t); % rna degradation
        v3 = w*Xrna(x,t); % protein synthesis
        v4 = Xp*Xprot(x,t); % protein degradation
        
        % RNA Laplacian
        del2rna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2rna = del2rna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2rna = del2rna + 0 + Xrna(x-1, t);
        else
            del2rna = del2rna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2rna = del2rna/(deltaX^2); % Laplacian

        % RNA derivative
        drna_dt = v1 - v2 + Dr*del2rna;

        % RNA update
        Xrna(x,t+1) = Xrna(x,t) + (drna_dt*deltaT);

        % protein Laplacian
        del2prot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2prot = del2prot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2prot = del2prot + 0 + Xprot(x-1, t);
        else
            del2prot = del2prot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2prot = del2prot/(deltaX^2); % Laplacian

        % protein derivative
        dprot_dt = v3 - v4 + Dp*del2prot;

        % protein update
        Xprot(x,t+1) = Xprot(x,t) + (dprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)))
plot(time, (Xprot(1,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna(2,:)))
plot(time, (Xprot(2,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Second Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(3,:)))
plot(time, (Xprot(3,:)))
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Third Position in Space')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southeast')
hold off

fig2 = figure(2);
subplot(1,2,1)
hold on
plot(space, (Xrna(:,1)))
plot(space, (Xprot(:,1)), '--')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southwest')
hold off

subplot(1,2,2)
hold on
plot(space, (Xrna(:,timesteps)))
plot(space, (Xprot(:,timesteps)))
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', 'Location','southwest')
hold off
