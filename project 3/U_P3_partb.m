%% Part B: A diffusing pair of mutually-inhibiting genes
%% 1-2. first position: initial X = 1mM. last position: initial Y = 1mM
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
totalT = 5; % s
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Yprot(spacesteps,1) = 1.0;
Yrna(spacesteps,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

%% 3. (a) total time = 50 s
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Yprot(spacesteps,1) = 1.0;
Yrna(spacesteps,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

%% 3. (b) total time = 100 s
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Yprot(spacesteps,1) = 1.0;
Yrna(spacesteps,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

%% 3. (c) total time = 200 s
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Yprot(spacesteps,1) = 1.0;
Yrna(spacesteps,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off


%% 3. (d) total time = 400 s
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Yprot(spacesteps,1) = 1.0;
Yrna(spacesteps,1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Xprot(spacesteps,1) = 1.0;
Xrna(spacesteps,1) = 1.0;
Yprot((ceil(spacesteps/2)),1) = 1.0;
Yrna((ceil(spacesteps/2)),1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Xprot(spacesteps,1) = 1.0;
Xrna(spacesteps,1) = 1.0;
Yprot((ceil(spacesteps/2)),1) = 1.0;
Yrna((ceil(spacesteps/2)),1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Xprot(spacesteps,1) = 1.0;
Xrna(spacesteps,1) = 1.0;
Yprot((ceil(spacesteps/2)),1) = 1.0;
Yrna((ceil(spacesteps/2)),1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
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
Yprot = (zeros(spacesteps, timesteps));
Yrna = (zeros(spacesteps, timesteps));
Xprot(1,1) = 1.0;
Xrna(1,1) = 1.0;
Xprot(spacesteps,1) = 1.0;
Xrna(spacesteps,1) = 1.0;
Yprot((ceil(spacesteps/2)),1) = 1.0;
Yrna((ceil(spacesteps/2)),1) = 1.0;

% loop over time and x
for t = 1:timesteps
    for x = 1:spacesteps
        v1 = u*(1-( (Yprot(x,t)^2) / ((K^2)+(Yprot(x,t)^2)) )); % X rna synthesis
        v2 = Xr*Xrna(x,t); % X rna degradation
        v3 = w*Xrna(x,t); % X protein synthesis
        v4 = Xp*Xprot(x,t); % X protein degradation
        v5 = u*(1-( (Xprot(x,t)^2) / ((K^2)+(Xprot(x,t)^2)) )); % Y rna synthesis
        v6 = Xr*Yrna(x,t); % Y rna degradation %document shows Xr*Xrna - typo?
        v7 = w*Yrna(x,t); % Y protein synthesis
        v8 = Xp*Yprot(x,t); % Y protein degradation
                
        % X RNA Laplacian
        del2Xrna = -2*Xrna(x,t); %reset value
        if (x==1)
            del2Xrna = del2Xrna + Xrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xrna = del2Xrna + 0 + Xrna(x-1, t);
        else
            del2Xrna = del2Xrna + Xrna(x+1,t) + Xrna(x-1,t);
        end
        del2Xrna = del2Xrna/(deltaX^2); % Laplacian

        % X RNA derivative
        dXrna_dt = v1 - v2 + Dr*del2Xrna;

        % X RNA update
        Xrna(x,t+1) = Xrna(x,t) + (dXrna_dt*deltaT);

        % X protein Laplacian
        del2Xprot = -2*Xprot(x,t); %reset value
        if (x==1)
            del2Xprot = del2Xprot + Xprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Xprot = del2Xprot + 0 + Xprot(x-1, t);
        else
            del2Xprot = del2Xprot + Xprot(x+1,t) + Xprot(x-1,t);
        end
        del2Xprot = del2Xprot/(deltaX^2); % Laplacian

        % X protein derivative
        dXprot_dt = v3 - v4 + Dp*del2Xprot;

        % X protein update
        Xprot(x,t+1) = Xprot(x,t) + (dXprot_dt*deltaT);



        % Y RNA Laplacian
        del2Yrna = -2*Yrna(x,t); %reset value
        if (x==1)
            del2Yrna = del2Yrna + Yrna(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yrna = del2Yrna + 0 + Yrna(x-1, t);
        else
            del2Yrna = del2Yrna + Yrna(x+1,t) + Yrna(x-1,t);
        end
        del2Yrna = del2Yrna/(deltaX^2); % Laplacian

        % Y RNA derivative
        dYrna_dt = v5 - v6 + Dr*del2Yrna;

        % Y RNA update
        Yrna(x,t+1) = Yrna(x,t) + (dYrna_dt*deltaT);

        % Y protein Laplacian
        del2Yprot = -2*Yprot(x,t); %reset value
        if (x==1)
            del2Yprot = del2Yprot + Yprot(x+1,t) + 0;
        elseif (x==spacesteps)
            del2Yprot = del2Yprot + 0 + Yprot(x-1, t);
        else
            del2Yprot = del2Yprot + Yprot(x+1,t) + Yprot(x-1,t);
        end
        del2Yprot = del2Yprot/(deltaX^2); % Laplacian

        % Y protein derivative
        dYprot_dt = v7 - v8 + Dp*del2Yprot;

        % Y protein update
        Yprot(x,t+1) = Yprot(x,t) + (dYprot_dt*deltaT);
    end
end

time = 0:deltaT:totalT;
space = 0:deltaX:Xmax;

fig1 = figure(1);
subplot(1,3,1)
hold on
plot(time, (Xrna(1,:)), 'b')
plot(time, (Xprot(1,:)), 'c')
plot(time, (Yrna(1,:)), 'r')
plot(time, (Yprot(1,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('First Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(time, (Xrna((ceil(spacesteps/2)),:)), 'b')
plot(time, (Xprot((ceil(spacesteps/2)),:)), 'c')
plot(time, (Yrna((ceil(spacesteps/2)),:)), '--r')
plot(time, (Yprot((ceil(spacesteps/2)),:)), '--m')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Middle Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(time, (Xrna(spacesteps,:)), 'b')
plot(time, (Xprot(spacesteps,:)), 'c')
plot(time, (Yrna(spacesteps,:)), 'r')
plot(time, (Yprot(spacesteps,:)), 'm')
xlabel('Time (s)')
ylabel('Concentration (mM)')
title('Last Position in Space')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

fig2 = figure(2);
subplot(1,3,1)
hold on
plot(space, (Xrna(:,1)),'b')
plot(space, (Xprot(:,1)), '--c')
plot(space, (Yrna(:,1)), 'r')
plot(space, (Yprot(:,1)), '--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('First Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,2)
hold on
plot(space, (Xrna(:,((timesteps/2)+1))), 'b')
plot(space, (Xprot(:,((timesteps/2)+1))),'--c')
plot(space, (Yrna(:,((timesteps/2)+1))),'r')
plot(space, (Yprot(:,((timesteps/2)+1))),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Middle Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off

subplot(1,3,3)
hold on
plot(space, (Xrna(:,timesteps)),'b')
plot(space, (Xprot(:,timesteps)),'--c')
plot(space, (Yrna(:,timesteps)),'r')
plot(space, (Yprot(:,timesteps)),'--m')
xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
title('Last Point in Time')
legend('[X_{rna}]', '[X_{protein}]', '[Y_{rna}]', '[Y_{protein}]')
hold off
