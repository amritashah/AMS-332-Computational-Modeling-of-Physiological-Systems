%% Part A: Poisson spike trains and basic measures of neural activity
clear
close all

%% 1-6
% 1. Spike trains
T = 5; % s
N = 50; % spikes
lambda = 10; % spikes/s
ISI = zeros((2*lambda*T),N);
S = zeros((2*lambda*T), N);
for i = 1:N
    u = rand((2*lambda*T),1);
    ISI(:, i) = -log(u) / lambda;
    S(:, i) = cumsum(ISI(:, i));
end
S(S>T)=NaN;
    
% 2. Raster plot
fig1 = figure(1);
plot(S, 1:N, '.k')
xlabel('Time (s)')
ylabel('Trial Number')
title('Raster Plot')

% 3. Firing rate
f = sum(~isnan(S), 'all') / (N*T);

% 4. PSTH
fTi = [];
Ti = 0.1:0.2:T;
for i = 0:0.2:(T-0.2)
    fTi(end+1) = sum((S>=i & S<(i+0.2)), 'all');
end
fTi = fTi ./ (N*0.2);
fig2 = figure(2);
bar(Ti, fTi, 'hist')
xlim([0 T])
xlabel('Time (s)')
ylabel('Firing Rate (spikes/s)')
title('Peri-Stimulus Time Histogram (PSTH) Plot')

% 5. CV
CVk = zeros(1, N);
for k = 1:N
    CVk(k) = std(ISI(:, k)) / mean(ISI(:, k));
end
k = 1:N;
fig3 = figure(3);
plot(k, CVk)
xlabel('Trial Number')
ylabel('Coefficient of Variability (CV)')
title('CV_{k}')
CV = mean(CVk);

% 6. Fano factor
spike_count = zeros(1, N);
for k = 1:N
    spike_count(k) = sum(~isnan(S(:, k)))';
end
FF = var(spike_count) / mean(spike_count);

f
CV
FF

%% 7. Repeat 5 & 6
durations = [5, 10, 100, 200, 400, 800, 1600, 3200, 6400];
CV_vector = zeros(size(durations));
FF_vector = zeros(size(durations));

for d = 1:length(durations)
    % 1. Spike trains
    T = durations(d); % s
    N = 50; % spikes
    lambda = 10; % spikes/s
    ISI = zeros((2*lambda*T),N);
    S = zeros((2*lambda*T), N);
    for i = 1:N
        u = rand((2*lambda*T),1);
        ISI(:, i) = -log(u) / lambda;
        S(:, i) = cumsum(ISI(:, i));
    end
    S(S>T)=NaN;
           
    % 3. Firing rate
    f = sum(~isnan(S), 'all') / (N*T);
    
    % 5. CV
    CVk = zeros(1, N);
    for k = 1:N
        CVk(k) = std(ISI(:, k)) / mean(ISI(:, k));
    end
    k = 1:N;
    CV = mean(CVk);
    CV_vector(d) = CV;
    
    % 6. Fano factor
    spike_count = zeros(1, N);
    for k = 1:N
        spike_count(k) = sum(~isnan(S(:, k)))';
    end
    FF = var(spike_count) / mean(spike_count);
    FF_vector(d) = FF;

end
fig4 = figure(4);
plot(durations, CV_vector)
hold on
plot(durations, FF_vector)
xlabel('Interval Duration (s)')
legend('CV', 'FF')
title('CV and FF as a function of Interval Duration')
hold off