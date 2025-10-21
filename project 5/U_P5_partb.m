%% Part B: Perceptron learning
%% 1. AND problem
clear

N = 2; % input units
eta = 1; % learning rate 
w = zeros(1, N+1); % initial synaptic weights
inputs = [1 -1 1 -1; 1 1 -1 -1; -1 -1 -1 -1]; % input patterns (with x3=-1)
presentations = 50;
performance = zeros(1, presentations);

for i = 1:presentations
    x = inputs(:, (randi(size(inputs, 2)))); % select an input pattern
    y = sign(dot(w, x)); % compute output

    if x(1)==1 && x(2)==1
        yt_and = 1;
    else
        yt_and = -1;
    end
    
    dw1 = eta*(yt_and-y)*x(1);
    dw2 = eta*(yt_and-y)*x(2);
    dw3 = eta*(yt_and-y)*x(3);

    w(1) = w(1) + dw1;
    w(2) = w(2) + dw2;
    w(3) = w(3) + dw3;

    if y == yt_and
        performance(i) = 1;
    else
        performance(i) = 0;
    end
end

% plot performance vs presentation number
fig1 = figure(1);
plot(1:presentations, performance, '.', 'MarkerSize', 10);
xlabel('Presentation number');
ylabel('Performance');
ylim([-0.2 1.2]);
title('AND problem Performance vs. Presentation Number')

%test
for i = 1:size(inputs, 2)
    x = inputs(:, i);
    y = sign(dot(w, x));
    fprintf('Input pattern: [%d %d %d], Output: %d\n', x(1), x(2), x(3), y);
end

%% OR problem
clear

N = 2; % input units
eta = 1; % learning rate 
w = zeros(1, N+1); % initial synaptic weights
inputs = [1 -1 1 -1; 1 1 -1 -1; -1 -1 -1 -1]; % input patterns (with x3=-1)
presentations = 50;
performance = zeros(1, presentations);

for i = 1:presentations
    x = inputs(:, (randi(size(inputs, 2)))); % select an input pattern
    y = sign(dot(w, x)); % compute output

    if x(1)==1 || x(2)==1
        yt_or = 1;
    else
        yt_or = -1;
    end
    
    dw1 = eta*(yt_or-y)*x(1);
    dw2 = eta*(yt_or-y)*x(2);
    dw3 = eta*(yt_or-y)*x(3);

    w(1) = w(1) + dw1;
    w(2) = w(2) + dw2;
    w(3) = w(3) + dw3;

    if y == yt_or
        performance(i) = 1;
    else
        performance(i) = 0;
    end
end

% plot performance vs presentation number
fig2 = figure(2);
plot(1:presentations, performance, '.', 'MarkerSize', 10);
xlabel('Presentation number');
ylabel('Performance');
ylim([-0.2 1.2]);
title('OR problem Performance vs. Presentation Number')

%test
for i = 1:size(inputs, 2)
    x = inputs(:, i);
    y = sign(dot(w, x));
    fprintf('Input pattern: [%d %d %d], Output: %d\n', x(1), x(2), x(3), y);
end

%% XOR problem
clear

N = 2; % input units
eta = 1; % learning rate 
w = zeros(1, N+1); % initial synaptic weights
inputs = [1 -1 1 -1; 1 1 -1 -1; -1 -1 -1 -1]; % input patterns (with x3=-1)
presentations = 2000;
performance = zeros(1, presentations);

for i = 1:presentations
    x = inputs(:, (randi(size(inputs, 2)))); % select an input pattern
    y = sign(dot(w, x)); % compute output
    
    if (x(1)==-1 && x(2)==1) || (x(1)==1 && x(2)==-1)
        yt_xor = 1;
    else
        yt_xor = -1;
    end
    
    dw1 = eta*(yt_xor-y)*x(1);
    dw2 = eta*(yt_xor-y)*x(2);
    dw3 = eta*(yt_xor-y)*x(3);

    w(1) = w(1) + dw1;
    w(2) = w(2) + dw2;
    w(3) = w(3) + dw3;

    % compute performance
    if y == yt_xor
        performance(i) = 1;
    else
        performance(i) = 0;
    end
end
% plot performance vs presentation number
fig3 = figure(3);
plot(1:presentations, performance, '.', 'MarkerSize', 10);
xlabel('Presentation number');
ylabel('Performance');
ylim([-0.2 1.2]);
title('XOR problem Performance vs. Presentation Number')

%test
for i = 1:size(inputs, 2)
    x = inputs(:, i);
    y = sign(dot(w, x));
    fprintf('Input pattern: [%d %d %d], XOR Output: %d\n', x(1), x(2), x(3), y);
end

%% 2. Random patterns
clear
M = 40;
N = 50;
inputs = randi(2,N,M) - 1;
inputs = inputs -1;
eta = 1; % learning rate 
w = zeros(N, 1); % initial synaptic weights
classes = [-ones(1,M/2) ones(1,M/2)];
classes = classes(randperm(length(classes)));
presentations = 1000;
performance = zeros(1, presentations);

for i = 1:presentations
    index = randi(size(inputs, 2));
    x = inputs(:,index);
    yt = classes(index);
    y = sign(dot(w, x)); % compute output
    
    if y == yt
        performance(i) = 1;
    else
        performance(i) = 0;
    end
    for j = 1:N
        dw = eta * (yt - y) * x(j);
        w(j) = w(j) + dw;
    end
end

fig4 = figure(4);
plot(1:presentations, performance, '.', 'MarkerSize', 10);
xlabel('Presentation number');
ylabel('Performance');
ylim([-0.2 1.2]);
title('Performance vs. Presentation Number')

%% 2d. 
clear
M = 90;
N = 50;
inputs = randi(2,N,M) - 1;
inputs = inputs -1;
eta = 1; % learning rate 
w = zeros(N, 1); % initial synaptic weights
classes = [-ones(1,M/2) ones(1,M/2)];
classes = classes(randperm(length(classes)));
presentations = 1000*M;
performance = zeros(1, presentations);

for i = 1:presentations
    index = randi(size(inputs, 2));
    x = inputs(:,index);
    yt = classes(index);
    y = sign(dot(w, x)); % compute output
    
    if y == yt
        performance(i) = 1;
    else
        performance(i) = 0;
    end
    for j = 1:N
        dw = eta * (yt - y) * x(j);
        w(j) = w(j) + dw;
    end

    if i>200 && all(performance(i-199 : i) == 1)
        disp(['n_{conv} =', num2str(i)]);
        break
    end
end

fig4 = figure(4);
plot(1:presentations, performance, '.', 'MarkerSize', 10);
xlabel('Presentation number');
ylabel('Performance');
ylim([-0.2 1.2]);
title('Performance vs. Presentation Number')
