clc;
clear all;
close all;
%% 

data = readNPY([path 'ai211.npy']);
inclinometer = readNPY([path 'ai217.npy']);
inclinometer = inclinometer(:,7);
x1 = data(:,4);
x2 = inclinometer;
t = 1:length(x2);
fs = 10000;
t = t / fs;

x2d = wdenoise(x2, 10,DenoisingMethod="BlockJS");

temp = zeros(length(x2d), 1);
lag = 1000;
for i = 1:length(x2d)-lag;
    temp(i+lag) = x2(i+lag) - x2d(i);
    
end
%%
x1 = x1;
x2 = temp;
pfas = [0.0000001];
close all;
factor = 5;
len = factor*fs;
figure()
y = [];
indices = [];
% lag = 0.2 * fs;
% x2 = x2(lag:end);
% x1 = x1(1:end - lag + 1);
norm_x1 = [];
norm_x2 = [];
for j = 1 : length(pfas)
    % subplot(1, 1, j)
    pfa = pfas(j);

    t = (1:295*fs) / fs;
    y = [];
    indices = [];
    mu = [0, 0]; % Mean vector [mu_x1, mu_x2]
    Sigma = [1, 0; 0, 1]; % Covariance matrix
    
    % Number of samples to generate
    n = 10000000;
    
    % Draw samples
    samples = mvnrnd(mu, Sigma, n);
    
    % Extract x1 and x2
    temp_x1 = samples(:, 1);
    temp_x2 = samples(:, 2);
    
    % transform to samples in z
    z = zeros(n,1);
    theta1 = 5;
    theta2 = 5;
    for i = 1:n
        z(i) = theta1 * temp_x1(i) + theta2 * temp_x2(i) - 0.5 * (theta1^2 + theta2^2);
    end
    N = n;
    
    z_sorted = sort(z);
    cdf_values = (1:N) / N;
    temp_pfas = 1-cdf_values;
    index = round((1-pfa) * N);
    threshold_T = z_sorted(index);
    for i = 1:295/factor
        % denoise the inclinometer data
        window1 = x1((i-1)*len+1: (i) * len);
        window2 = x2((i-1)*len+1: (i) * len);
        window1 = window1 - mean(window1);
        window2 = window2 - mean(window2);
        sigma1 = std(window1);
        sigma2 = std(window2);
        window1 = window1./std(window1);
        window2 = window2./std(window2);
        norm_x1 = [norm_x1; window1];
        norm_x2 = [norm_x2; window2];
        window = theta1 * window1 + theta2 * window2 - 0.5 * (theta1^2 + theta2^2);
        index = find(window > threshold_T);
        % if (length(index) < 10)
        %     index= [];
        % end
        y = [y; window];
        shift = (i-1) *len;
        indices = [indices; index+shift];

    end

end
%% come back to original signal
close all;
figure()
subplot(3,1,1);
plot(t, y, 'k')
xlim([0 t(end)]);
set(gca, 'FontName', "Times New Roman", "FontSize", 10);
hold on;
yline(threshold_T, "color", "m", "LineWidth", 2);
plot(t(indices), y(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r')
grid on;
% ylabel("z");
% xlabel("Time [sec]");

subplot(3,1,2);
plot(t, x1(1:295*fs), 'color', 'k');
hold on;
grid on;
set(gca, 'FontName', "Times New Roman", "FontSize", 10);
plot(t(indices), x1(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
xlim([0 t(end)]);

subplot(3,1,3)
plot(t, inclinometer(1:295*fs), 'color', 'b', "LineWidth", 1.5);
hold on;
set(gca, 'FontName', "Times New Roman", "FontSize", 10);
xlim([0 t(end)]);
grid on;
plot(t(indices), inclinometer(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
xlabel("Time [sec]")

print_plot("1.png", 5, 4, 800)