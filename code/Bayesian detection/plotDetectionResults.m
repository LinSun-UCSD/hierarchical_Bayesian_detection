clc;
clear all;
close all;

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
%%
close all;
Ctheta = [   1.6861   -0.0197;
    -0.0197    0.0880];
mu = [10.8202; 0.8962];
Cw = [1 0; 0 1];

y = [];
indices = [];
pfas = [0.000001];
factor = 5;
len = fs*factor;
for j = 1 : length(pfas)
    subplot(2, 1, j)
    pfa = pfas(j);


    y = [];
    indices = [];
    % get thrshold
    N = 10000000;
    x = mvnrnd([0; 0], Cw, N);
    x = x';
    z = zeros(length(x), 1);
    for ii = 1 : length(x)
        z(ii) = getMapping(x(:, ii), mu, Cw, Ctheta);
    end
    z_sorted = sort(z);
    cdf_values = (1:N) / N;
    temp_pfas = 1-cdf_values;
    index = round((1-pfa) * N);

    z_value = z_sorted(index);
    for i = 1:296/factor
        % denoise the inclinometer data
        window1 = x1((i-1)*len+1: (i) * len);
        window2 = x2((i-1)*len+1: (i) * len);
        [window1, window2] = normalize(window1, window2);
        z = zeros(length(window1),1);
        for ii = 1:length(z)
            z(ii) = getMapping([window1(ii); window2(ii)], mu, Cw, Ctheta);
        end
        index = find(z> z_value);
        if (length(index) < 10)
            index= [];
        end
        y = [y; z];
        shift = (i-1) *len;
        indices = [indices; index+shift];
    end
end
%%
close all;
figure()
subplot(3,1,1)
t = (1:length(y)) / fs;
plot(t, y, 'k')
hold on;
plot(t(indices), y(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
% plot(t, T, 'r');
xlim([t(1) t(end)])
grid on;
set(gca, 'FontName', "Times New Roman", "FontSize", 10);
% come back to original dataset
path = 'D:\UCSD Post-doc\pythonCode\SooLocks\signal processing\00272281_00272575\'
data = readNPY([path 'ai211.npy']);
fs = 10000;

% range = 156*fs : 250 * fs;
range = 1:length(data);
x = data(:, 4);
x = x - mean(x);

% load inclinometer
inclinometer = readNPY([path 'ai217.npy']);
inclinometer = inclinometer(:,7);
x1 = x;
x2 = inclinometer;
t = (1:length(x1))/fs;

subplot(3,1,2)
plot(t, x1, 'k','LineWidth',1.2);
grid on;
hold on;
xlim([t(1) t(end)])
plot(t(indices), x1(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');

set(gca, 'FontName', "Times New Roman", "FontSize", 10);

subplot(3,1,3)
plot(t, x2,'b','LineWidth',1.2);
grid on;
hold on;
xlim([t(1) t(end)])
plot(t(indices), x2(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
xlabel("Time [sec]");

set(gca, 'FontName', "Times New Roman", "FontSize", 10);
print_plot("1.png",5,4,800)


%% utility functions
function z = getMapping(x, mu, Cw, Ctheta)

Sigma = Cw + Ctheta;  % Combined covariance
log_det_Sigma = log(det(Sigma));
log_det_Cw = log(det(Cw));
z = -0.5 * (log_det_Sigma) - 0.5*(x-mu)' *(Sigma \ (x-mu));
z = z + 0.5 * (log_det_Cw) + 0.5*x' * (Cw \ x);

end
function [window1, window2] = normalize(window1, window2)
window1 = window1 - mean(window1);
window2 = window2 - mean(window2);
sigma1 = std(window1);
sigma2 = std(window2);
max_spike1 = max(abs(window1));
max_spike2 = max(abs(window1));

% Step 2: Scale both channels to have the same maximum spike magnitude
desired_spike_magnitude = 1; % Use the larger spike as reference
scaled_x1 = window1 / max_spike1 * desired_spike_magnitude;
scaled_x2 = window2 / max_spike2 * desired_spike_magnitude;

% Step 3: Normalize the variance to 1
window1 = scaled_x1 / std(scaled_x1);
window2 = scaled_x2 / std(scaled_x2);


end