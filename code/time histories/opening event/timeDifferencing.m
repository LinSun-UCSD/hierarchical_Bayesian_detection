clc;
clear all;
close all;
%% 
addpath("D:\UCSD Post-doc\paper\paper 3\figures\code\Bayesian detection");

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
close all;
figure()
subplot(2,1,1)
plot(t, x1, "k");
grid minor;
xlim([35, 37])
title("Sector Gear Bearing Block")
ylim([-0.3 0.3])
set(gca, "FontName", "Times New Roman", "FontSize", 10)


subplot(2,1,2);
plot(t, x2, 'b');
hold on;
plot(t, x2d, "color", "r")
title("Sector Gear Inclinometer Time History")
set(gca, "FontName", "Times New Roman", "FontSize", 10);
grid minor;
xlim([35, 37]);
print_plot("1.png",4,4,800);
%%
close all;
figure()
plot(t, temp,'b');
title("First-order Time-differenced Angular Displacement Signal")
set(gca, "FontName", "Times New Roman", "FontSize", 10);
grid on;
xlim([0 296])
ylim([-0.6 0.6])
print_plot("1.png", 6, 2, 800)
