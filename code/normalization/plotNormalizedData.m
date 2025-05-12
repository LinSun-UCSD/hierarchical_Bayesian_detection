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
figure();
plot(x1);
figure()
plot(x2)
%% plot the 
close all;
figure()
count = 0;
factor = 0.01;
len = fs*factor;
norm_x1 = [];
norm_x2 = [];

for i = 1:296/factor
    window1 = x1((i-1)*len+1: (i) * len);
    window2 = x2((i-1)*len+1: (i) * len);
    [w1, w2] = normalize(window1, window2);
    norm_x1 = [norm_x1; w1];
    norm_x2 = [norm_x2; w2];
   
end
%%
close all;
figure()
subplot(2,1,1)
plot((1:length(norm_x1))/fs, norm_x1, "Color", "k");
grid on;
xlim([0 length(norm_x1)/fs]);
xlabel("Time [sec]")
ylabel("x_{1}")
set(gca, "FontSize", 10, "FontName", "Times New Roman");
subplot(2,1,2)
plot((1:length(norm_x1))/fs, norm_x2, "color", "b");
grid on;
xlim([0 length(norm_x1)/fs]);
ylabel("x_{2}");
xlabel("Time [sec]")
set(gca, "FontSize", 10, "FontName", "Times New Roman");
print_plot("1.png", 5,4,800)


%%
function [window1, window2] = normalize(window1, window2)
window1 = window1 - mean(window1);
window2 = window2 - mean(window2);
window1 = window1/std(window1);
window2 = window2 /std(window2);
end