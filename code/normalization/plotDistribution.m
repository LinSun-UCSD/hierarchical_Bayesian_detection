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
len = fs;
for i = 80
    window1 = x1((i-1)*len+1: (i) * len);
    window2 = x2((i-1)*len+1: (i) * len);
    % [window1, window2] = normalize(window1, window2);
    hist3([window1, window2], 'Nbins', [20 20]); % 2D histogram
    xlabel('x_{1}');
    ylabel('x_{2}');
    zlabel('Count');
    view(3);
    box on;
    set(gca, "FontName", "Times New Roman", "FontSize", 10);
   
end
print_plot("1.png", 4, 3, 800)
%% plot the PDFs of x1, x2

for i = 80
    window1 = x1((i-1)*len+1: (i) * len);
    window2 = x2((i-1)*len+1: (i) * len);
    % [window1, window2] = normalize(window1, window2);
    figure()
    [f,xi] = ksdensity(window1);
    plot(xi, f, "LineWidth", 1.2, "LineStyle", "-")
    hold on;
    [f,xi] = ksdensity(window2);
    plot(xi, f, "LineWidth", 1.2, "LineStyle", "-")
    legend({"x_{1}", "x_{2}"});
    ylabel("PDF");
    xlabel("Data")
    set(gca, "FontName", "Times New Roman", "FontSize", 10);
    
    % xticks(-5:1:5);
    grid on;
end
print_plot("1.png", 4, 2, 800)
%%
function [window1, window2] = normalize(window1, window2)
window1 = window1 - mean(window1);
window2 = window2 - mean(window2);
window1 = window1/std(window1);
window2 = window2 /std(window2);
end