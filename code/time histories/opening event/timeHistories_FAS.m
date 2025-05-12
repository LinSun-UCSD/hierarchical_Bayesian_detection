clear all;
close all;
clc;


Fs = 10000;
t = (1 : length(data)) / Fs;
%%
titleName = {"Motor", "Speed Reducer Bearing Block",...
    "Jack Shaft Bearing Block", "Sector Gear Bearing Block"}
figure()
for i = 1:4
    subplot(5 ,1,i)
    plot(t, data(:,i),'Color', 'k');
    grid on;
    xlim([t(1) t(end)]);
    % xlabel("Time [sec]");
    % ylabel("Accel. [g]");
    title(titleName{i});
    set(gca, "FontName", "Times New Roman", "FontSize", 10);
end
print_plot("1.png", 6, 6.5, 800)
%%
subplot(5,1,5)
inclinometer = readNPY([path 'ai217.npy']);
inclinometer = inclinometer(:,7);
plot(t, inclinometer,'Color', 'b', "LineWidth", 1.5);
grid on;
xlim([t(1) t(end)]);
% xlabel("Time [sec]");
% ylabel("Accel. [g]");
title("Sector Gear Inclinometer Time History");
set(gca, "FontName", "Times New Roman", "FontSize", 10);
print_plot("1.png", 6, 8, 800)
%%

 % Normalize the magnitude spectrum
% s_mag = abs(s);                 % Take the magnitude
% s_norm = s_mag / max(s_mag(:)); % Normalize to range [0, 1]
% 
% % Plot normalized STFT magnitude spectrum
% figure;
% imagesc(t_stft, f, s_norm);     % Time vs frequency with normalized magnitude
% axis xy; % Correct axis orientation
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% ylim([0 5000])
% % title('Normalized STFT Magnitude Spectrum of Sector Gear Bearing Block');
% colorbar;
% caxis([0 0.1]);                 % Set color range from 0 to 0.5
% set(gca, "FontName", "Times New Roman", "FontSize", 10);
% print_plot("1.png", 6, 4, 800)
%%
close all
figure()
subplot(2,1,1)
plot(t, data(:, 4), 'k')
title("Sector Gear Bearing Block")
ylim([-0.3 0.3])
set(gca, "FontName", "Times New Roman", "FontSize", 10)
xlim([32 55]);
grid minor;

% ylabel("Accel. [g]")
subplot(2,1,2)
plot(t, inclinometer,'b', "LineWidth", 1);
xlim([32 55]);
title("Sector Gear Inclinometer Time History")
set(gca, "FontName", "Times New Roman", "FontSize", 10);
grid minor ;
% ylabel("Degree [^{o}]")
% ylabel("Accel. [g]")
print_plot("1.png",6,4,800);
%%
figure()
subplot(2,1,1)
plot(t, data(:, 4), 'k')
title("Sector Gear Bearing Block")
ylim([-0.3 0.3])
set(gca, "FontName", "Times New Roman", "FontSize", 10)
xlim([43  45]);
grid minor;

% ylabel("Accel. [g]")
subplot(2,1,2)
plot(t, inclinometer,'b', "LineWidth", 1);
xlim([43 45]);
title("Sector Gear Inclinometer Time History")
set(gca, "FontName", "Times New Roman", "FontSize", 10);
grid minor ;

xd = wdenoise(inclinometer, 10,DenoisingMethod="BlockJS");
hold on;
plot(t, xd,'r', "LineWidth", 1);
xlim([43 45]);
% ylabel("Degree [^{o}]")
% ylabel("Accel. [g]")
print_plot("1.png",4,4,800);
