clc;
close all;
clear all;
%%
sigma1 = 1; % Standard deviation for x1
sigma2 = 1; % Standard deviation for x2
A1 = 5; % Mean offset for x1
A2 = 5; % Mean offset for x2
% Define the grid
x1 = linspace(-5, 10, 100); % Range for x1
x2 = linspace(-5, 10, 100); % Range for x2
[X1, X2] = meshgrid(x1, x2);

% Define the first joint PDF (centered at (0, 0))
joint_pdf1 = @(x1, x2) (1 / (sqrt(2 * pi) * sigma1)) .* exp(-((x1).^2) / (2 * sigma1^2)) .* ...
                       (1 / (sqrt(2 * pi) * sigma2)) .* exp(-((x2).^2) / (2 * sigma2^2));

% Evaluate the first joint PDF
Z1 = joint_pdf1(X1, X2);

joint_pdf2 = @(x1, x2) (1 / (sqrt(2 * pi) * sigma1)) .* exp(-((x1 - A1).^2) / (2 * sigma1^2)) .* ...
                       (1 / (sqrt(2 * pi) * sigma2)) .* exp(-((x2 - A2).^2) / (2 * sigma2^2));

% Evaluate the second joint PDF
Z2 = joint_pdf2(X1, X2);

figure;
contour(X1, X2, Z1, 15, 'LineWidth', 1.5); % Contour with 15 levels
hold on;

% Plot contour for the second joint PDF
contour(X1, X2, Z2, 15, 'LineWidth', 1.5);

% Add labels and legend
% title('Contours of Joint PDFs');
xline(0,"color",'black')
yline(0,"color",'black')
colorbar()
xlabel('x_1');
ylabel('x_2');
grid on;
hold off;
set(gca, "FontSize", 10, "FontName", "Times New Roman")
print_plot("1.png",4,3,800)