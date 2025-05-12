clear all;
close all;
close all;

%%
% Parameters
theta1 = 5;
theta2 = 5;
x1 = -50:10:50; % Define x1
x2 = -50:10:50; % Define x2

% Create meshgrid for x1 and x2
[X1, X2] = meshgrid(x1, x2);

% Calculate z as a function of x1 and x2
Z = theta1 * X1 + theta2 * X2 - 0.5 * (theta1^2 + theta2^2);

% Plot the 3D surface
figure;
surf(X1, X2, Z); % Create a 3D surface plot
xlabel('x_1');
ylabel('x_2');
zlabel('z');
% title('3D Plot of z as a Function of x_1 and x_2');
colormap('jet'); % Apply a colormap
set(gca, "FontSize", 10, "FontName", "Times New Roman");
box on;
print_plot("1.png", 4, 3, 800)