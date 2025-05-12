clc;
close all;
clear all;
%%
close all
Cw = [1 0; 0 1];
Ctheta = [   1.6861   -0.0197;
   -0.0197    0.0880];
mu = [10.8202; 0.8962];
% Parameters of the Gaussian distribution
x1 = linspace(-5, 20, 100); % Range for x
x2 = linspace(-5, 20, 100); % Range for y
[X1, X2] = meshgrid(x1, x2);

% Compute the joint Gaussian probability density
pos = [X1(:), X2(:)];
pdf = mvnpdf(pos, mu', Cw + Ctheta); % MATLAB's multivariate normal PDF function
pdf = reshape(pdf, size(X1)); % Reshape to grid format

% Plot the contour
contour(x1, x2, pdf, 10, 'LineWidth', 1.2); % 20 contour levels

hold on;
pdf = mvnpdf(pos, [0 0], Cw); % MATLAB's multivariate normal PDF function
pdf = reshape(pdf, size(X1)); % Reshape to grid format
xlim([-5 15]);
% yticks(-5:2:5)
ylim([-5 15]);
grid on;
contour(x1, x2, pdf, 10, 'LineWidth', 1.2); % 20 contour levels
xlabel('x_1');
ylabel('x_2');
xline(0, "k");
yline(0, "k");
% title('Contour Plot of Joint Gaussian Distribution');
colorbar; % Optional: Adds a color bar
set(gca, "FontName", "Times New Roman", "FontSize", 10);
print_plot("1.png",4,3,800)