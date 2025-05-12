close all;
clear all;
% Define range for x1 and x2
x1 = linspace(-5, 30, 100); 
x2 = linspace(-5, 30, 100); 

% Initialize PDF matrix
pdf = zeros(length(x1), length(x2));

% Define covariance matrix (identity) and mean
Cw = [1 0; 0 1]; 
mu = [10, 5]; 

for i = 1:length(x1)
    for j = 1:length(x2)
        % Compute Gaussian PDF manually
        pdf(i,j) = (1 / (2 * pi)) * exp(-0.5 * ((x1(i) - mu(1))^2 + (x2(j) - mu(2))^2));
    end
end

% Create meshgrid for plotting
[X1, X2] = meshgrid(x1, x2);

% Plot contour
figure;
contour(X1, X2, pdf', 20, 'LineWidth', 1.5); % 20 contour levels
colorbar; % Add colorbar for reference

figure()
surf(X1, X2, pdf')
