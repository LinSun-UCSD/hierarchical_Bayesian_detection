close all;
clear all;
theta = load("theta.mat");
x1 = linspace(-5, 15, 100); % Range for x
x2 = linspace(-5, 15, 100); % Range for y
[X1, X2] = meshgrid(x1, x2);
Cw = [1 0; 0 1];
% Compute the joint Gaussian probability density
pos = [X1(:), X2(:)];
denominator = zeros(length(x1), length(x2),1);
nominator = zeros(length(x1), length(x2),1);
        temp1 = 0;
        temp2 = 0;
z = zeros(length(x1), length(x2),1);
theta = mvnrnd([10; 1], [1 0; 0 1;], 100)
for i = 1:length(x1)
    for j = 1:length(x2)

        temp1 = sum((1 / (sqrt(2 * pi * Cw(1,1)))) * exp(-0.5 * ((x1(i) - theta(:,1)) / sqrt(Cw(1,1))).^2)...
            *(1 / (sqrt(2 * pi * Cw(2,2)))) .* exp(-0.5 * ((x2(j) - theta(:,2)) / sqrt(Cw(2,2))).^2));
        temp2 = sum((1 / (sqrt(2 * pi * Cw(1,1)))) * exp(-0.5 * ((x1(i)) / sqrt(Cw(1,1))).^2)...
            *(1 / (sqrt(2 * pi * Cw(2,2)))) .* exp(-0.5 * ((x2(j)) / sqrt(Cw(2,2))).^2));
        denominator(i,j) = temp1;
        nominator(i,j) = temp2;
        z(i,j) = log(temp1) - log(temp2);
    end
end

figure;
contour(X1, X2, denominator', 20); % Adjust '20' for the number of contour levels
xlabel('x_1');
ylabel('x_2');
nominator = reshape(nominator, size(X1));
% title('Contour Plot of Denominator');
grid on;
hold on;
box on;
contour(X1, X2, nominator', 10, 'LineWidth', 1.2);
xline(0, "k");
yline(0, "k");
colorbar; % Adds a colorbar for reference
set(gca,'FontSize',10, "FontName", "Times New Roman");
print_plot("1.png",4.3,3,800);