
close all;
clear all;
theta = load("theta.mat");
x1 = linspace(-20, 20, 20); % Range for x
x2 = linspace(-20, 20, 20); % Range for y
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

%% plot the mapping between x1, x2, and z
figure;

% Use surf for a 3D surface plot
surf(X1, X2, z');
xlabel('x_1');
ylabel('x_2');
zlabel('z');
box on;
% title('Mapping between x_1, x_2 and z');
% shading interp;  % Smooth shading
set(gca, "FontName", "Times New Roman", "FontSize", 10);
print_plot("1.png", 4, 3, 800)