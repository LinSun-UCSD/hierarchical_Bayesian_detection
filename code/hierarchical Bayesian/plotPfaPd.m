close all;
clear all;
theta = load("theta.mat").theta;
x1 = linspace(-15, 15, 100); % Range for x
x2 = linspace(-15, 15, 100); % Range for y
[X1, X2] = meshgrid(x1, x2);
Cw = [1 0; 0 1];
% Compute the joint Gaussian probability density
pos = [X1(:), X2(:)];
denominator = zeros(length(x1), length(x2),1);
nominator = zeros(length(x1), length(x2),1);
temp1 = 0;
temp2 = 0;
z = zeros(length(x1), length(x2),1);

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
%%
pfa = 0.000001;
N = 1000000;
x = mvnrnd([0; 0], Cw, N);
temp_z = findZ(X1, X2, z, x);
z_sorted = sort(temp_z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;
index = round((1-pfa) * N);
z_value = z_sorted(index);

plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "r")
hold on;
plot(z_sorted(index), temp_pfas(index), 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r','HandleVisibility', 'off');
%% get pd
N = 100;
l = 100;
x = zeros(l, N,2);
for i = 1:l
    x(i,:,:) = mvnrnd(theta(i,:), Cw, N);
end
x = reshape(x,[],2);

temp_z = zeros(N,1);
for i = 1:N*l
    % Interpolate to find z at the specified (x1_query, x2_query)
    x1_query = x(i,1);
    x2_query = x(i,2);
    temp_z(i) = interp2(X1, X2, z, x1_query, x2_query, 'linear'); % 'linear' is the default method
end

z_sorted = sort(temp_z);
cdf_values = (1:length(z_sorted)) / length(z_sorted);
pds = 1-cdf_values;
%%
clc;
valid_idx = ~isnan(z_sorted) ;
z_sorted = z_sorted(valid_idx);
pds = pds(valid_idx);

plot(z_sorted, pds, '-','Color','b', 'LineWidth',2);
xlabel('z ');
ylabel('Pd');

temp = interp1(z_sorted, pds', z_value);
hold on;
plot(z_value, temp, 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b','HandleVisibility', 'off');
xlabel("z");
ylabel("Probability");
grid on;
legend({"p_{FA}", "p_{D}"});
set(gca, "FontSize", 10, "FontName", "Times New Roman");
%%
grid on;
print_plot("1.png", 4, 3, 800);

%%
close all;
figure();
scatter(x(:,1), x(:,2), 10,'k','.');
box on;
xlabel("x_1");
ylabel("x_2");
grid on;

set(gca, "FontSize", 10, "FontName", "Times New Roman");
grid on;
print_plot("2.png", 4, 3, 800);
%%
function temp_z = findZ(X1, X2, z, x)
N = length(x);
temp_z = zeros(N,1);
for i = 1:N
    x1_query = x(i,1);
    x2_query = x(i,2);
    temp_z(i) = interp2(X1, X2, z, x1_query, x2_query, 'linear'); % 'linear' is the default method
end
end