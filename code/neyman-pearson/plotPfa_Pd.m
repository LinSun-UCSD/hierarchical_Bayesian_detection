clear all;
close all;
clc;
%%
% Parameters of the bivariate Gaussian distribution
mu = [0, 0]; % Mean vector [mu_x1, mu_x2]
Sigma = [1, 0; 0, 1]; % Covariance matrix

% Number of samples to generate
n = 10000000;

% Draw samples
samples = mvnrnd(mu, Sigma, n);

% Extract x1 and x2
x1 = samples(:, 1);
x2 = samples(:, 2);

% transform to samples in z
z = zeros(n,1);
theta1 = 5;
theta2 = 5;
for i = 1:n
    z(i) = theta1 * x1(i) + theta2 * x2(i) - 0.5 * (theta1^2 + theta2^2);
end
N = n;
pfa = 0.000001;
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;
index = round((1-pfa) * N);
z_value = z_sorted(index);
figure()
plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "r");
hold on;
plot(z_sorted(index), temp_pfas(index), 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r','HandleVisibility', 'off');
% N = length(samples);
%% plot pd


% Parameters of the bivariate Gaussian distribution
mu = [5, 5]; % Mean vector [mu_x1, mu_x2]
Sigma = [1, 0; 0, 1]; % Covariance matrix

% Draw samples
samples = mvnrnd(mu, Sigma, n);

% Extract x1 and x2
x1 = samples(:, 1);
x2 = samples(:, 2);

% transform to samples in z
z = zeros(n,1);
theta1 = mu(1);
theta2 = mu(2);
for i = 1:n
    z(i) = theta1 * x1(i) + theta2 * x2(i) - 0.5 * (theta1^2 + theta2^2);
end
N = n;
% pfa = 0.1;
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pds = 1-cdf_values;

% z_value = z_sorted(index);
% figure()
plot(z_sorted, temp_pds, "LineWidth", 2, "Color", "b");
hold on;
pd = interp1(z_sorted, temp_pds, z_value);
plot(z_value, pd, 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b','HandleVisibility', 'off');
% N = length(samples);
xlabel("z");
ylabel("Probability");
legend({"p_{FA}", "p_{D}"});
set(gca, "FontSize", 10, "FontName", "Times New Roman");
grid on;
print_plot("1.png", 4, 3, 800)
