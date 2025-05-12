clear all;
close all;
clc;

data = readNPY([path 'ai211.npy']);
fs = 10000;


%%
close all;
% plot pfa
Ctheta = [   1.6861   -0.0197;
    -0.0197    0.0880];
mu = [10.8202; 0.8962];
Cw = [1 0; 0 1];
figure()
pfa = 0.000001;
N = 10000000;
x = mvnrnd([0; 0], Cw, N);
x = x';
z = zeros(length(x), 1);
for ii = 1 : length(x)
    z(ii) = getMapping(x(:, ii), mu, Cw, Ctheta);
end
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;
index = round((1-pfa) * N);

z_value = z_sorted(index);

plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "r")
hold on;
plot(z_sorted(index), temp_pfas(index), 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r','HandleVisibility', 'off');
% plot pD
x = mvnrnd(mu, Cw+Ctheta, N);
x = x';
z = zeros(length(x), 1);
for ii = 1 : length(x)
    z(ii) = getMapping(x(:, ii), mu, Cw, Ctheta);
end
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pds = 1-cdf_values;


plot(z_sorted, temp_pds, "LineWidth", 2, "Color", "b")
hold on;
pd = interp1(z_sorted, temp_pds, z_value,  'linear', 'extrap');
plot(z_value, pd, 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b','HandleVisibility', 'off');
grid on;
xlabel("z");
ylabel("Probability");
legend({"p_{FA}", "p_{D}"});
set(gca, "FontSize", 10, "FontName", "Times New Roman");
%%
print_plot("1.png",4,3,800)

%% utility functions
function z = getMapping(x, mu, Cw, Ctheta)

Sigma = Cw + Ctheta;  % Combined covariance
log_det_Sigma = log(det(Sigma));
log_det_Cw = log(det(Cw));
z = -0.5 * (log_det_Sigma) - 0.5*(x-mu)' *(Sigma \ (x-mu));
z = z + 0.5 * (log_det_Cw) + 0.5*x' * (Cw \ x);

end
function [window1, window2] = normalize(window1, window2)
window1 = window1 - mean(window1);
window2 = window2 - mean(window2);
sigma1 = std(window1);
sigma2 = std(window2);
max_spike1 = max(abs(window1));
max_spike2 = max(abs(window1));

% Step 2: Scale both channels to have the same maximum spike magnitude
desired_spike_magnitude = 1; % Use the larger spike as reference
scaled_x1 = window1 / max_spike1 * desired_spike_magnitude;
scaled_x2 = window2 / max_spike2 * desired_spike_magnitude;

% Step 3: Normalize the variance to 1
window1 = scaled_x1 / std(scaled_x1);
window2 = scaled_x2 / std(scaled_x2);


end