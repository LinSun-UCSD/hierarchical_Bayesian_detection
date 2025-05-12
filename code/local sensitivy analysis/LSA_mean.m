clear all;
close all;
clc;

%%

m = mean(data);
c = cov(data);
res = zeros(5, 100, 2);
for ii = 1:5
    val = m(ii);
    lower = val - 1.5 * val;
    upper = val + 1.5 * val;
    range = linspace(lower, upper, 100);
    for jj = 1:length(range);
        mean_temp = m;
        mean_temp(ii) = range(jj);
        delta = [-abs(val) abs(val)]/10000;
        for kk = 1:2
            n = 1000; % Set the number of samples you want to generate
            samples = mvnrnd(mean_temp+delta(kk), c, n);
            theta = zeros(n, n, 2);
            for i = 1:n
                mu = samples(i,1:2);
                sigma = samples(i, 3:4);
                rho = samples(i, 5);
                C = [sigma(1) ^ 2 sigma(1)*sigma(2)*rho;
                    sigma(1)*sigma(2)*rho sigma(2) ^ 2];
                eigVals = eig(C); % Compute eigenvalues

                if all(eigVals >= 0)
                    temp = mvnrnd(mu, C, n);
                else
                    [V, D] = eig(C);  % Eigen decomposition
                    D(D < 0) = 1e-6;         % Set negative eigenvalues to 0
                    Sigma_psd = V * D * V'; % Reconstruct PSD matrix
                    % Sigma_psd = C + 1e-6 * eye(size(C));
                    temp = mvnrnd(mu, Sigma_psd, n);
                end

                %
                theta(i, :, :) = temp;
            end
            theta = reshape(theta,[], 2);

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

            pfa = 0.000001;
            N = 1000000;
            x = mvnrnd([0; 0], Cw, N);
            temp_z = findZ(X1, X2, z, x);
            z_sorted = sort(temp_z);
            cdf_values = (1:N) / N;
            temp_pfas = 1-cdf_values;
            index = round((1-pfa) * N);
            z_value = z_sorted(index);
            res(ii, jj, kk) = z_value;
        end

    end
end
save("mean_LSA.mat", "res")

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