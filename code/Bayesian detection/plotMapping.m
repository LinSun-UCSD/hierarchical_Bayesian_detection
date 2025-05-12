clc;
close all;
clear all;
Cw = [1 0; 0 1];
Ctheta = [   1.6861   -0.0197;
   -0.0197    0.0880];
mu = [10.8202; 0.8962];
detCw = det(Cw);                  % Determinant of covariance matrix
invCw = inv(Cw);  

x1_range = -20:2:20;             % Range for x1
x2_range = -20:2:20;             % Range for x2
[x1Grid, x2Grid] = meshgrid(x1_range, x2_range);  % Create a grid

% Evaluate the region function on the grid
region_values = zeros(size(x1Grid));  % Preallocate for efficiency
for i = 1:size(x1Grid, 1)
    for j = 1:size(x1Grid, 2)
        % Construct x as a 2D vector [x1; x2]
        x = [x1Grid(i, j); x2Grid(i, j)];
        % Compute the region function value
        region_values(i, j) = getMapping(x, mu, Cw, Ctheta);
    end
end
figure;
surf(x1Grid, x2Grid, region_values);  % 3D surface plot
xlabel('x_1');
ylabel('x_2');
zlabel('z');
box on;
% title('Mapping between x_1, x_2 and z');
% shading interp;  % Smooth shading
set(gca, "FontName", "Times New Roman", "FontSize", 10);
print_plot("1.png", 4, 3, 800)
%%
function z = getMapping(x, mu, Cw, Ctheta)

    Sigma = Cw + Ctheta;  % Combined covariance
    log_det_Sigma = log(det(Sigma));
    log_det_Cw = log(det(Cw));
    z = -0.5 * (log_det_Sigma) - 0.5*(x-mu)' *(Sigma \ (x-mu));
    z = z + 0.5 * (log_det_Cw) + 0.5*x' * (Cw \ x);

end