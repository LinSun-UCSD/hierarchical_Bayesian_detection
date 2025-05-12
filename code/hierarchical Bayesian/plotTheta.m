close all;
clear all;
clc;


%% read the directory folders


% Get the list of all items in the directory
allItems = dir(path);

% Filter only directories (ignoring files)
isFolder = [allItems.isdir];  % Logical array indicating folders
folderNames = {allItems(isFolder).name};  % Get names of folders

% Remove '.' and '..' (current and parent directories)
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

% Display the folder names
disp('Folders:');
disp(folderNames);
%%
colors = {'k', 'b', 'r','c','green'};
legendNames = {"May", "June", "July", "Sept.", "Nov."};
monthlyValuesMean = cell(1, 12); % 12 months
monthlyValuesSTD = cell(1, 12); % 12 months
figure()
for i = 1:5
    norm_window1_mean = [];
    norm_window1_std = [];
    norm_window2_mean = [];
    norm_window2_std = [];
    folderPath = [path '\' folderNames{i}];
    % Get a list of all .mat files in the folder
    matFiles = dir(fullfile(folderPath, '*.mat'));
    
    % Loop through each file and load the data
    for k = 1:length(matFiles)
        % Get the full path of the .mat file
        matFilePath = fullfile(folderPath, matFiles(k).name);

        % Load the .mat file
        data = load(matFilePath);
        dataToStore = data.dataToStore;
        datetimeString = dataToStore.time_stamp;

        % Convert to datetime object
        dt = datetime(datetimeString, 'InputFormat', 'yyyy-MM-dd''T''HH_mm_ss');
        
        % Get the month name
        monthNumber = month(dt);

        % Display the result
        disp(['Month Number: ', num2str(monthNumber)]);
        
       
        if strcmp(dataToStore.event_label, 'full-open')
            for l = 1:length(dataToStore.spikeTime)
                clear temp;
                if iscell(dataToStore.data1)
                    if length(dataToStore.spikeTime) > 1
                        temp1 = dataToStore.data1(1,l);
                        temp1 = temp1{1};
                        temp2 = dataToStore.data2(1,l);
                        temp2 = temp2{1};
                    else
                        temp1 = dataToStore.data1;
                        temp2 = dataToStore.data2;
                    end
                    positiveValues = temp1(temp1 > 0);
                    if length(positiveValues) == 0
                        continue;
                    end
                    % norm_window1_mean = [norm_window1_mean; mean(positiveValues)];
                    % norm_window1_std = [norm_window1_std; std(positiveValues)];
                    monthlyValuesSTD{monthNumber} = [monthlyValuesSTD{monthNumber}, std(positiveValues)];
                    monthlyValuesMean{monthNumber} = [monthlyValuesMean{monthNumber}, mean(positiveValues)];
                    % positiveValues = temp2(temp2 > 0);
                    % norm_window2_mean = [norm_window2_mean; mean(positiveValues)];
                  
                end
            end
        end
    end


end


%%

colors = {'k', 'b', 'r','c','green'};
legendNames = {"May", "June", "July", "Sept.", "Nov."};
monthlyValuesOpen = zeros(12,1); 
monthlyValuesTotal = zeros(12,1);
figure()
for i = 1:5
    norm_window1_mean = [];
    norm_window1_std = [];
    norm_window2_mean = [];
    norm_window2_std = [];
    folderPath = [path '\' folderNames{i}];
    % Get a list of all .mat files in the folder
    matFiles = dir(fullfile(folderPath, '*.mat'));
    
    % Loop through each file and load the data
    for k = 1:length(matFiles)
        % Get the full path of the .mat file
        matFilePath = fullfile(folderPath, matFiles(k).name);

        % Load the .mat file
        data = load(matFilePath);
        dataToStore = data.dataToStore;
        datetimeString = dataToStore.time_stamp;

        % Convert to datetime object
        dt = datetime(datetimeString, 'InputFormat', 'yyyy-MM-dd''T''HH_mm_ss');
        
        % Get the month name
        monthNumber = month(dt);
        if strcmp(dataToStore.event_label, 'full-open') || strcmp(dataToStore.event_label, 'open') 
            monthlyValuesOpen(monthNumber) = monthlyValuesOpen(monthNumber)  + 1;

        end
         monthlyValuesTotal(monthNumber) = monthlyValuesTotal(monthNumber)  + 1;
    end


end 
save("mean.mat", "monthlyValuesMean");
save("std.mat", "monthlyValuesSTD")
%%
path = 'D:\UCSD Post-doc\matlabCode\collect statistics\data';

% Get the list of all items in the directory
allItems = dir(path);

% Filter only directories (ignoring files)
isFolder = [allItems.isdir];  % Logical array indicating folders
folderNames = {allItems(isFolder).name};  % Get names of folders

% Remove '.' and '..' (current and parent directories)
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

% Display the folder names
disp('Folders:');
disp(folderNames);
%%
figure()
for i = 1
    subplot(2,3,i)
    norm_window1_mean = [];
    norm_window1_std = [];
    norm_window2_mean = [];
    norm_window2_std = [];
    folderPath = [path '\' folderNames{i}];
    % Get a list of all .mat files in the folder
    matFiles = dir(fullfile(folderPath, '*.mat'));
    
    % Loop through each file and load the data
    for k = 1:length(matFiles)
        % Get the full path of the .mat file
        matFilePath = fullfile(folderPath, matFiles(k).name);

        % Load the .mat file
        data = load(matFilePath);
        dataToStore = data.dataToStore;

        if strcmp(dataToStore.event_label, 'full-open')
            for l = 1:length(dataToStore.spikeTime)
                clear temp;
                if iscell(dataToStore.data1)
                    if length(dataToStore.spikeTime) > 1
                        temp1 = dataToStore.data1(1,l);
                        temp1 = temp1{1};
                        temp2 = dataToStore.data2(1,l);
                        temp2 = temp2{1};
                    else
                        temp1 = dataToStore.data1;
                        temp2 = dataToStore.data2;
                    end
                    positiveValues = temp1(temp1 > 0);
                    if length(positiveValues) == 0
                        continue;
                    end
                    norm_window1_mean = [norm_window1_mean; mean(positiveValues)];
                    norm_window1_std = [norm_window1_std; std(positiveValues)];
 
                    positiveValues = temp2(temp2 > 0);
                    norm_window2_mean = [norm_window2_mean; mean(positiveValues)];
                    norm_window2_std = [norm_window2_std; std(positiveValues)];
                end
            end
        end


        % You can now process the loaded data as needed
    end

end
% print_plot("1.png",7,4.5,800)
%%


%%
%% plot the correlation coefficient


for i = 1

    correlation = [];
    folderPath = [path '\' folderNames{i}];
    % Get a list of all .mat files in the folder
    matFiles = dir(fullfile(folderPath, '*.mat'));

    % Loop through each file and load the data
    for k = 1:length(matFiles)
        % Get the full path of the .mat file
        matFilePath = fullfile(folderPath, matFiles(k).name);

        % Load the .mat file
        data = load(matFilePath);
        dataToStore = data.dataToStore;

        if strcmp(dataToStore.event_label, 'full-open')
            for l = 1:length(dataToStore.spikeTime)
                clear temp;
                if iscell(dataToStore.data1)
                    if length(dataToStore.spikeTime) > 1
                        temp1 = dataToStore.data1(1,l);
                        temp1 = temp1{1};
                        temp2 = dataToStore.data2(1,l);
                        temp2 = temp2{1};
                    else
                        temp1 = dataToStore.data1;
                        temp2 = dataToStore.data2;
                    end
                    % positiveValues1 = temp1(temp1>0);
                    if length(temp1) == 0
                        continue;
                    end
                   

                    % positiveValues2 = temp2(temp2>0);
                    correlation = [correlation; corr(temp1', temp2')];
                end
            end
        end
    end


  
end
%%
close all;
figure()
subplot(1,3,1)
scatter(norm_window1_mean, norm_window2_mean);
subplot(1,3,2)
scatter(norm_window1_std, norm_window2_std);
subplot(1,3,3)
histogram(correlation);
%%
data = [norm_window1_mean norm_window2_mean norm_window1_std norm_window2_std correlation];
data(any(isnan(data), 2), :) = [];
save("data.mat", "data")
%%
c = cov(data);
m = mean(data);
n = 1000; % Set the number of samples you want to generate
samples = mvnrnd(m, c, n);
theta = zeros(n, n, 2);
for i = 1:n
    mu = samples(i,1:2);
    sigma = samples(i, 3:4);
    rho = samples(i, 5);
    C = [sigma(1) ^ 2 sigma(1)*sigma(2)*rho;
         sigma(1)*sigma(2)*rho sigma(2) ^ 2];
    temp = mvnrnd(mu, C, n);
    theta(i, :, :) = temp;
end
theta = reshape(theta,[], 2);
%%
save("theta.mat", "theta")
%%
figure()
scatter(theta(:,1), theta(:,2));

figure()
histogram(theta(:,1))
%%
close all;
% Create a grid for the contour plot
x = linspace(min(theta(:,1)), max(theta(:,1)), 100);
y = linspace(min(theta(:,2)), max(theta(:,2)), 100);
[X, Y] = meshgrid(x, y);

% Estimate the density using kernel density estimation
data = [theta(:,1), theta(:,2)];
density = mvksdensity(data, [X(:), Y(:)]); % MATLAB's multivariate KDE function
density = reshape(density, size(X)); % Reshape density to match the grid

% Plot the scatter plot
scatter(theta(:,1), theta(:,2), 10,'k'); % Scatter plot of points
hold on;
box on;
% Overlay the contour plot
contour(X, Y, density, 20, 'LineWidth', 1.5); % Contour plot with 20 levels
hold off;
grid on;
xlabel("\theta_1");
ylabel("\theta_2")
set(gca,'FontSize',10, "FontName", "Times New Roman");
print_plot("1.png",4,3,800)
%%

%%
