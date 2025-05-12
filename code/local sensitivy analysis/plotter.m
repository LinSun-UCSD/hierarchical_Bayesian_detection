clc;
clear all;
close all;
%%

m = mean(data2);
c = cov(data2);
colors = {'k', 'b', 'r','c','green'};
figure();
for ii = 1:5
    val = m(ii);
    lower = val - 1.5 * val;
    upper = val + 1.5 * val;
    range = linspace(lower, upper, 100);
    subplot(5,1,ii)
    plot(range, data1(ii, :, 2) - data1(ii, :, 1), 'color', colors{ii}, "LineWidth", 1.2);
    xlim([min(lower, upper), max(lower, upper)]);
    ylabel("LSA");
    xlabel(['\mu_' num2str(ii)])
    set(gca, "FontSize", 10, "FontName", "Times New Roman");
    grid on;
    ylim([-4 6])
end
% print_plot("1.png", 5, 6, 800)
%%

m = mean(data2);
c = cov(data2);
for ii = 1:5
    val = sqrt(c(ii,ii));
    lower = val * 0.1;
    upper = val * 10;
    range = linspace(lower, upper, 100);
    subplot(5,1,ii)
    plot(range, data1(ii, :, 2) - data1(ii, :, 1), 'color', colors{ii}, "LineWidth", 1.2);
    xlim([min(lower, upper), max(lower, upper)]);
    ylabel("LSA");
    xlabel(['\sigma_' num2str(ii)])
    set(gca, "FontSize", 10, "FontName", "Times New Roman");
    grid on;
    ylim([-4 6])
end
% print_plot("2.png", 5, 6, 800)
%%

close all;
figure();
count = 1;
for ii = 1:5
    for mm = ii+1:5
        range = linspace(-0.99, 0.99, 100);
        if count>5
        subplot(5,1,count-5)
        temp2 = reshape(data1(ii,mm,:,2), 1, []);
        temp1 = reshape(data1(ii,mm,:,1), 1, []);
        plot(range, temp2- temp1, 'Color', 'k', "LineWidth", 1.2)
        
        ylabel("LSA");
        xlabel(['\rho_' num2str(ii) '_' num2str(mm)])
        set(gca, "FontSize", 10, "FontName", "Times New Roman");
        grid on;
        ylim([-4 6])
        end
        count = count + 1 ;
    end
end
print_plot("1.png", 5, 6, 800)