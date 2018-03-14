% for her to make errorbars
close all

all_means = [7.96, 9.77, 9.15];
all_se    = [0.051, 0.049, 0.049];
mycolor = {'r', 'g', 'b'};
figure
for ii = 1:3    
    hold on
    bar(ii, all_means(ii), mycolor{ii});
    
end

h = errorbar(1:3,all_means,all_se, 'k.', 'linewidth', 1);
h.CapSize = 25;

ylim([7 10.5]);
xlabel('Conditions');
ylabel('Offers');

legend({'DG', 'UG social','UG non-social'})
set(gca, 'XTick', 1:3)
set(gca,'xticklabel',{'DG', 'UG social','UG non-social'})

hold off

% h=errorbar(x,y,err,'ro');
% h.CapSize = 12;