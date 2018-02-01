% for quick visualization of subject data
% close all

% should open subject mat, and will have sub_data with 
% column 1: offer
% column 2: accepted
% column 3: payment
% column 4: social
% column 5: opponent

% standard error function
ste = @(x) (std(x))/(sqrt(length(x)));


figure
% % plot_titles = {'SocImp', 'NonSocImp', 'SocExp', 'NonSocExp'};
plot_leg = {'0 endow', '10 endow', '20 endow'};
all_leg = {'r-^', 'b-.o', 'k:s'};

% %
% % for pp = 1:4
% %
% %     subplot(2, 2, pp)

    hold on
for ii = 0:2
    
    %     plot(1:ntr, mean(squeeze((ugm(:, 2, ii, :)))), all_leg{ii}, 'markers',5, 'linewidth', 1)
    
    %         errorbar(mean(squeeze((ugm(:, pp, ii, :)))), ste(squeeze((ugm(:, pp, ii, :)))), all_leg{ii}, 'markers',8, 'linewidth', 1)
    
%     ugm_tmp = squeeze((ugm(:, pp, ii, :))); % subject X trial matrix
    
% %     
% %     umeans = [];
% %     ustes  = [];
% %     for nn = 1:size(ugm_tmp, 2)
% %         ugm_clean = reshape(ugm_tmp(:,nn), [], 1);
% %         
% %         ugm_clean(any(isnan(ugm_clean), 2), :) = [];
% %         
% %         if ~isempty(ugm_clean)
% %             umeans(nn) = mean(ugm_clean);
% %             ustes(nn)  = ste(ugm_clean);
% %         end
% %     end
    
% %     errorbar(umeans, ustes, all_leg{ii}, 'markers',8, 'linewidth', 1)
    
    tmp = (sub_data(sub_data(:, 5) == ii, 1));
    
    plot((tmp), all_leg{ii+1}, 'markers',8, 'linewidth', 1)
    
end
ylim([0 15])
xlim([0 24])
legend(plot_leg, 'location', 'northwest')
% %     title(plot_titles{pp})
    hold off
% % end
