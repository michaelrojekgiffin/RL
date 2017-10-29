% this script gets the responses we collected from the responders and
% extracts their parameters, i.e. the slope and intercepts of their logit
% response functions

% clear all
close all
% clc
logitp          = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

offers          = 0:20;

cur_dir         = pwd;
txtnm           = fullfile(cur_dir,'ultimatum_responders.xlsx');
[num,txt,raw]   = xlsread(txtnm);

resp_params     = NaN(5, 2);

endow_freq      = NaN(5, 20);
endow_count     = 0;
figure
for ee = 1:5 % because we have 5 different conditions going up to 20 in steps of 5
    endow_lines = find(num(:,6)==endow_count); % get rows for the conditions (specified by starting endowment)
    for oo = offers % get frequency of each offer being accepted
        endow_freq(ee, oo+1)  = sum(num(endow_lines,4) == 1 & num(endow_lines,7) == oo) / sum(num(endow_lines,7) == oo);
    end
    subplot(5, 1, ee)
    
    endow_sig = endow_freq(ee, :);
    for ss = 2:length(endow_sig)
        if endow_sig(ss) < endow_sig(ss-1)
            endow_sig(ss) = endow_sig(ss-1);
        end
    end
    [xData, yData]      = prepareCurveData( offers, endow_sig);
    
    % Set up fittype and options.
    % ft = fittype( 'a/(1+exp(-b*x))', 'independent', 'x', 'dependent', 'y' );
    ft                  = fittype( 'exp(intercept+slope.*(x))/(1+exp(intercept+slope.*(x)))', 'independent', 'x', 'dependent', 'y' );
    
    opts                = fitoptions( 'Method', 'NonlinearLeastSquares' );
    % opts.Display = 'On';
    opts.StartPoint     = [-15 5];
    
    % Fit model to data.
    [fitresult, ~]      = fit( xData, yData, ft, opts );
    resp_params(ee, 1)  = fitresult.intercept;
    resp_params(ee, 2)  = fitresult.slope;
    
    
    plot(logitp([fitresult.intercept,fitresult.slope],offers),'-r', 'Linewidth', 1.5);
    %     plot(offers, endow_freq(ee, :), '-r', 'Linewidth', 1.5);
        hold on
    bar(endow_freq(ee, :), 'b');
    xlim([offers(1), offers(end)]);
    title(sprintf('%d euro starting endowment', endow_count));
    xlabel('Proposer offer');
    ylabel('Acc freq');
    hold off
    endow_count     =   endow_count+5;
end

save('resp_params.mat', 'resp_params');
% %%
% [xData, yData] = prepareCurveData( offers, endow_freq(1, :) );
% logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
% % Set up fittype and options.
% % ft = fittype( 'a/(1+exp(-b*x))', 'independent', 'x', 'dependent', 'y' );
% ft = fittype( 'exp(intercept+slope.*(x))/(1+exp(intercept+slope.*(x)))', 'independent', 'x', 'dependent', 'y' );
% 
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % opts.Display = 'On';
% opts.StartPoint = [-12 2];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% plot(logitp([fitresult.intercept,fitresult.slope],offers));
