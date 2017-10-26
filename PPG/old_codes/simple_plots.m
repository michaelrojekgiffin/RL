% this script is only used to make some pretty plots for the expository
% ppt.
clear
close all
% define functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));


% define parameters
priors = [-9.5, 2];
offers = 0:10;

plot(offers,logitp(priors(1:2), offers),'r', 'Linewidth', 1.5)
% xlabel('Offer (euro)')
% ylabel('P(Success)')

% hold on;
% plot(x1, y1, '--k', 'Linewidth', 1.5)
% hold on;
% plot(x2, y2, '--k', 'Linewidth', 1.5)



%% Using own built in function to draw arrow
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} ); 
x   = [5, 5];                       % both x1 and x2 coordinates
y   = [0,logitp(priors(1:2), 5)];   % both y1 and y2 coordinates

x1   = [5, 5];                       % both x1 and x2 coordinates
y1   = [0,logitp(priors(1:2), 5)];   % both y1 and y2 coordinates

x2   = [3, 5.5];
y2   = [logitp(priors(1:2), x2(1)+.5), logitp(priors(1:2), x2(2)+.5)];

plot(offers,logitp(priors(1:2), offers),'r', 'Linewidth', 1.5)
hold on
drawArrow(x,y, 'linewidth', 1,'color','k', 'MaxHeadSize',.5); 
% hold on
% drawArrow(x2,y2,'linewidth',3,'color','r')




%% Using annotation for arrows
% qpt1   = [5/offers(end), 0/1];
% qpt2   = [5/offers(end), logitp(priors(1:2), 5)/1];
% dpt    = qpt2 - qpt1;
x   = [mean(offers)/offers(end), 5/offers(end)];
y   = [0/1,logitp(priors(1:2), 5)/1];

plot(offers,logitp(priors(1:2), offers),'r', 'Linewidth', 1.5)

ha          = annotation('doublearrow');    % store the arrow information in ha
ha.Parent   = gca;                          % associate the arrow the the current axes
ha.X        = [5 5];                        % the location in data units of x1 and x2
ha.Y        = [logitp(priors(1:2), 5) 0];   % the location in data units of y1 and y2

hb          = annotation('doublearrow');   
hb.Parent   = gca;                         
hb.X        = [3.5 5.5];                   
hb.Y        = [logitp(priors(1:2), 3.5) logitp(priors(1:2), 5.5)];   


% hold on
% quiver(qpt1(1), qpt1(2), dpt(1), dpt(2), 'MaxHeadSize', 3)
% legend('acceptance function \phi')
xlabel('Offer (euro)')
ylabel('P(Success)')
