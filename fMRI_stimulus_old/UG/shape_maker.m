% making of shapes for ultimatum game
% Author: Michael Giffin, January, 2018

%% Create  Filled Circle
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);
% 
% xlim([1 300]);
% ylim([1 300]);
h = patch;
axis off
colormap(gray)

% Set Coordinates
t = linspace(0,2*pi,120);
h.XData = 256 + 128*cos(t);
h.YData = 256 + 128*sin(t);
h.FaceVertexCData = h.XData;
% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/circle_filled.png')

%% Create Empty Circle
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);

h = line;
axis off
colormap(gray)

% Set Coordinates
t = linspace(0,2*pi,120);
h.XData = 256 + 128*cos(t);
h.YData = 256 + 128*sin(t);
set(h(1),'linewidth',10);

% Save Image
% img = getframe(gcf);
% imwrite(img.cdata,'images/circle_frame.png')

%% Create filled Square
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [100 300 300 100 100];
y = [100 100 300 300 100];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/square_filled.png')


%% Create empty Square
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

% h = patch;
axis off
colormap(gray)

x = [100 300 300 100 100];
y = [100 100 300 300 100];

patch(x, y, 'white', 'linewidth', 10);

% h.XData = x;
% h.YData = y;
% h.FaceVertexCData = h.XData;

% 
% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/square_frame.png')


%% Create filled Triangle
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [100 300 200 100];
y = [100 100 300 100];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/triangle_filled.png')


%% Create empty Triangle
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [100 300 200 100];
y = [100 100 300 100];

patch(x, y, 'white', 'linewidth', 10);

% h.XData = x;
% h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/triangle_frame.png')


%% Create filled Diamond
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

hyp = 200/(sqrt(2)); % hypotenuse or diagonal of square to keep size the same
% hyp = 100;

x = [55 55+hyp 55+hyp+hyp 55+hyp 55];
y = [55+hyp 55 55+hyp 55+hyp+hyp 55+hyp];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/diamond_filled.png')


%% Create empty Diamond
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

hyp = 200/(sqrt(2)); % hypotenuse or diagonal of square to keep size the same

x = [55 55+hyp 55+hyp+hyp 55+hyp 55];
y = [55+hyp 55 55+hyp 55+hyp+hyp 55+hyp];

patch(x, y, 'white', 'linewidth', 10);


% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/diamond_frame.png')


%% Create filled trapazoid
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [160 240 330 70 160];
y = [100 100 300 300 100];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/trapazoid_filled.png')


%% Create empty trapazoid
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [160 240 330 70 160];
y = [100 100 300 300 100];

patch(x, y, 'white', 'linewidth', 10);

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/trapazoid_frame.png')

%% Create filled star
% I had to fill it with a circle, it took longer than it should have to get
% this right
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

xCenter = 200;
yCenter = 200;

h = patch;
hold on
% to fill interior with circle
% Set Coordinates
t = linspace(0,2*pi,120);
h.XData = xCenter + 57*cos(t);
h.YData = yCenter + 57*sin(t);
h.FaceVertexCData = h.XData;


% Define parameters
numberOfPoints = 5;
rotationAngle = 0;
xCenter = 200;
yCenter = 200;
% Determine the angles that the arm tips are at
theta = (0 : (numberOfPoints-1)/numberOfPoints*pi : (numberOfPoints-1)*pi) + rotationAngle;
% Define distance from the arm tip to the center of star.
amplitude = 150;
% Get x and y coordinates of the arm tips.
x = amplitude .* cos(theta) + xCenter;
y = amplitude .* sin(theta) + yCenter;


h = patch;
axis off
colormap(gray)


h.XData = x;
h.YData = y;
h.FaceVertexCData = h.XData;

% % the distance between two points on the star
% ptdis = sqrt(((y(4) - y(1))^2) + ((x(4) - x(1))^2));
% 
% % the angles of the resulting triangle are 36, 36, and 108
% % the length of the edge of one of the stars triangles is
% triedge = (ptdis * sind(36))/(sind(108));
% 
% % and the length of the interior of the star's pentagon
% pentedge = (x(1) - x(2)) - (2*triedge);


% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/star_filled.png')


%% Create empty star
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);


% Define parameters
numberOfPoints = 5;
rotationAngle = 0;
xCenter = 200;
yCenter = 200;
% Determine the angles that the arm tips are at
theta = (0 : (numberOfPoints-1)/numberOfPoints*pi : (numberOfPoints-1)*pi) + rotationAngle;
% Define distance from the arm tip to the center of star.
amplitude = 150;
% Get x and y coordinates of the arm tips.
x = amplitude .* cos(theta) + xCenter;
y = amplitude .* sin(theta) + yCenter;

h = patch;
axis off
colormap(gray)

patch(x, y, 'white', 'linewidth', 10);

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/star_frame.png')



%% Create Filled Clover
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);


xCenter = 256;
yCenter = 256;

% to fill the inner circle space
h = patch;
hold on

axis off
colormap(gray)

% Set Coordinates
t = linspace(0,2*pi,120);

h.XData = xCenter + 45*cos(t);
h.YData = yCenter + 45*sin(t);

h = patch;

% original circle code
% h.XData = 256 + 128*cos(t);
% h.YData = 256 + 128*sin(t);

h.XData = [210 + 80*cos(t), 302 + 80*cos(t)];
h.YData = [210 + 80*sin(t), 302 + 80*sin(t)];

h = patch; 
h.XData = 310 + 80*cos(t);
h.YData = 190 + 80*sin(t);


hold off

% h.FaceVertexCData = [h.XData];

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/clover_filled.png')



%% Create Empty Clover
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);


xCenter = 256;
yCenter = 256;
t = linspace(0,2*pi,120);

h = line;
hold on
% original circle code
% h.XData = 256 + 128*cos(t);
% h.YData = 256 + 128*sin(t);

% h.XData = [210 + 80*cos(t), 302 + 80*cos(t)];
% h.YData = [210 + 80*sin(t), 302 + 80*sin(t)];

h.XData = [210 + 80*cos(t)];
h.YData = [210 + 80*sin(t)];

set(h(1),'linewidth',10);

h = line; 

h.XData = [302 + 80*cos(t)];
h.YData = [302 + 80*sin(t)];

set(h(1),'linewidth',10);

h = line; 
h.XData = 310 + 80*cos(t);
h.YData = 190 + 80*sin(t);

set(h(1),'linewidth',10);

hold off

% h.FaceVertexCData = [h.XData];

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/clover_frame.png')



%% Create Filled helix
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);

xCenter = 256;
yCenter = 256;

h = patch;

axis off
colormap(gray)

% Set Coordinates
t = linspace(0,2*pi,25);

h.XData = [xCenter + 20*t, xCenter - 20*t];
h.YData = [yCenter + 128*sin(t), yCenter - 128*sin(t)];


% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/helix_filled.png')


%% Create Empty helix
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);

xCenter = 256;
yCenter = 256;

h = line;
hold on
axis off
colormap(gray)

% Set Coordinates
t = linspace(0,2*pi,25);

% h.XData = [xCenter + 20*t, xCenter - 20*t];
% h.YData = [yCenter + 128*sin(t), yCenter - 128*sin(t)];

h.XData = [xCenter + 20*t];
h.YData = [yCenter + 128*sin(t)];

set(h(1),'linewidth',10);

h = line;

h.XData = [xCenter - 20*t];
h.YData = [yCenter - 128*sin(t)];
set(h(1),'linewidth',10);

h = line;

h.XData = [xCenter/2 xCenter+xCenter/2];
h.YData = [yCenter yCenter];
set(h(1),'linewidth',10);


% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/helix_frame.png')


%% Create Filled Crescent
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);
% 
% xlim([1 300]);
% ylim([1 300]);
h = patch;
axis off
colormap(gray)

% Set Coordinates
% t = linspace(0,2*pi,120);
t = linspace(0,2*pi,120);

x = 256 + 128*cos(t);
y = 256 - 128*sin(t);

x = x((length(x)/2):end);
y = y((length(y)/2):end);

x = x((length(x)/7):end);
y = y((length(y)/7):end);



h.XData = x;
% h.YData = 256 + 128*sin(t);
h.YData = y;

h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/crescent_filled.png')


%% Create Empty Crescent
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);
% 
% xlim([1 300]);
% ylim([1 300]);
h = line;
axis off
colormap(gray)

% Set Coordinates
% t = linspace(0,2*pi,120);
t = linspace(0,2*pi,120);

x = 256 + 128*cos(t);
y = 256 - 128*sin(t);

x = x((length(x)/2):end);
y = y((length(y)/2):end);

x = x((length(x)/7):end);
y = y((length(y)/7):end);

h.XData = x;
% h.YData = 256 + 128*sin(t);
h.YData = y;

set(h(1),'linewidth',10);

hold on

h = line;
h.XData = [x(1), x(end)];
% h.YData = 256 + 128*sin(t);
h.YData = [y(1) y(end)];
set(h(1),'linewidth',10);


% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/crescent_frame.png')


%% Create filled Tetris F
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)


x = [100 100 150 150 300 300 200 200 150 150 100];
y = [150 200 200 300 300 250 250 100 100 150 150];


h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/tetris_F_filled.png')



%% Create empty Tetris F
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = line;
axis off
colormap(gray)


x = [100 100 150 150 300 300 200 200 150 150 100];
y = [150 200 200 300 300 250 250 100 100 150 150];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

set(h(1),'linewidth',10);

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/tetris_F_frame.png')


%% Create filled bird
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

% x = [100 300 300 100 100];
% y = [100 100 300 300 100];

% x = [100 200 220 275 300 220 200 150 100];
% y = [100 200 180 300 275 100 160 100 100];

x = [100 200 220 275 300 220 200 150 100];
y = [100 200 180 300 275 100 160 100 100];


h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;
% needs to be rotated once to the right in the png file
% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/bird_filled.png')


%% Create empty bird
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = line;
axis off
colormap(gray)

% x = [100 300 300 100 100];
% y = [100 100 300 300 100];

% x = [100 200 220 275 300 220 200 150 100];
% y = [100 200 180 300 275 100 160 100 100];

x = [100 200 220 275 300 220 200 150 100];
y = [100 200 180 300 275 100 160 100 100];

set(h(1),'linewidth',10);


h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;
% needs to be rotated once to the right in the png file
% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/bird_frame.png')



%% Create filled Staircase
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

hold on

h = patch;
axis off
colormap(gray)

% x = [100 300 300 100 100];
% y = [100 100 300 300 100];

x = [100 150 150 100 100];
y = [100 100 150 150 100];

h.XData = x;
h.YData = y;

h = patch;
axis off
colormap(gray)

x = x+50;
% y = y+50;

h.XData = x;
h.YData = y;


h = patch;
axis off
colormap(gray)

% x = x+50;
y = y+50;

h.XData = x;
h.YData = y;


h = patch;
axis off
colormap(gray)

x = x+50;
% y = y+50;

h.XData = x;
h.YData = y;


h = patch;
axis off
colormap(gray)

% x = x+50;
y = y+50;

h.XData = x;
h.YData = y;


h = patch;
axis off
colormap(gray)

% x = [100 300 300 100 100];
% y = [100 100 300 300 100];

x = [200 250 250 200 200];
y = [100 100 150 150 100];

h.XData = x;
h.YData = y;


% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/stairs_filled.png')



%% Create Empty Staircase
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

hold on

h = line;
axis off
colormap(gray)

% x = [100 300 300 100 100];
% y = [100 100 300 300 100];

x = [100 150 150 100 100];
y = [100 100 150 150 100];

h.XData = x;
h.YData = y;
set(h(1),'linewidth',10);

h = line;
axis off
colormap(gray)

x = x+50;
% y = y+50;

h.XData = x;
h.YData = y;
set(h(1),'linewidth',10);

h = line;
axis off
colormap(gray)

% x = x+50;
y = y+50;

h.XData = x;
h.YData = y;
set(h(1),'linewidth',10);

h = line;
axis off
colormap(gray)

x = x+50;
% y = y+50;

h.XData = x;
h.YData = y;
set(h(1),'linewidth',10);

h = line;
axis off
colormap(gray)

% x = x+50;
y = y+50;

h.XData = x;
h.YData = y;
set(h(1),'linewidth',10);

h = line;
axis off
colormap(gray)

% x = [100 300 300 100 100];
% y = [100 100 300 300 100];

x = [200 250 250 200 200];
y = [100 100 150 150 100];

h.XData = x;
h.YData = y;
set(h(1),'linewidth',10);
hold off
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/stairs_frame.png')

