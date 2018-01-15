function [] = viz_gauge(gauge)

% This function display rating_gauge.
% Input is a gauge structure with the following fields

% gauge.coord =[X,Y] (e.g. O,O)
% gauge.width =[width] (e.g. 200)
% gauge.minmax = [minNumber maxnumber] (e.g 0.25)
% gauge.Mgrad = [number of Major gradation] i.e. with number
% gauge.mgrad = [number of minor gradation] i.e. without number
% gauge.color = [R,G,B] color of the gauge
% gauge.curs = position of the cursor
% gauge.curscol = [R,G,B] color of the cursor


% Get general coord
Y= gauge.coord(2);
Xmin = gauge.coord(1)-gauge.width(1)./2; Xmax = gauge.coord(1)+gauge.width(1)./2;
YCmin = gauge.coord(2)+ 50; YCmax = gauge.coord(2)+ 100;

% FILL THE CONFIDENCE AND ADD THE LOTTERY GAUGE
%================================================%
if ~isempty(gauge.lot)
    
    % Filling Gauges
    xconf = ((gauge.conf-.5)./.5).*gauge.width;
    cgpencol(.9,.9,.9); cgrect(Xmin+xconf./2, YCmin+25 , xconf, 50)
    % recalling Confidence proba
    cgpencol(gauge.color)
    cgfont('Arial',38); cgtext(strcat('Your confidence: '),gauge.coord(1)-65,gauge.coord(2)+135)
    cgtext(strcat(num2str(gauge.conf*100),'%'),gauge.coord(1)+100,gauge.coord(2)+135)
    
    % Drawing lottery Gauge
    YLmin = gauge.coord(2)- 50; YLmax = gauge.coord(2)- 100;
    xlot = ((gauge.lot-.5)./.5).*gauge.width;
    % Filling Lottery Gauges
    cgpencol(.3,.3,.3); cgrect(Xmin+xlot./2, YLmin-25 , xlot, 50)
    % Add lottery cursor
    cgpenwid(5);cgpencol(0,0,1)
    cgdraw(Xmin+xlot,YLmin,Xmin+xlot,YLmax);
    % drawing gauge outline
    cgpenwid(3);cgpencol(gauge.color)
    cgdraw(Xmin,YLmin,Xmax,YLmin); cgdraw(Xmin,YLmax,Xmax,YLmax); cgdraw(Xmin,YLmin,Xmin,YLmax); cgdraw(Xmax,YLmin,Xmax,YLmax);
    % recalling Lottery proba
    cgpencol(gauge.color)
    cgfont('Arial',38); cgtext(strcat('The lottery: '),gauge.coord(1)-45,gauge.coord(2)-135)
    cgtext(strcat(num2str(gauge.lot*100),'%'),gauge.coord(1)+80,gauge.coord(2)-135)
    
end
%=============================%

% Add cursor
cgpenwid(5);cgpencol(gauge.curscol)
cgdraw(gauge.curs,YCmin,gauge.curs,YCmax);

% ADD THE CONFIDENCE GAUGE
%=============================%
% drawing gauge outline
cgpenwid(3); cgpencol(gauge.color)
cgdraw(Xmin,YCmin,Xmax,YCmin); cgdraw(Xmin,YCmax,Xmax,YCmax); cgdraw(Xmin,YCmin,Xmin,YCmax); cgdraw(Xmax,YCmin,Xmax,YCmax);
cgpencol(gauge.color)
cgfont('Arial',38); cgtext(strcat('Your confidence: '),gauge.coord(1)-65,gauge.coord(2)+135)
%=============================%

% ADD THE GRADATIONS
%=============================%
% Add big Gradation
for k_grad = 1:gauge.Mgrad
    
    numb = gauge.minmax(1) + (k_grad - 1)* (gauge.minmax(2) - gauge.minmax(1))./(gauge.Mgrad-1);
    Xpos = Xmin(1) + (k_grad - 1)* round(gauge.width./(gauge.Mgrad-1));
    
    cgpenwid(3);cgpencol(gauge.color);cgdraw(Xpos,YCmin,Xpos,YCmin-10);        % conf grad
    if ~isempty(gauge.lot); cgdraw(Xpos,YLmin,Xpos,YLmin+10); end               % lot grad
    cgfont('Arial',38); cgtext(strcat(num2str(numb),'%'),Xpos,gauge.coord(2))   % numbers
    
end

% Add small Gradation
for k_grad = 1:gauge.mgrad
    Xpos = Xmin(1) + (k_grad - 1)* round(gauge.width./(gauge.mgrad-1));
    cgpenwid(2);cgpencol(gauge.color); cgdraw(Xpos,YCmin,Xpos,YCmin-5);         % conf grad
    if ~isempty(gauge.lot);cgdraw(Xpos,YLmin,Xpos,YLmin+5);end;                 % lot grad
end
%=============================%


if ~isempty(gauge.type)
    cgfont('Arial',38);
    cgpenwid(5);
    XXmin = Xmin - 20; XXmax = Xmax+20;
    YYmin = Y + gauge.type*75 - 40 ; YYmax =  Y + gauge.type*75 + 40;   
    if gauge.type == 1
        cgpencol(gauge.curscol);
        cgtext(strcat('Your confidence: '),gauge.coord(1)-65,gauge.coord(2)+135)
        cgtext(strcat(num2str(gauge.conf*100),'%'),gauge.coord(1)+100,gauge.coord(2)+135)
        cgdraw(XXmin,YYmin,XXmax,YYmin); cgdraw(XXmin,YYmax,XXmax,YYmax); cgdraw(XXmin,YYmin,XXmin,YYmax); cgdraw(XXmax,YYmin,XXmax,YYmax);
    elseif gauge.type == -1
        cgpencol(0,0,1);
        cgtext(strcat('The lottery: '),gauge.coord(1)-45,gauge.coord(2)-135)
        cgtext(strcat(num2str(gauge.lot*100),'%'),gauge.coord(1)+80,gauge.coord(2)-135)
        cgdraw(XXmin,YYmin,XXmax,YYmin); cgdraw(XXmin,YYmax,XXmax,YYmax); cgdraw(XXmin,YYmin,XXmin,YYmax); cgdraw(XXmax,YYmin,XXmax,YYmax);
    end
end

cgpencol(gauge.color)
