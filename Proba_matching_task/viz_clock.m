function [] = viz_clock(clockwatch)

% This function display lotteries.\
% Input is a clockwatch structure with the following fields

% clockwatch.coord =[X,Y] (e.g. O,O)
% clockwatch.width =[heigh,width] (e.g. 200,50)
% clockwatch.proba = [clockwatch.proba] (e.g 0.25)
% clockwatch.spin = clockwatch.spin angle if any
% clockwatch.text = display test for lotttery carac (0/1)
% clockwatch.ptext_coord = [X,Y] for text for clockwatch proba
% clockwatch.win = [], 0 or 1 (draw a red (0 -loose) or green (1 - win) box, with text)

if ~isempty(clockwatch.coord)
    % display clockwatch caracs as text
    cgtext('Please wait...',clockwatch.text_coord(1),clockwatch.text_coord(2));
    % Display wheel areas, and manage the 90 angle
    cgpencol(clockwatch.color) % winning area
    cgarc(clockwatch.coord(1),clockwatch.coord(2),clockwatch.width(1),clockwatch.width(1),0,360,'S')
    
    for k_tick = 1:12;
        ang = k_tick *2*pi/12;
        cgpencol(0,0,0)
        cgpenwid(3)
        cgdraw(clockwatch.coord(1)-(clockwatch.width(1)-10)*cos(ang)/2,...
            clockwatch.coord(2)+(clockwatch.width(1)-10)*sin(ang)/2,...
            clockwatch.coord(1)-clockwatch.width(1)*cos(ang)/2,...
            clockwatch.coord(2)+clockwatch.width(1)*sin(ang)/2)
    end
    
    
end
cgpencol(0,0,0)
cgpenwid(3)
cgarc(clockwatch.coord(1),clockwatch.coord(2),clockwatch.width(1),clockwatch.width(1),0,0)

% Display the clockwatch hand (refereshing will make it spin !)
if clockwatch.spin
    cgpencol(0,0,0)
    cgpenwid(3)
    cgdraw(clockwatch.coord(1),...
        clockwatch.coord(2),...
        clockwatch.coord(1)-clockwatch.width(1)*cos(-(clockwatch.spin+pi/2))/2,...
        clockwatch.coord(2)+clockwatch.width(1)*sin((-clockwatch.spin+pi/2))/2)
end

% % Display clockwatch result
% X = 0; Y = -50;
% W = 400; H = 200;
% X1 = X - W./2; X2 = X + W./2;
% Y1 = Y - H./2; Y2 = Y + H./2;
% if clockwatch.win==0
%     cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
%      cgpenwid(10);cgpencol(1,0,0);
%     cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
%      cgpenwid(3)
%     cgpencol(1,0,0)
%     cgtext('ANSWER INCORRECT !',0,Y+25)
%     cgtext(clockwatch.FBloss{clockwatch.inc},0,Y-25);
% elseif clockwatch.win==1
%     cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
%      cgpenwid(10);cgpencol(0,1,0);
%     cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
%      cgpenwid(3);cgpencol(0,1,0)
%     cgtext('ANSWER CORRECT !!',0,Y+25)
%     cgtext(clockwatch.FBwin{clockwatch.inc},0,Y-25);
% end


end

