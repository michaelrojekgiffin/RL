function [] = viz_wheels(lottery)

% This function display lotteries.\
% Input is a lottery structure with the following fields

% lottery.coord =[X,Y] (e.g. O,O)
% lottery.width =[heigh,width] (e.g. 200,50)
% lottery.proba = [lottery.proba] (e.g 0.25)
% lottery.spin = lottery.spin angle if any
% lottery.text = display test for lotttery carac (0/1)
% lottery.ptext_coord = [X,Y] for text for lottery proba
% lottery.win = [], 0 or 1 (draw a red (0 -loose) or green (1 - win) box, with text)

X=90-lottery.proba*360;

if ~isempty(lottery.coord)
    
    % display lottery caracs as text
    cgtext(strcat([num2str(lottery.proba*100) ' %']),lottery.ptext_coord(1),lottery.ptext_coord(2));
    
    % Display wheel areas, and manage the 90 angle
    if lottery.proba <0.25
        cgpencol(lottery.winning_color) % winning area
        cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),X,90,'S')
        cgpencol(lottery.losing_color) % losing area
        cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),90,360,'S')
        cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),0,X,'S')
    elseif lottery.proba <1
        cgpencol(lottery.winning_color) % winning area
        cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),0,90,'S')
        cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),X,360,'S')
        cgpencol(lottery.losing_color) % losing area
        cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),90,X,'S')
    else
        cgpencol(lottery.winning_color) % winning area
        cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),90,90,'S')
    end
    cgpencol(0,0,0)
    cgpenwid(3)
    cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),0,0)
    
    % Display the lottery hand (refereshing will make it spin !)
    if lottery.spin
        cgpencol(0,0,0)
        cgpenwid(3)
        cgdraw(lottery.coord(1),...
            lottery.coord(2),...
            lottery.coord(1)-lottery.width(1)*cos(-(lottery.spin+pi/2))/2,...
            lottery.coord(2)+lottery.width(1)*sin((-lottery.spin+pi/2))/2)
    end
    
%     % Display Lottery result
%     X = 0; Y = -50;
%     W = 400; H = 200;
%     X1 = X - W./2; X2 = X + W./2;
%     Y1 = Y - H./2; Y2 = Y + H./2;
%     if lottery.win==0
%         cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
%         cgpenwid(10); cgpencol(1,0,0);
%         cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
%          cgpenwid(3); cgpencol(1,0,0)
%         cgtext('LOSING LOTTERY!',0,Y+25)
%         cgtext(lottery.FBloss{lottery.inc},0,Y-25);
%     elseif lottery.win==1
%         cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
%          cgpenwid(10); cgpencol(0,1,0);
%         cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
%          cgpenwid(3); cgpencol(0,1,0)
%         cgtext('WINNING LOTTERY!',0,Y+25)
%         cgtext(lottery.FBwin{lottery.inc},0,Y-25);
%     end
    
    
end

