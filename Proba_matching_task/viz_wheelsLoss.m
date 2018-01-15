function [] = viz_wheelsLoss(lottery)

% This function display lotteries.\
% Input is a lottery structure with the following fields

% lottery.coord =[X,Y] (e.g. O,O)
% lottery.width =[heigh,width] (e.g. 200,50)
% lottery.proba = [lottery.proba] (e.g 0.25)
% lottery.value = [lottery.value] (e.g 1)
% lottery.box = logical (0: absent, 1: present)
% lottery.spin = lottery.spin angle if any
% lottery.text = display test for lotttery carac (0/1)
% lottery.vtext_coord = [X,Y] for text for value
% lottery.ptext_coord = [X,Y] for text for value
% lottery.win = [], 0 or 1 (draw a red (0 -loose) or green (1 - win) box, with text)

X=90-lottery.proba*360;

if ~isempty(lottery.coord)
    
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
    cgpencol(1,1,1)
    cgpenwid(3)
    cgarc(lottery.coord(1),lottery.coord(2),lottery.width(1),lottery.width(1),0,0)
    
    % display box around chosen lottery
    cgpencol(1,1,1)
    if lottery.box==1
        cgpenwid(3)
        Ymin = lottery.coord(2)-200;Ymax = lottery.coord(2)+250;
        Xmin = lottery.coord(1)-250;Xmax = lottery.coord(1)+250;
        cgdraw(Xmax,Ymax,Xmin,Ymax)
        cgdraw(Xmax,Ymin,Xmin,Ymin)
        cgdraw(Xmin,Ymax,Xmin,Ymin)
        cgdraw(Xmax,Ymin,Xmax,Ymax)
        
    end
    
    
    % Display the rectangles for both magnitudes
    height_g = lottery.width * lottery.gp;   % height of magnitude 1
    height_l = lottery.width * lottery.lp;   % height of magnitude 2
    width_rec = 40;                            % width for both magnitudes
    coord_rec_g = [lottery.coord(1) + (lottery.width(1)./2 + 50),lottery.coord(2) - lottery.width(1)./2];
    coord_rec_l = [lottery.coord(1) - (lottery.width(1)./2 + 50),lottery.coord(2) - lottery.width(1)./2];
    
    cgpencol(lottery.winning_color) % blue rectangle
    cgrect(coord_rec_g(1), coord_rec_g(2) + height_g / 2, width_rec, height_g)
    cgpencol(1, 1, 1) % white around the rectangle
    cgpenwid(3)
    lc = coord_rec_g(1) - width_rec / 2; % left_cood
    rc = coord_rec_g(1) + width_rec / 2; % right_cood
    dc = coord_rec_g(2); % down_cood
    uc = coord_rec_g(2) + height_g; % up_cood
    cgdraw(lc, dc, lc, uc) % left bar
    cgdraw(rc, dc, rc, uc) % right bar
    cgdraw(lc, dc, rc, dc) % dow bar
    cgdraw(lc, uc, rc, uc) % up bar
    
    cgpencol(lottery.losing_color)% yellow rectangle
    cgrect(coord_rec_l(1), coord_rec_l(2) + height_l / 2, width_rec, height_l)
    cgpencol(1, 1, 1) % white around the rectangle
    cgpenwid(3)
    lc = coord_rec_l(1) - width_rec / 2; % left_cood
    rc = coord_rec_l(1) + width_rec / 2; % right_cood
    dc = coord_rec_l(2); % down_cood
    uc = coord_rec_l(2) + height_l; % up_cood
    cgdraw(lc, dc, lc, uc) % left bar
    cgdraw(rc, dc, rc, uc) % right bar
    cgdraw(lc, dc, rc, dc) % dow bar
    cgdraw(lc, uc, rc, uc) % up bar
    
    % display lottery caracs as text
    if lottery.text==1
        cgpencol(lottery.winning_color)
        cgtext(strcat([num2str(lottery.proba*100) '% chance to win: ' num2str(lottery.gain) ' €']),lottery.gtext_coord(1),lottery.gtext_coord(2));
        cgpencol(lottery.losing_color)
        cgtext(strcat([num2str((1-lottery.proba)*100) '% to lose: ' num2str(lottery.loss) ' €']),lottery.ltext_coord(1),lottery.ltext_coord(2));
        % cgtext(strcat([num2str(lottery.proba*100) ' %']),lottery.ptext_coord(1),lottery.ptext_coord(2));
    end
    
end

