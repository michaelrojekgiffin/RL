
function [] = viz_outcome(trial_type)


% Display Lottery result
X = 0; Y = -50;
W = 400; H = 200;
X1 = X - W./2; X2 = X + W./2;
Y1 = Y - H./2; Y2 = Y + H./2;
if trial_type.win==1
    cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
    cgpenwid(10); cgpencol(1,0,0);
    cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
    cgpenwid(3); cgpencol(1,0,0)
    cgtext('LOSING LOTTERY!',0,Y+25)
    cgtext(trial_type.FBloss{trial_type.inc},0,Y-25);
elseif trial_type.win==2
    cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
    cgpenwid(10); cgpencol(0,1,0);
    cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
    cgpenwid(3); cgpencol(0,1,0)
    cgtext('WINNING LOTTERY!',0,Y+25)
    cgtext(trial_type.FBwin{trial_type.inc},0,Y-25);
elseif trial_type.win==3
    cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
    cgpenwid(10);cgpencol(1,0,0);
    cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
    cgpenwid(3)
    cgpencol(1,0,0)
    cgtext('ANSWER INCORRECT !',0,Y+25)
    cgtext(trial_type.FBloss{trial_type.inc},0,Y-25);
elseif trial_type.win==4
    cgpencol(.3,.3,.3); cgrect(X, Y , W, H)
    cgpenwid(10);cgpencol(0,1,0);
    cgdraw(X1,Y1,X2,Y1); cgdraw(X1,Y2,X2,Y2); cgdraw(X1,Y1,X1,Y2); cgdraw(X2,Y1,X2,Y2);
    cgpenwid(3);cgpencol(0,1,0)
    cgtext('ANSWER CORRECT !!',0,Y+25)
    cgtext(trial_type.FBwin{trial_type.inc},0,Y-25);
end
