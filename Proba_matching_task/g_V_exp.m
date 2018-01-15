function [gx] = g_V_exp(x,Phi,u,in)

% variables
G1 = u(1,:);
L1 = u(2,:);
P1 = u(3,:);
G2 = u(4,:);
L2 = u(5,:);
P2 = u(6,:);

% params
pphi = Phi(1);
pr = Phi(2);
plam = Phi(3);
pgam = Phi(4);

% get values
O1 = get_v(G1,L1,P1,pr,plam,pgam);
O2 = get_v(G2,L2,P2,pr,plam,pgam);

% get choices
gx = 1./(1+exp(-(O1-O2)*pphi));

prec=10^-9;
gx(gx<prec)=prec; gx(gx>(1-prec))=1-prec;

end