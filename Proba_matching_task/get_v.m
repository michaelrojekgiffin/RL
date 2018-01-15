    function [v] = get_v(G,L,P,fr,flam,fgam)
        % compute Expected value according to Prospect Theory
        v = proba_weight(P,fgam)*(G.^fr) - flam.*proba_weight(1-P,fgam)*(L.^fr);    
    end