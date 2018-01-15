function [ wp ] = proba_weight(p,gamma_PT )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

wp = (p.^gamma_PT)./((((p.^gamma_PT)+(1-p).^gamma_PT)).^(1/gamma_PT));
end

