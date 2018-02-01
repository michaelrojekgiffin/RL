all_conds = NaN(100, 2);
cond_count = 0;
for ii = 1:length(all_conds)
    cond_count = cond_count +1;
    all_conds(ii, 1) = ii;
    all_conds(ii, 2) = cond_count;
    if cond_count == 4
        cond_count = 0;
    end
end
