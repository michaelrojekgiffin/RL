% script just creates the conditions for my log form
num_sub = 100;

sub_log = NaN(num_sub, 2);

for ii = 1:num_sub
    sub_log(ii, 1) = ii;
    condition = mod(ii, 4);
    if condition == 0
        condition = 4;
    end
    sub_log(ii, 2) = condition;
end