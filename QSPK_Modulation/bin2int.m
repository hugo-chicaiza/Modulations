%% bin2int.m function definition
function [y] = bin2int(b)
% function to convert binary numbers to integers.
y=0;
for i=1:length(b)
    y=y+b(i)*2^(length(b)-i);
end;
% end of bin2int.m