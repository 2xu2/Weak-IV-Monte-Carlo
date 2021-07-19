function [AR,tF] = ARtFtest(stage2t, stage1F)
%AR and tF test for IV regress
if (stage2t > 1.96) & (stage1F > 104.3)
    AR = 1;
    tF = 1;
elseif (stage2t > 1.96) & (stage1F < 104.3)
    AR = 1;
    tF = 0;
else
    AR = 0;
    tF = 0;
end
end

