X1 = 1.5e-3:1e-4:2.5e-3;
Y1 = [];
Y2 = [];
for param = 0.0015:0.0001:0.0025
    result = sprintf('Guided_RGB_RID_%.4f.mat',param);
    eval(['load ' num2str(result)]);
    Y1 = [Y1 mPSNR(end)];
    Y2 = [Y2 mSSIM(end)];
end
createfigure(X1, Y1, Y2);