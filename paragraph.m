
para = 1:1:100;
CZ = stopCZ(para, 1);
C = stopC(para, 1);
Z = stopZ(para, 1);
% set(gca, 'YScale', 'log')
plot(para,CZ,'--s','LineWidth',1.5,...
    'MarkerFaceColor','k',...
    'MarkerSize',3);
hold all;
plot(para,C,'-.s','LineWidth',1.5,...
    'MarkerFaceColor','k',...
    'MarkerSize',3);
plot(para,Z,'-s','LineWidth',1.5,...
    'MarkerFaceColor','k',...
    'MarkerSize',3);
axis([1 100 0 2e0]);
legend('|C_{k+1}-Z_{k+1}|','|C_{k+1}-C_{k}|','|Z_{k+1}-Z_{k}|');
xlabel({'Iteration k'});
ylabel('Maximal Error');
grid on;
% set(gca, 'YScale', 'log')


