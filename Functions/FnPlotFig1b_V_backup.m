function FnPlotFig1b

load  Ransa_Movahadi_100trial.mat
for i=1:size(MixingVectorerror_Movahedi,1)
    for j=1:1:size(MixingVectorerror_Movahedi,2)
        nhat_mov(i,j)=length(find(MixingVectorerror_Movahedi{i,j}<.1));
        nhat_ours(i,j)=length(find(MixingVectorerror_Ours{i,j}<.1));

    
    end
end
% 
% b = bar(1:10,[mean(nhat_mov,1) mean(nhat_ours,1)]);
% b.FaceColor = 'g';
% figure,stem(1:5,mean(nhat_mov,1),'--k>')
% hold on
% stem(1:5,mean(nhat_ours,1),'-*rs')
% bpcombined = [mean(nhat_mov,1), mean(nhat_ours,1)];
% hb = bar(1:10, bpcombined, 'grouped');
close all
XTick={'[3,5]', '[3,6]','[3,8]','[3,9]','[3,10]'};

% title('Blocking Probability vs Routing Level');
Legend={'\sigma_{off} = 10^{-6}','\sigma_{off}=10^{-4}','\sigma_{off}=10^{-3}'};
% plot(BAS)
%% Plot Config.

figure,
plot(1:5,log10(mean(BAS_Detected_Ours,1)),'-.rs','MarkerSize',8,'LineWidth',2);
LogMeanBAS_Ours=log10(mean(BAS_Detected_Ours,1));
LogMeanBAS_Mov=log10(mean(BAS_Detected_Movahedi,1));

xlim([.85 5.32])

text([1 1.6 2.7 3.7 4.7],[-2.9,-2.78, -2.4, -2.24, -2.10 ] , '$$\hat{n}$$','Interpreter','Latex','fontsize', 17,'fontweight','bold' )
% text([1.1 1.8 3 4 4.9],LogMeanBAS_Ours+.1 , '$$\hat{n}$$','Interpreter','Latex','fontsize', 17,'fontweight','bold' )

hold on
plot(1:5,log10(mean(BAS_Detected_Movahedi,1)),'--k>','MarkerSize',8,'LineWidth',2);
text([1.1 1.8 3 4 4.9],LogMeanBAS_Mov+.1 , '$$\hat{n}$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' )

set(gca, 'XTick', 1:5,...
    'XTickLabel', XTick);
% xticklabels(XTick,'Interpreter','Latex' )
% set(gca, 'fontsize', 10,'fontweight','bold','MarkerSize',15);

%    set(gca,'YTick',[-pi 0 pi], 'YTickLabel', {'-\pi','0','\pi'}, 'fontsize', 18);
%     set(gca, 'XTickLabel', , 'fontsize', 6);

xlabel('$$[m , n]$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' );
ylabel('$$\log_{10}(\textrm{BAS})$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' )
legend({'Proposed Algorithm','Algorithim in [16]'},'location','Best')
set(gca, 'fontsize', 10,'fontweight','bold');
box on
