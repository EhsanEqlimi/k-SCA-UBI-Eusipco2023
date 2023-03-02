function FnPlotFig1a
m_Set=[3,4,5]; % The number of sensors a.k.a the ambient dimension
n_Set=[4,5,6,8;5,6,7,8;6,7,8,9];% The number of sources (ith row correspond ith elemnet in m_set
% Prod=cartprod(m_Set,n_Set(1,:)');
% BAS=rand(4,12);

R1=load ('Result_Sigma_1e-6_SNR_Inf.mat');
R2=load ('Result_Sigma_1e-4_SNR_Inf.mat');
R3=load ('Result_Sigma_1e-3_SNR_Inf.mat');
R4=load ('Result_Sigma_1e-3_SNR_Inf_4_8.mat');
R5=load ('Result_Sigma_1e-3_SNR_Inf_m_5.mat');
R6=load('Result_81Trial_40_30.mat');

BAS=[R1.BAS_Set;R2.BAS_Set;R3.BAS_Set R4.BAS_Set R5.BAS_Set(1:2)];
NMSE=[R1.NMSE_Set;R2.NMSE_Set;R3.NMSE_Set R4.NMSE_Set R5.NMSE_Set(1:2)];
Norm2Error=[R1.Norm2Error_Set;R2.Norm2Error_Set;R3.Norm2Error_Set R4.Norm2Error_Set R5.Norm2Error_Set(1:2)];
Legend={'\sigma_{off} = 10^{-6}','\sigma_{off}=10^{-4}','\sigma_{off}=10^{-3}'};
% plot(BAS)
%% Plot Config.
% Color={'r','b',[1 .33 .99],[.49 1 .63],'n','k'};
Color={'r','b','k'};
Marker={'--rs','--*','<-','--x','--o','->'};
f = figure;
f.Name = 'First figure';
Marker2={'o','+','<','x','<','>'};

for j=1:3
    hold on
    Indici=log10(BAS(j,:));
    %     Indici=(Norm2Error(j,:));
    
    
    
    %'MarkerFaceColor',Color{j},'MarkerEdgeColor',Color{j},'MarkerSize',8,...
    %         PL=plot(1:4,Indici(1:4),5:8,Indici(5:8),9:12,Indici(9:12));
    PL=plot(1:4,Indici(1:4),5:8,Indici(5:8),9:10,Indici(9:10));
    %           PL=plot(1:4,Indici(1:4),5:7,Indici(5:7));
    %     PL=plot(1:4,Indici(1:4),5:7,Indici(5:7));
    
    
    
    PLCurrent(j)=PL(1);
    for t=1:3
        PL(t).Marker=Marker2{j};
        PL(t).MarkerSize=6;
        PL(t).MarkerEdgeColor=Color{j};
        PL(t).MarkerFaceColor=Color{j};
        PL(t).LineWidth=2;
        PL(t).Color=Color{j};
    end
    %         ylabel([Indici '-' CohInds{i}]);
    ylabel('Log_{10}(BAS)');
    %           ylabel('NMSE (dB)');
    %               ylabel('Norm2-Error');
    
    
    Alln=R1.n_Set;
    Ct=0;
    for w=1:3
        for t=1:4
            Ct=Ct+1;
            XTickLs{Ct}=['[' num2str(R1.m_Set(w)) ',' num2str(R1.n_Set(w,t)) ']'];
            
        end
    end
    xticklabels(XTickLs(1:10))
    set(gca, 'fontsize', 11,'fontweight','bold');
    
    %    set(gca,'YTick',[-pi 0 pi], 'YTickLabel', {'-\pi','0','\pi'}, 'fontsize', 18);
    %     set(gca, 'XTickLabel', , 'fontsize', 6);
    
    xlabel('[m , n]');
    
    legend(PLCurrent,Legend,'location','Best','fontsize', 11,'fontweight','bold')
    
end
