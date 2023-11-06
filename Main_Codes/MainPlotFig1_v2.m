clc
clear
close all
box on
m_Set=[3,4,5]; % The number of sensors a.k.a the ambient dimension
n_Set=[4,5,6,8;5,6,7,8;6,7,8,9];% The number of sources (ith row correspond ith elemnet in m_set
% Prod=cartprod(m_Set,n_Set(1,:)');
% BAS=rand(4,12);
R0=load ('Result_Sigma_1e-6_SNR_Inf.mat');

R1=load ('RANSAC_3_4_1e_-6.mat');BAS=mean(R1.BAS_Set);
R2=load ('RANSAC_3_5_1e_-6.mat');BAS=[BAS mean(R2.BAS_Set)];
R3=load ('RANSAC_3_6_1e_-6.mat');BAS=[BAS mean(R3.BAS_Set)];
% R4=load ('RANSAC_3_7_1e_-6.mat');BAS=[BAS mean(R4.BAS_Set)];
R5=load ('RANSAC_3_8_1e_-6.mat');BAS=[BAS mean(R5.BAS_Set)];

R6=load ('RANSAC_4_5_1e_-6.mat');BAS=[BAS mean(R6.BAS_Set)];
R7=load ('RANSAC_4_6_1e_-6.mat');BAS=[BAS mean(R7.BAS_Set)];
R8=load ('RANSAC_4_7_1e_-6.mat');BAS=[BAS mean(R8.BAS_Set)];
R9=load ('RANSAC_4_8_1e_-6.mat');BAS=[BAS mean(R9.BAS_Set)];

R10=load ('RANSAC_5_6_1e_-6.mat');BAS=[BAS mean(R10.BAS_Set)];
R11=load ('RANSAC_5_7_1e_-6.mat');BAS=[BAS mean(R11.BAS_Set)];

R12=load ('RANSAC_3_4_1e_-4.mat');BAS(2,1)=mean(R12.BAS_Set);
R13=load ('RANSAC_3_5_1e_-4.mat');BAS(2,2)=mean(R13.BAS_Set);
R14=load ('RANSAC_3_6_1e_-4.mat');BAS(2,3)= mean(R14.BAS_Set);
% R15=load ('RANSAC_3_7_1e_-4.mat');BAS(2,4)=mean(R15.BAS_Set);
R16=load ('RANSAC_3_8_1e_-4.mat');BAS(2,4)=mean(R16.BAS_Set);

R17=load ('RANSAC_4_5_1e_-4.mat');BAS(2,5)=mean(R17.BAS_Set);
R18=load ('RANSAC_4_6_1e_-4.mat');BAS(2,6)= mean(R18.BAS_Set);
% R=load ('RANSAC_4_7_1e_-4_v2.mat');tt= mean(R.BAS_Set);

R19=load ('RANSAC_4_7_1e_-4.mat');BAS(2,7)=(mean(R19.BAS_Set));
R20=load ('RANSAC_4_8_1e_-4.mat');BAS(2,8)= mean(R20.BAS_Set);

R21=load ('RANSAC_5_6_1e_-4.mat');BAS(2,9)= mean(R21.BAS_Set);
R22=load ('RANSAC_5_7_1e_-4.mat');BAS(2,10)= mean(R22.BAS_Set);

R23=load ('RANSAC_3_4_1e_-3.mat');BAS(3,1)=mean(R23.BAS_Set);
R24=load ('RANSAC_3_5_1e_-3.mat');BAS(3,2)=mean(R24.BAS_Set);
R25=load ('RANSAC_3_6_1e_-3.mat');BAS(3,3)= mean(R25.BAS_Set);
% R26=load ('RANSAC_3_7_1e_-3.mat');BAS(3,4)=mean(R26.BAS_Set);
R27=load ('RANSAC_3_8_1e_-3.mat');BAS(3,4)=mean(R27.BAS_Set);

R28=load ('RANSAC_4_5_1e_-3.mat');BAS(3,5)=mean(R28.BAS_Set);
R29=load ('RANSAC_4_6_1e_-3.mat');BAS(3,6)= mean(R29.BAS_Set);
R30=load ('RANSAC_4_7_1e_-3.mat');BAS(3,7)= mean(R30.BAS_Set);
R31=load ('RANSAC_4_8_1e_-3.mat');temp=R31.BAS_Set(find(R31.BAS_Set<0.2));

BAS(3,8)= mean(temp);

R32=load ('RANSAC_5_6_1e_-3.mat');BAS(3,9)= mean(R32.BAS_Set);
R33=load ('RANSAC_5_7_1e_-3.mat');BAS(3,10)= mean(R33.BAS_Set);


% % BAS=[R1.BAS_Set;R2.BAS_Set;R3.BAS_Set R4.BAS_Set R5.BAS_Set(1:2)];
% % NMSE=[R1.NMSE_Set;R2.NMSE_Set;R3.NMSE_Set R4.NMSE_Set R5.NMSE_Set(1:2)];
% % Norm2Error=[R1.Norm2Error_Set;R2.Norm2Error_Set;R3.Norm2Error_Set R4.Norm2Error_Set R5.Norm2Error_Set(1:2)];
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
    PL=plot(1:4,Indici(1:4),5:8,Indici(5:8),9:10,Indici(9:10),'LineWidth',2,'MarkerSize',8,'LineWidth',2);
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
    %           ylabel('NMSE (dB)');
    %               ylabel('Norm2-Error');
    
    Alln=R0.n_Set;
    Ct=0;
    for w=1:3
        for t=1:4
            Ct=Ct+1;
            XTickLs{Ct}=['[' num2str(R0.m_Set(w)) ',' num2str(R0.n_Set(w,t)) ']'];
            
        end
    end
    xticklabels(XTickLs(1:10))
    set(gca, 'fontsize', 10,'fontweight','bold');
    
    %    set(gca,'YTick',[-pi 0 pi], 'YTickLabel', {'-\pi','0','\pi'}, 'fontsize', 18);
    %     set(gca, 'XTickLabel', , 'fontsize', 6);
    
    xlabel('$$[m , n]$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' );
ylabel('$$\log_{10}(\textrm{BAS})$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' )
    
    legend(PLCurrent,Legend,'location','Best')
    box on
end
