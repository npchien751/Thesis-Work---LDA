function [  ] = SFP_FinalFigures( class_results, fig )


% Read in SFP_classification output variables
headers_used = class_results.headers_used;
rawdata.endmembers = class_results.rawdata.endmembers;
lowGW_observed = class_results.lowGW_observed;
lowGW_synthetic = class_results.lowGW_synthetic;
highGW_observed = class_results.highGW_observed;
highGW_synthetic.chemistry = class_results.highGW_synthetic.chemistry;
highGW_synthetic.type = class_results.highGW_synthetic.type;
constants = class_results.constants;
coefficients = class_results.coefficients;
conf_matrix = class_results.conf_matrix;
conf_matrix_percents = class_results.conf_matrix_percents;
corr_synth_data_and_scores = class_results.corr_synth_data_and_scores;
bckwd_feat_sel = class_results.bckwd_feat_sel;
SWIFT_type.type_and_probabilities = class_results.SWIFT_type.type_and_probabilities;
SWIFT_type.chem = class_results.SWIFT_type.chem;
SWIFT_type.scores = class_results.SWIFT_type.scores;
SWIFT_type.names = class_results.SWIFT_type.names;
SWIFT_type.exp = class_results.SWIFT_type.exp;

lowGW = class_results.lowGW;
highGW = class_results.highGW;
GW_FW_high = class_results.GW_high1{1};
GW_RS_high = class_results.GW_high1{2};
GW_SEP_high = class_results.GW_high1{3};
scores = class_results.scores;
cont_gw = class_results.cont_gw;
endmembers = class_results.endmembers_for_figures;
FW = endmembers{1};
RS = endmembers{2};
SEP = endmembers{3};

% Known sources of contamination
fwGW = SWIFT_type.chem((SWIFT_type.exp()==1),:);
rsGW = SWIFT_type.chem((SWIFT_type.exp()==2),:);
sepGW = SWIFT_type.chem((SWIFT_type.exp()==3),:);



% Set solute variables and number
%I = 1; 
%Na = 1; K = 2; Mg = 3; Ca = 4; Cl = 5; Br = 6; Sr = 7; Ba = 8; SO4 = 9;
Na = 1; K = 2; Mg = 3; Ca = 4; Cl = 5; Br = 6; Sr = 7; Ba = 8;
n_sol = 8;

% THESE LINES DETERMINES FIGURE output
% Plotting commands
warning('off','MATLAB:Axes:NegativeDataInLogAxis')

% Set colors for categories  
high_c=[0 1 0];  
syn_c=[0 0 0];  % black
FW_c=[1 0 0];  % red
RS_c=[1 0.65 0];  % orange
SEP_c=[0.5 0.5 0.25];  % brown (dull)


color(1,:)=FW_c;
color(2,:)=RS_c;
color(3,:)=SEP_c;

panel_axislabeltext_size=20;
panel_numtext_size=14;
p=panel_axislabeltext_size;
p2=panel_numtext_size;

SWIFT_FW=SWIFT_type.chem((SWIFT_type.type_and_probabilities(:,1)==1),:);
SWIFT_RS=SWIFT_type.chem((SWIFT_type.type_and_probabilities(:,1)==2),:);
SWIFT_SEP=SWIFT_type.chem((SWIFT_type.type_and_probabilities(:,1)==3),:);
SWIFT_FW_s=SWIFT_type.scores((SWIFT_type.type_and_probabilities(:,1)==1),:);
SWIFT_RS_s=SWIFT_type.scores((SWIFT_type.type_and_probabilities(:,1)==2),:);
SWIFT_SEP_s=SWIFT_type.scores((SWIFT_type.type_and_probabilities(:,1)==3),:);
SWIFT_FW_p = SWIFT_type.type_and_probabilities((SWIFT_type.type_and_probabilities(:,1)==1),:);
SWIFT_RS_p = SWIFT_type.type_and_probabilities((SWIFT_type.type_and_probabilities(:,1)==2),:);
SWIFT_SEP_p = SWIFT_type.type_and_probabilities((SWIFT_type.type_and_probabilities(:,1)==3),:);
SWIFT_FW_e = SWIFT_type.exp((SWIFT_type.type_and_probabilities(:,1)==1),:);
SWIFT_RS_e = SWIFT_type.exp((SWIFT_type.type_and_probabilities(:,1)==2),:);
SWIFT_SEP_e = SWIFT_type.exp((SWIFT_type.type_and_probabilities(:,1)==3),:);
figno=0;



% Posterior probabilities (FIG 1)
% Organic Waste trial synthetic data
% Matlab structure: organic_class
if any(fig==1)
figno=figno+1;
figure('units','inches','position',[.1 .1 8 4]);
    data_9sol = [13 2; 19 2; 22 2; 24 3; 25 5; 29 10; 31 10];
    class_accur = bar(data_9sol, 'stacked');
    set(class_accur,{'FaceColor'},{[0 0 .8];[.95 0 0]});
    set(gca,'xticklabel', {'99%','95%','90%','80%','70%','50%','All'});
    set(gca,'fontsize',35)
    xlabel('Posterior Probability Ranges','FontSize',40);
    ylabel('Number of Samples Classified','FontSize',40);
    legend_1 = legend('Correct','Incorrect','location','Northwest');
    set(legend_1,'Fontsize',30);
end

% LDA scores for organic waste trial with symbolized misclasses (FIG 2)
if any(fig==2)
figno=figno+1;
figure('units','inches','position',[.1 .1 8 4]);
    % Synthetic data
    le1 = gscatter(scores(:,1),scores(:,2),highGW_synthetic.type,color,'oooooo',6,'off');
    hold on; 
    % Correctly classified samples
    le2 = plot(SWIFT_FW_s(1:7,1),SWIFT_FW_s(1:7,2),'o','MarkerSize',15,'MarkerFaceColor',FW_c,'MarkerEdgeColor','k');
    le3 = plot(SWIFT_RS_s(1:15,1),SWIFT_RS_s(1:15,2),'o','MarkerSize',15,'MarkerFaceColor',RS_c,'MarkerEdgeColor','k');
    le4 = plot(SWIFT_SEP_s(3:11,1),SWIFT_SEP_s(3:11,2),'o','MarkerSize',15,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','k');
    % Organic waste confused for formation water
    le5 = plot(SWIFT_FW_s(8:13,1),SWIFT_FW_s(8:13,2),'s','MarkerSize',18,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','b','LineWidth',1.5);
    % organic waste confused for road salt
    le6 = plot(SWIFT_RS_s(16:17,1),SWIFT_RS_s(16:17,2),'d','MarkerSize',18,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','b','LineWidth',1.5);
    % road salt confused for organic waste
    le7 = plot(SWIFT_SEP_s(1:2,1),SWIFT_SEP_s(1:2,2),'h','MarkerSize',18,'MarkerFaceColor',RS_c,'MarkerEdgeColor','b','LineWidth',1.5);
    % Axis Options
    plot([0 0],[-20 20],'--k');
    plot([-20 20],[0 0],'--k');
    axis([-20 20 -20 20]);
    set(gca,'FontSize',35);
    xlabel('Score1','FontSize',p,'FontSize',40);
    ylabel('Score2','FontSize',p,'FontSize',40);
    fig2_lgnd = legend([le5 le6 le7],{'Organic Waste misclassified as Formation Water','Organic Waste misclassified as Road Salt','Road Salt misclassified as Organic Waste'});
    fig2_lgnd_title = get(fig2_lgnd,'Title');
    set(fig2_lgnd_title,'string','Symbology for Misclassified Validation Data');
    set(fig2_lgnd,'Fontsize',20);
    set(fig2_lgnd,'Location','Southeast');
end

% Bivariate plot of LDA coefficients (FIG 3)
if any(fig==3)
figno=figno+1;
figure('units','inches','position',[.1 .1 8 4]);
    % Synthetic data
    gscatter(scores(:,1),scores(:,2),highGW_synthetic.type,color,'oooooo',6,'off');
    hold on;
    % Set score-coefficient vectors
    coeff_labels = {'Na','K','Mg','Cl','Br','Sr','Ba','SO4'};
    coeff_s1 = [-0.89,-1.07,-0.04,-1.07,6.39,1.39,-0.14,0.09];
    coeff_s2 = [-0.3,-1.28,-0.85,3.29,-1.76,0.44,1.36,0.07];
    quiver(zeros(1,length(coeff_s1)),zeros(1,length(coeff_s2)),coeff_s1,coeff_s2,'LineWidth',1.5,'MarkerFaceColor','b','MarkerEdgeColor','b')
    hold on;
    % Add solute labels
    for i=1:length(coeff_labels)
        a = text(coeff_s1(1,i),coeff_s2(1,i),coeff_labels(1,i));
        set(a,'FontSize',20);
    end
    % Axis Options - zoomed in
    plot([0 0],[-4 4],'--k');
    plot([-2 8],[0 0],'--k');
    axis([-2 8 -4 4]);
%     % Axis Options - full scope
%     plot([0 0],[-20 20],'--k');
%     plot([-20 20],[0 0],'--k');
%     axis([-20 20 -20 20]);
    set(gca,'FontSize',30);
    xlabel('Score1','FontSize',p,'FontSize',30);
    ylabel('Score2','FontSize',p,'FontSize',30);
end

% Chloride-Bromide against scores panels(FIG 4)
if any(fig==4)
figno=figno+1;
figure('units','inches','position',[.1 .1 8 4]);
subplot(2,2,2);
    % Synthetic data
    gscatter(scores(:,1),scores(:,2),highGW_synthetic.type,color,'oooooo',6,'off');
    hold on;
    plot(SWIFT_FW_s(:,1),SWIFT_FW_s(:,2),'o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k');
    plot(SWIFT_RS_s(:,1),SWIFT_RS_s(:,2),'o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k');
    plot(SWIFT_SEP_s(:,1),SWIFT_SEP_s(:,2),'o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k');
    % Axis options
    plot([0 0],[-20 20],'--k');
    plot([-20 20],[0 0],'--k');
    axis([-20 20 -20 20]);
    set(gca,'FontSize',p2);
    xlabel('Score1','FontSize',p);
    ylabel('Score2','FontSize',p);
subplot(2,2,4);
    % Synthetic data
    gscatter(scores(:,1),scores(:,2),highGW_synthetic.type,color,'oooooo',6,'off');
    hold on; 
    % Correctly classified samples
    plot(SWIFT_FW_s(1:7,1),SWIFT_FW_s(1:7,2),'o','MarkerSize',10,'MarkerFaceColor',FW_c,'MarkerEdgeColor','k');
    plot(SWIFT_RS_s(1:15,1),SWIFT_RS_s(1:15,2),'o','MarkerSize',10,'MarkerFaceColor',RS_c,'MarkerEdgeColor','k');
    plot(SWIFT_SEP_s(3:11,1),SWIFT_SEP_s(3:11,2),'o','MarkerSize',10,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','k');
    % Organic waste confused for formation water
    plot(SWIFT_FW_s(8:13,1),SWIFT_FW_s(8:13,2),'o','MarkerSize',10,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','k');
    % organic waste confused for road salt
    plot(SWIFT_RS_s(16:17,1),SWIFT_RS_s(16:17,2),'o','MarkerSize',10,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','k');
    % road salt confused for organic waste
    plot(SWIFT_SEP_s(1:2,1),SWIFT_SEP_s(1:2,2),'o','MarkerSize',10,'MarkerFaceColor',RS_c,'MarkerEdgeColor','k');
    % Axis Options
    plot([0 0],[-20 20],'--k');
    plot([-20 20],[0 0],'--k');
    axis([-20 20 -20 20]);
    set(gca,'FontSize',p2);
    xlabel('Score1','FontSize',p);
    ylabel('Score2','FontSize',p);
% Calculate Chloride:Bromide ratios for different end members
lowGW_ClBr = lowGW(:,5)./lowGW(:,6);
highGW_ClBr = highGW(:,5)./highGW(:,6);
FW_ClBr = FW(:,5)./FW(:,6);
RS_ClBr = RS(:,5)./RS(:,6);
SEP_ClBr = SEP(:,5)./SEP(:,6);
fwGW_ClBr = fwGW(:,5)./fwGW(:,6);
rsGW_ClBr = rsGW(:,5)./rsGW(:,6);
sepGW_ClBr = sepGW(:,5)./sepGW(:,6);
% Cl:Br plots with unclassified groundwater
subplot(2,2,1);
    loglog(lowGW(:,Cl),lowGW_ClBr,'o','MarkerSize',6,'MarkerFaceColor',syn_c,'MarkerEdgeColor','k');
    hold on;
    loglog(highGW(:,Cl),highGW_ClBr,'o','MarkerSize',6,'MarkerEdgeColor','k');  
    loglog(FW(:,Cl),FW_ClBr,'o','MarkerSize',6,'MarkerFaceColor',FW_c,'MarkerEdgeColor','k');
    loglog(RS(:,Cl),RS_ClBr,'o','MarkerSize',6,'MarkerFaceColor',RS_c,'MarkerEdgeColor','k');
    loglog(SEP(:,Cl),SEP_ClBr,'o','MarkerSize',6,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','k');
    hold on;
    xlabel('Chloride (ppm)','FontSize',15);
    ylabel('Chloride:Bromide Ratio','FontSize',15);
    axis([.1 1000000 10 100000]);
    legend('Observed Low Salinity GW','Observed High Salinity GW','Formation Water','Road Salt','Organic Waste','Location','northwest');
% Cl:Br plots with known sources of contamination
subplot(2,2,3);
    loglog(lowGW(:,Cl),lowGW_ClBr,'o','MarkerSize',6,'MarkerFaceColor',syn_c,'MarkerEdgeColor','k');
    hold on;
    loglog(fwGW(:,Cl),fwGW_ClBr,'*','MarkerSize',6,'MarkerEdgeColor',FW_c);
    loglog(rsGW(:,Cl),rsGW_ClBr,'*','MarkerSize',6,'MarkerEdgeColor',RS_c);
    loglog(sepGW(:,Cl),sepGW_ClBr,'*','MarkerSize',6,'MarkerEdgeColor',SEP_c);
    loglog(FW(:,Cl),FW_ClBr,'o','MarkerSize',6,'MarkerFaceColor',FW_c,'MarkerEdgeColor','k');
    loglog(RS(:,Cl),RS_ClBr,'o','MarkerSize',6,'MarkerFaceColor',RS_c,'MarkerEdgeColor','k');
    loglog(SEP(:,Cl),SEP_ClBr,'o','MarkerSize',6,'MarkerFaceColor',SEP_c,'MarkerEdgeColor','k');
    hold on; 
    xlabel('Chloride (ppm)','FontSize',15);
    ylabel('Chloride:Bromide Ratio','FontSize',15);
    axis([.1 1000000 10 100000]);
end

% Bivariate plot of LDA coefficients with multiplier (FIG 5)
if any(fig==5)
figno=figno+1;
figure('units','inches','position',[.1 .1 8 4]);
    % Synthetic data
    gscatter(scores(:,1),scores(:,2),highGW_synthetic.type,color,'oooooo',6,'off');
    hold on;
    % Set score-coefficient vectors
    coeff_labels = {'Na','K','Mg','Cl','Br','Sr','Ba','SO4'};
    coeff_s1 = 2.5*[-0.89,-1.07,-0.04,-1.07,6.39,1.39,-0.14,0.09];
    coeff_s2 = 2.5*[-0.3,-1.28,-0.85,3.29,-1.76,0.44,1.36,0.07];
    quiver(zeros(1,length(coeff_s1)),zeros(1,length(coeff_s2)),coeff_s1,coeff_s2,'LineWidth',2.6,'MarkerFaceColor','b','MarkerEdgeColor','b')
    hold on;
    % Add solute labels
    for i=1:length(coeff_labels)
        a = text(coeff_s1(1,i),coeff_s2(1,i),coeff_labels(1,i));
        set(a,'FontSize',20);
    end

    % Axis Options - full scope
    plot([0 0],[-20 20],'--k');
    plot([-20 20],[0 0],'--k');
    axis([-20 20 -20 20]);
    set(gca,'FontSize',30);
    xlabel('Score1','FontSize',p,'FontSize',30);
    ylabel('Score2','FontSize',p,'FontSize',30);
end

%
%
%
%%% Defense Presentation Extra Plots
%
%
%
% CDF plots of observed and synthetic low salinity GW data sets (FIG 6)
if any(fig==6);
figno=figno+1;
figure;
subplot(2,4,1);
    [f,x]=ecdf(lowGW(:,Cl));   
    semilogx(x,f,'Color','r','LineWidth',2);
      set(gca,'FontSize',p2);
      title('Empirical CDFs of Low Salinity Data','FontSize',20);
      xlabel('Chloride (ppm)','FontSize',p);
      ylabel('F(x)','FontSize',p);
    hold on;
    [f,x]=ecdf(lowGW_synthetic(:,Cl));
    semilogx(x,f,'Color','b','LineWidth',2);
    legend('Observed','Synthetic');
    axis([0.1 100 0 1]);
subplot(2,4,2);
    [f,x]=ecdf(lowGW(:,Br));   
    semilogx(x,f,'Color','r','LineWidth',2);
      set(gca,'FontSize',p2);
      xlabel('Bromide (ppm)','FontSize',p);
      ylabel('F(x)','FontSize',p);
    hold on;
    [f,x]=ecdf(lowGW_synthetic(:,Br));
    semilogx(x,f,'Color','b','LineWidth',2);
    axis([0.001 1 0 1]);
subplot(2,4,3);
    [f,x]=ecdf(lowGW(:,Na));   
    semilogx(x,f,'Color','r','LineWidth',2);
      set(gca,'FontSize',p2);
      xlabel('Sodium (ppm)','FontSize',p);
      ylabel('F(x)','FontSize',p);
    hold on;
    [f,x]=ecdf(lowGW_synthetic(:,Na));
    semilogx(x,f,'Color','b','LineWidth',2);
    axis([0.1 1000 0 1]);
% subplot(2,4,4);
%     [f,x]=ecdf(lowGW(:,SO4));   
%     semilogx(x,f,'Color','r','LineWidth',2);
%       set(gca,'FontSize',p2);
%       xlabel('Iodide (ppb)','FontSize',p);
%       ylabel('F(x)','FontSize',p);
%     hold on;
%     [f,x]=ecdf(lowGW_synthetic(:,SO4));
%     semilogx(x,f,'Color','b','LineWidth',2);
%     axis([1 100 0 1]);
subplot(2,4,5);
    [f,x]=ecdf(lowGW(:,Ca));   
    semilogx(x,f,'Color','r','LineWidth',2);
      set(gca,'FontSize',p2);
      xlabel('Calcium (ppm)','FontSize',p);
      ylabel('F(x)','FontSize',p);
    hold on;
    [f,x]=ecdf(lowGW_synthetic(:,Ca));
    semilogx(x,f,'Color','b','LineWidth',2);
    axis([1 100 0 1]);
subplot(2,4,6);
    [f,x]=ecdf(lowGW(:,Mg));   
    semilogx(x,f,'Color','r','LineWidth',2);
      set(gca,'FontSize',p2);
      xlabel('Magnesium (ppm)','FontSize',p);
      ylabel('F(x)','FontSize',p);
    hold on;
    [f,x]=ecdf(lowGW_synthetic(:,Mg));
    semilogx(x,f,'Color','b','LineWidth',2);
    axis([0.1 100 0 1]);
subplot(2,4,7);
    [f,x]=ecdf(lowGW(:,K));   
    semilogx(x,f,'Color','r','LineWidth',2);
      set(gca,'FontSize',p2);
      xlabel('Potassium (ppm)','FontSize',p);
      ylabel('F(x)','FontSize',p);
    hold on;
    [f,x]=ecdf(lowGW_synthetic(:,K));
    semilogx(x,f,'Color','b','LineWidth',2);
    axis([0.1 10 0 1]);
subplot(2,4,8);
    [f,x]=ecdf(lowGW(:,Sr));   
    semilogx(x,f,'Color','r','LineWidth',2);
      set(gca,'FontSize',p2);
      xlabel('Strontium (ppb)','FontSize',p);
      ylabel('F(x)','FontSize',p);
    hold on;
    [f,x]=ecdf(lowGW_synthetic(:,Sr));
    semilogx(x,f,'Color','b','LineWidth',2);
    axis([10 10000 0 1]);
end
    
% CDF plots of Cl concentrations in high salinity GW data sets (FIG 7)
if any(fig==7);
figno=figno+1;
figure;
[f,x]=ecdf(highGW(:,Cl));   
semilogx(x,f,'Color',high_c,'LineWidth',2);
  set(gca,'FontSize',p2);
  title('Empirical CDFs of Chloride Data','FontSize',20);
  xlabel('Chloride (ppm)','FontSize',p);
  ylabel('F(x)','FontSize',p);
hold on;
[f,x]=ecdf(GW_FW_high(:,Cl));
semilogx(x,f,'Color',FW_c,'LineWidth',2);
[f,x]=ecdf(GW_RS_high(:,Cl));
semilogx(x,f,'Color',RS_c,'LineWidth',2);
[f,x]=ecdf(GW_SEP_high(:,Cl));
semilogx(x,f,'Color',SEP_c,'LineWidth',2);
  legend('Observed High Salinity GW','Formation Water','Road Salt','Organic Waste','Location','best');
end

% Scatter plot matrix of observed low salinity groundwater chemistry; log
% transformed  (FIG 9)
if any(fig==9);
figno=figno+1;
mini(1,:)=min(lowGW);
maxi(1,:)=max(lowGW);
mini(2,:)=min(SWIFT_type.chem);
maxi(2,:)=max(SWIFT_type.chem);
mini(3,:)=min(highGW);
maxi(3,:)=max(highGW);
mini(3,:)=min(lowGW_synthetic);
maxi(3,:)=max(lowGW_synthetic);
mini(4,:)=min(cont_gw);
maxi(4,:)=max(cont_gw);
mini=min(mini);
maxi=max(maxi);
figure('units','inches','position',[.1 .1 6 9]);
subplot(3,2,1);
    gscatter(cont_gw(:,Cl),cont_gw(:,7),highGW_synthetic.type,color,'oooooo',5,'off');
    set(gca,'FontSize',p2,'xscale','log','yscale','log');
    xlabel('Chloride (ppm)','FontSize',p);
    ylabel('Bromide (ppm)','FontSize',p);
    axis([mini(:,Cl) maxi(:,Cl) mini(:,7) maxi(:,7)]);
subplot(3,2,2);
    gscatter(cont_gw(:,5),cont_gw(:,4),highGW_synthetic.type,color,'oooooo',5,'off');
    set(gca,'FontSize',p2,'xscale','log','yscale','log');
    xlabel('Calcium (ppm)','FontSize',p);
    ylabel('Magnesium (ppm)','FontSize',p);
    legend('Formation Water','Road Salt','Septic Effluent','Animal Waste');
    axis([mini(:,5) maxi(:,5) mini(:,4) maxi(:,4)]);
subplot(3,2,3);
    gscatter(cont_gw(:,Cl),cont_gw(:,2),highGW_synthetic.type,color,'oooooo',5,'off');
    set(gca,'FontSize',p2,'xscale','log','yscale','log');
    xlabel('Chloride (ppm)','FontSize',p);
    ylabel('Sodium (ppm)','FontSize',p);
    axis([mini(:,Cl) maxi(:,Cl) mini(:,2) maxi(:,2)]);
subplot(3,2,4);
    gscatter(cont_gw(:,5),cont_gw(:,3),highGW_synthetic.type,color,'oooooo',5,'off');
    set(gca,'FontSize',p2,'xscale','log','yscale','log');
    xlabel('Calcium (ppm)','FontSize',p);
    ylabel('Potassium (ppm)','FontSize',p);
    axis([mini(:,5) maxi(:,5) mini(:,3) maxi(:,3)]);
subplot(3,2,5);
    gscatter(cont_gw(:,Cl),cont_gw(:,1),highGW_synthetic.type,color,'oooooo',5,'off');
    set(gca,'FontSize',p2,'xscale','log','yscale','log');
    xlabel('Chloride (ppm)','FontSize',p);
    ylabel('Iodine (ppb)','FontSize',p);
    axis([mini(:,Cl) maxi(:,Cl) mini(:,1) maxi(:,1)]);
subplot(3,2,6);
    gscatter(cont_gw(:,5),cont_gw(:,8),highGW_synthetic.type,color,'oooooo',5,'off');
    set(gca,'FontSize',p2,'xscale','log','yscale','log');
    xlabel('Calcium (ppm)','FontSize',p);
    ylabel('Strontium (ppb)','FontSize',p);
    axis([mini(:,5) maxi(:,5) mini(:,8) maxi(:,8)]);
end

end 
