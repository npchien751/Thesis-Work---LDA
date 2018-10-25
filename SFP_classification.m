function [ output ] = SFP_classification( param_input )

% This tool creates a linear discriminant analysis model to classify
% the most probable sources of salinity in shallow groundwater using the
% methods outlined in Lautz et al (2014) and Chien & Lautz (2018). The 
% code is used to create synthetic chemical data for low salinity 
% groundwater and potential salinity sources (e.g. formation water, 
% road salt runoff). The synthetic end-members (low salinity groundwater
% and salinity sources) are then mixed using a simple two-component mixing
% model to generate synthetic datasets of high salinity groundwater created
% by mixing with the various sources of salinity.  The synthetic high
% salinity datasets are then used to train a discriminant analysis model to
% discriminate between the various populations (e.g. water mixed with
% formation water, road salt, animal waste, or septic effluent).  The
% discriminant analysis model is then used to determine the most probable
% sources of salinity in observations of high salinity shallow groundwater.
% More details can be found in Lautz et al (2014) and Chien & Lautz (2018).
%
% INPUT VARIABLES:
% param_input: Use of this model is dependent on the script, 
%               SFP_parameters.m which reads an properly formatted excel 
%               file and outputs the appropriate variables for use by this 
%               model.
%
% OUTPUT VARIABLES:
% output.headers        header rows for solute concentration output arrays
% output.rawdata.ALL    all raw data entered in the model
% output.rawdata.GW     all shallow groundwater data entered in the model
% output.rawdata.FW     all formation water data entered in the model
% output.rawdata.RS     all road salt runoff data entered in the model
% output.rawdata.SEP    all septic effluent data entered in the model
% output.rawdata.ANIM   all animal waste data entered in the model
% output.lowGW_observed     all observations of low salinity groundwater
% output.lowGW_synthetic    all synthetic low salinity groundwater samples
%                               generated by the model
% output.highGW_observed    all observations of high salinity groundwater
% output.highGW_synthetic.chemistry     combined training data set; includes
%                                       all synthetic samples of high 
%                                       salinity groundwater
% output.highGW_synthetic.type  array indicating classification of
%                                       training data
% output.constants          Constants for linear classifiers
% output.coefficients       Coefficients for linear classifiers
% output.conf_matrix        Confusion matrix showing classification error
%                               for training dataset
% output.conf_matrix_percents       Confusion matrix as percentages
% output.corr_synth_data_and_scores     Correlation coefficients for the
%                                           training data solute 
%                                           concentrations and the 
%                                           discriminant analysis scores
% output.bckwd_feat_sel      Results of backward feature selection;
%                            indicates which solutes were selected for final model
% output.SWIFT_type.type_and_probability    Classification and associated
%                                               probabilities for classified 
%                                               unknowns
% output.SWIFT_type.chem    Solute concentrations in classified samples
% output.SWIFT_type.scores  Discriminant analysis scores for classified
%                               samples
% output.SWIFT_type.name    Sample ID numbers for classified high sal. GW
% output.SWIFT_type.exp     Integer value representing which endmember a
%                               specific sample belongs to


% Resets random number generator seed so the program always returns the 
% same result when rerun.
rng('default')
warning('off','all');

% Read in datafile salinity and solute parameters
test_size = param_input.test_size;     % Sets # of samples for synthetic datasets
Cl_cutoff = param_input.Cl_cutoff;     % Sets cutoff between high and low salinity datasets, in mg/L
Cl_mdl = param_input.Cl_mdl;           % Sets detection limit for Cl data, in mg/L
n_sol = param_input.n_sol;             % Sets # of solutes in the input data
Cl = param_input.Cl;                   % Location of Cl in solute dataset
prior_prob = param_input.prior;       % Set prior probabilities
% Create arrays for GW (observed drinking wells), and endmembers.
headers = param_input.headers;
endmembers = param_input.endmembers;
GW = param_input.GW;
GW_name = param_input.GW_name;
GW_exp = param_input.GW_exp;
maxI = param_input.maxI;               % Max. iterations for cov matrix
mixcap = param_input.mixcap;

% CREATE SYNTHETIC DATASET OF LOW SALINITY GROUNDWATER
% Separate GW observations above and below the set Cl cutoff value
lowGW = GW(GW(:,Cl)<Cl_cutoff,1:n_sol);
highGW = GW(GW(:,Cl)>Cl_cutoff,1:n_sol);
highGW_name = GW_name(GW(:,Cl)>Cl_cutoff,1);
if size(GW_exp, 1) == size(GW, 1)
    highGW_exp = GW_exp(GW(:,Cl)>Cl_cutoff,1);
else
    highGW_exp = NaN(size(GW,1),1);
end

% Compute mean and covariance matrix of the log-transformed low salinity GW
% data for output.  
lowGW_log = log(lowGW);
if any(isnan(lowGW_log(:)))
    [lowGW_log_mean, lowGW_log_C] = ecmnmle(lowGW_log,'diagonal',100);
    lowGW_log_mean = lowGW_log_mean';
else
    lowGW_log_mean = mean(lowGW_log);
    lowGW_log_C = cov(lowGW_log);
end

% Create two random sets of synthetic low salinity GW samples (n = test size).
% Use a random covariance matrix pulled from the Inverse Wishart distribution.
[dim1, dim2] = size(lowGW_log_C);
C_wish_GW = NaN(dim1, dim2, length(test_size));
lowGW_rand_log1 = NaN(length(test_size), length(lowGW_log_mean));
lowGW_rand_log2 = NaN(length(test_size), length(lowGW_log_mean));
for n = 1:test_size
    C_wish_GW(:,:,n) = iwishrnd(lowGW_log_C,size(lowGW_log,1))*(size(lowGW_log,1)-size(lowGW_log_C,1)-1);
    lowGW_rand_log1(n,:) = mvnrnd(lowGW_log_mean,C_wish_GW(:,:,n),1);
    lowGW_rand_log2(n,:) = mvnrnd(lowGW_log_mean,C_wish_GW(:,:,n),1);
end
lowGW_rand_untrans1 = exp(lowGW_rand_log1);
lowGW_rand_untrans2 = exp(lowGW_rand_log2);

% Use this loop to remove synthetic samples with Cl concentrations above
% the Cl cutoff threshold and below the Cl detection limit. Replace with 
% new random samples.  Repeat until you have the test size of samples, 
% with none above or below the cutoffs.
for n = 1:test_size;
    anySaline = any(lowGW_rand_untrans1(n,Cl)>Cl_cutoff)+any(lowGW_rand_untrans1(n,Cl)<Cl_mdl); 
    while anySaline; 
        lowGW_rand_log1(n,:) = mvnrnd(lowGW_log_mean,C_wish_GW(:,:,n),1);
        lowGW_rand_untrans1(n,:) = exp(lowGW_rand_log1(n,:));
        anySaline = any(lowGW_rand_untrans1(n,Cl)>Cl_cutoff) + any(lowGW_rand_untrans1(n,Cl)<Cl_mdl); 
    end;
end
for n = 1:test_size;
    anySaline = any(lowGW_rand_untrans2(n,Cl)>Cl_cutoff)+any(lowGW_rand_untrans2(n,Cl)<Cl_mdl); 
    while anySaline; 
        lowGW_rand_log2(n,:) = mvnrnd(lowGW_log_mean,C_wish_GW(:,:,n),1);
        lowGW_rand_untrans2(n,:) = exp(lowGW_rand_log1(n,:));
        anySaline = any(lowGW_rand_untrans2(n,Cl)>Cl_cutoff) + any(lowGW_rand_untrans2(n,Cl)<Cl_mdl); 
    end; 
end


% MIX THE LOW SALINITY SYNTHETIC GW SAMPLES WITH VARIABLE AMOUNTS
% OF RANDOM SAMPLES OF END-MEMBERS
% Determine the range of percent contamination for each end-member that
% will generate a high salinity dataset with similar mean and SD of Cl in
% observed high salinity GW.
lowGW_logmean_Cl = nanmean(log(lowGW(:,Cl)));
highGW_logmean_Cl = nanmean(log(highGW(:,Cl)));
lowGW_logSD_Cl = nanstd(log(lowGW(:,Cl)));
highGW_logSD_Cl = nanstd(log(highGW(:,Cl)));

r = randn(test_size,1);
logmean_Cl = cell(1,length(endmembers));
p_mean = cell(1,length(endmembers));
p_SD = cell(1,length(endmembers));
p_log = cell(1,length(endmembers));

for n = 1:length(endmembers)
    logmean_Cl{n} = nanmean(log(endmembers{n}(:,Cl)));
    p_mean{n} = (highGW_logmean_Cl - lowGW_logmean_Cl) / (logmean_Cl{n} - lowGW_logmean_Cl);
    p_SD{n} = (highGW_logSD_Cl - lowGW_logSD_Cl) / (logmean_Cl{n} - lowGW_logSD_Cl);
    p_log{n} = ((r * p_SD{n}) + p_mean{n});
end

% Cap the possible percent contamination for different endmembers
% These numbers are set by trial and error
% Might need to revise if end member order or type changes
% See parameters for specifics
for n = 1:length(endmembers)
    if isfinite(mixcap{n})
        p_log{n}(p_log{n}>mixcap{n}) = mixcap{n};
    end
end

% Create random samples of end-members (based on statistics of observed
% data) to mix with low salinity groundwater. First, determine a mean and
% covariance matrix for the end-member populations using observed data.
EM_log = cell(1, length(endmembers));
log_mean = cell(1, length(endmembers));
log_C = cell(1, length(endmembers));
for n = 1:length(endmembers)
    EM_log{n} = log(endmembers{n});
    if any(isnan(EM_log{n}(:)))
        [log_mean{n}, log_C{n}] = ecmnmle(EM_log{n}, 'diagonal', maxI{n});
        log_mean{n} = log_mean{n}';
    else
        log_mean{n} = mean(EM_log{n});
        log_C{n} = cov(EM_log{n}, 'partialrows');       
    end
end

% Create two random sets of synthetic end-member samples (n = test size).
% Use a random covariance matrix pulled from the Inverse Wishart
% distribution. Use second set of random numbers to increase sample size
% for some populations (e.g. septic effluent, see below).
[C_wish, rand_log1, rand_log2] = deal(cell(1,length(endmembers)));
[rand_untrans1, rand_untrans2] = deal(cell(1,length(endmembers)));
for n = 1:length(endmembers)
    for m = 1:test_size
        C_wish{n}(:,:,m) = iwishrnd(log_C{n},size(EM_log{n},1))*(size(EM_log{n},1)-size(log_C{n},1) - 1);
        rand_log1{n}(m,:) = mvnrnd(log_mean{n},C_wish{n}(:,:,m),1);
        rand_log2{n}(m,:) = mvnrnd(log_mean{n},C_wish{n}(:,:,m),1);
    end
    rand_untrans1{n} = exp(rand_log1{n});
    rand_untrans2{n} = exp(rand_log2{n});
end

% Mix random low salinity GW with the random end members. Use
% random percentages created above.  Remove any synthetic high salinity 
% samples with Cl concentrations below the Cl cutoff. (Final arrays of 
% contaminated high salinity groundwater may have <3000 samples)
[GW_log1,GW_log2, GW_mix1, GW_mix2] = deal(cell(1, length(endmembers)));
[pct_mix1, pct_mix2,GW_high1,GW_high2] = deal(cell(1, length(endmembers)));
[pct_high1, pct_high2] = deal(cell(1, length(endmembers)));
for n = 1:length(endmembers)
    for m = 1:test_size
        GW_log1{n}(m,Cl) = p_log{n}(m,1)*rand_log1{n}(m,Cl)+(1-p_log{n}(m,1))*log(lowGW_rand_untrans1(m,Cl));
        GW_log2{n}(m,Cl) = p_log{n}(m,1)*rand_log2{n}(m,Cl)+(1-p_log{n}(m,1))*log(lowGW_rand_untrans2(m,Cl));
    end
    GW_mix1{n} = exp(GW_log1{n});
    GW_mix2{n} = exp(GW_log2{n});
    for m = 1:test_size
        pct_mix1{n}(m,1) = (GW_mix1{n}(m,Cl)-lowGW_rand_untrans1(m,Cl))/(rand_untrans1{n}(m,Cl)-lowGW_rand_untrans1(m,Cl));
        pct_mix2{n}(m,1) = (GW_mix2{n}(m,Cl)-lowGW_rand_untrans2(m,Cl))/(rand_untrans2{n}(m,Cl)-lowGW_rand_untrans2(m,Cl));
    end
    for p = 1:n_sol
        for m = 1:test_size
            GW_mix1{n}(m,p) = pct_mix1{n}(m,1)*rand_untrans1{n}(m,p)+(1-pct_mix1{n}(m,1))*lowGW_rand_untrans1(m,p);
            GW_mix2{n}(m,p) = pct_mix2{n}(m,1)*rand_untrans2{n}(m,p)+(1-pct_mix2{n}(m,1))*lowGW_rand_untrans2(m,p);
        end
    end
    GW_high1{n} = GW_mix1{n}(GW_mix1{n}(:,Cl)>Cl_cutoff,:);
    GW_high2{n} = GW_mix2{n}(GW_mix2{n}(:,Cl)>Cl_cutoff,:);
    pct_high1{n} = pct_mix1{n}(GW_mix1{n}(:,Cl)>Cl_cutoff,:);
    pct_high2{n} = pct_mix2{n}(GW_mix2{n}(:,Cl)>Cl_cutoff,:);
end


% TRAIN DISCRIMINANT ANALYSIS MODEL TO DISTINGUISH HIGHLY SALINE SAMPLES
% CREATED BY VARIOUS END-MEMBERS
% Combine synthetic high salinity datasets for Discriminant Analysis
[EM_size1, EM_size2] = deal(cell(1, length(endmembers)));
for n = 1:length(endmembers)
    EM_size1{n} = size(GW_high1{n},1);
    EM_size2{n} = size(GW_high2{n},1);
end
for n = 1:length(endmembers)
    if n == 1
        start_size = 1;
        end_size = EM_size1{n};
    else
        start_size = 1;
        for m = (n-1):-1:1
            start_size = EM_size1{m} + start_size;
        end
        end_size = start_size + EM_size1{n} - 1;
    end
    cont_gw(start_size : end_size, :) = GW_high1{n};   % Create ann array of chemistry
    cont_gw_type(start_size : end_size, 1) = n;        % Create an array for category    
end

% These lines determine if any of the synthetic high salinity groundwater
% endmember datasets are too small (defined as n < 1500 ) due to the Cl 
% threshold. If they are, additional data points are added
for n = 1:length(endmembers)
    if size(GW_high1{n},1) < 1500
        temp_size = size(cont_gw,1);
        cont_gw(temp_size + 1:temp_size + EM_size2{n}, :) = GW_high2{n};
        cont_gw_type(temp_size + 1:temp_size + EM_size2{n}, 1) = n;
    end
end


% Use backward sequential feature selection to identify the parameters
% needed for the discriminant analysis.  Log transform data for Discriminant Analysis.
Vec=log(cont_gw);
Label=cont_gw_type;
f=@(TrainVec,TrainLabel,TestVec,TestLabel) sum(TestLabel ~= classify(TestVec,TrainVec,TrainLabel));
[bfs, historybfs] = sequentialfs(f,Vec,Label,'direction','backward');

% Write new data vectors that have only the solutes identified using
% sequential selection. 
v=0;
for n=1:n_sol;
    if historybfs.In(end,n)==1; 
        v=v+1;
        cont_gw_sel(:,v)=cont_gw(:,n);
        highGW_sel(:,v)=highGW(:,n);
    end
end

n_sol_sel=nnz(historybfs.In(end,:));  % Reset number of solutes being used

% Standardize log-transformed solute concentrations for synthetic high
% salinity GW
cont_gw_log=log(cont_gw_sel);
cont_gw_z=zscore(cont_gw_log);
% Standardize log-transformed solute concentrations for observed high
% salinity GW, using mean and SD of synthetic data, so that scores are on
% same scale for observed and synthetic data
highGW_z = NaN(size(highGW_sel,1),n_sol_sel);
for n=1:n_sol_sel
    highGW_z(:,n)=((log(highGW_sel(:,n)))-nanmean(cont_gw_log(:,n)))./nanstd(cont_gw_log(:,n));
end

% Create Discriminant Analysis classifier and then apply to training data
% to generate confusion matrix. Use a 10-fold cross-validation.
c = cvpartition(cont_gw_type,'k',10);
order = unique(cont_gw_type); % Order of the group labels
f = @(xtr,ytr,xte,yte)confusionmat(yte,classify(xte,xtr,ytr),'order',order);
conf_matrix = crossval(f,cont_gw_z,cont_gw_type,'partition',c);
conf_matrix = reshape(sum(conf_matrix),length(endmembers),length(endmembers));
conf_matrix_per = NaN(length(endmembers),length(endmembers));
for n=1:length(endmembers)
    for m=1:length(endmembers)
        conf_matrix_per(n,m)=100*conf_matrix(n,m)/sum(conf_matrix(n,:));
    end
end
[prdctn,err,P,logp,coeff] = classify(cont_gw_z,cont_gw_z,cont_gw_type);

% Apply Discriminant Analysis classifier to observed high salinity GW data,
% which are all unknowns
[prdctn2,err2,P2,logp2,coeff2] = classify(highGW_z,cont_gw_z,cont_gw_type,'linear',prior_prob);

% Compute "scores" for all synthetic high salinity data
set_size=size(cont_gw_z,1);
for t=1:(length(endmembers) - 1)
    Kof = coeff(t,t+1).const;
    L = coeff(t,t+1).linear; 
    L=L';
    loadings(:,t)=L;    % to output loadings and coefficients used to compute scores
    coefficients(:,t)=Kof;
    % Compute scores as Kof (a constant) + L*[all solutes] 
    score_part(1:set_size,1:n_sol_sel)=zeros;  % Preallocate matrix
    scores(1:set_size,t)=zeros;  % Preallocate matrix
    for n=1:set_size
        for r=1:n_sol_sel;
            score_part(n,r)=L(1,r).*(cont_gw_z(n,r));
        end
        scores(n,t)=Kof+sum(score_part(n,1:n_sol_sel));
    end
    % Compute correlations between scores and synthetic data variables
    corr_data_scores(t,1:n_sol_sel)=zeros;  % Preallocate matrix
    for n=1:n_sol_sel;
        corr_data_scores(t,n)=corr(cont_gw_sel(:,n),scores(:,t));
    end
end

% Compute "scores" for all observed high salinity data 
set_size=size(highGW,1);
for t=1:(length(endmembers) - 1)
    Kof = coeff(t,t+1).const;
    L = coeff(t,t+1).linear; 
    L=L';
    % Compute scores as Kof (a constant) + L*[all solutes] 
    score_part2(1:set_size,1:n_sol_sel)=zeros;  % Preallocate matrix
    scores2(1:set_size,t)=zeros;  % Preallocate matrix
    for n=1:set_size
        for r=1:n_sol_sel;
            score_part2(n,r)=L(1,r).*(highGW_z(n,r));
        end
        scores2(n,t)=Kof(1,1)+sum(score_part2(n,1:n_sol_sel));
    end
end

% Create an array reporting results for observed high salinity GW data
set_size=size(highGW,1);
g=0;
for n=1:set_size
    if P2(n,1)>0
        g=g+1;
        SWIFT_type(g,1)=prdctn2(n,1);                   % predicted class
        SWIFT_type(g,2:1+size(P2,2))=P2(n,1:end);       % probabilities
        SWIFT_type(g,4+size(P2,2))=max(P2(n,:));        % max probability
        SWIFT_type_chem(g,:)=highGW(n,:);               % chemistry
        SWIFT_type_scores(g,:)=scores2(n,:);            % classifier scores
        SWIFT_type_name(g,1)=highGW_name(n,1);          % sample names
        SWIFT_type_exp(g,1)=highGW_exp(n,1);
    end
end

% Compute Pearson Correlation Coeff between scores and synthetic high 
% salinity groundwater, for each contaminant.
for n = 1:length(endmembers)
    for t = 1:(length(endmembers) - 1)
        corr_EM_scores{n}(t,1:n_sol_sel)=zeros;  % Preallocate matrix
        for m=1:n_sol_sel;
            corr_EM_scores{n}(t,m)=corr(cont_gw_sel(cont_gw_type==n,m),scores(cont_gw_type==n,t));
        end
    end
    corr_EM_scores{n} = corr_EM_scores{n}';
end
    

% Output model results of interest
output.headers_used=headers;
output.rawdata.endmembers=endmembers;
output.lowGW_observed=lowGW;
output.lowGW_synthetic=lowGW_rand_untrans1;
output.highGW_observed=highGW;
output.highGW_synthetic.chemistry=cont_gw;
output.highGW_synthetic.type=cont_gw_type;
output.constants=loadings;
output.coefficients=coefficients;
output.conf_matrix=conf_matrix;
output.conf_matrix_percents=conf_matrix_per;
output.corr_synth_data_and_scores=corr_data_scores';
output.corr_synth_data_and_scores_EM = corr_EM_scores;
output.bckwd_feat_sel=bfs;
output.SWIFT_type.type_and_probabilities=SWIFT_type;
output.SWIFT_type.chem=SWIFT_type_chem;
output.SWIFT_type.scores=SWIFT_type_scores;
output.SWIFT_type.names=SWIFT_type_name;
output.SWIFT_type.exp=SWIFT_type_exp;

% Extra output for figures; can disregard
output.lowGW = lowGW;
output.highGW = highGW;
output.GW_high1 = GW_high1;
output.scores = scores;
output.cont_gw = cont_gw;
output.endmembers_for_figures = param_input.endmembers;

end