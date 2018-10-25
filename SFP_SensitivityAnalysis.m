% Nathaniel Chien
% Start Date: April 20th, 2017
% Linear Discriminant Analysis Project
% Sensitivity Analysis Code

% Objectives
%  - Test the sensitivity of model error to variable numbers of solutes
%  - Test the sensitivity of model error to prior probabilities
%  - Test the sensitivity of model error to mixing caps

%   Params: array of parameters in the following order
%       1)  Solute starting index
%       2)  Solute ending index
%       3)  Location of chloride in list of solutes
%       4)  Cell array of covariance matrix iterations
%       5)  Cell array of mixing percentages capped
%       6)  Synthetic data test size
%       7)  Chloride Cutoff
%       8)  Chloride detection limit

% Name the excel data file being used and pick out end member
data_file = 'RevisedData2_OrganicWaste.xlsx';
FW = xlsread(data_file,'FW','D:M');
RS = xlsread(data_file,'RS','D:M');
ORG = xlsread(data_file,'SEP','D:M');
% Read in excel data file for groundwater and endmember datasets
[ ~, headers ] = xlsread(data_file, 'HEADERS', 'A1:N1');
[ GW_exp, ~, GW_name ] = xlsread(data_file,'GW','A:E');
GW_name = GW_name(2:62,:);
GW = xlsread(data_file,'GW','F:O');   
% Endmember data array
em_Data = {headers,GW,GW_name,GW_exp,FW,RS,ORG};


% Sensitivity Analysis Below

% Parameters I have been using previously
param1 = {2,10,5,{1000,30,200},{NaN,0.38,0.26},3000,20,0.2,{0.33,0.33,0.33}};
test_run = SFP_parameters(em_Data,param1);
test_out = SFP_classification(test_run);
% Calculate error
test_error = 0;
for i=1:length(test_out.SWIFT_type.exp)
    if (test_out.SWIFT_type.exp(i) ~= test_out.SWIFT_type.type_and_probabilities(i,1))
        test_error = test_error + 1;
    end
end
test_error = test_error / length(test_out.SWIFT_type.exp) * 100;

% Test error related to prior probabilities
prior_temp = NaN(1000,4);
% Determine random prior probabilities that sum to 1
rand_prior = rand(1000,3);
rand_prior_sum = rand_prior(:,1) + rand_prior(:,2) + rand_prior(:,3);
rand_prior(:,1) = rand_prior(:,1) ./ rand_prior_sum;
rand_prior(:,2) = rand_prior(:,2) ./ rand_prior_sum;
rand_prior(:,3) = rand_prior(:,3) ./ rand_prior_sum;

for i=1:length(rand_prior)
    % Set parameters
    param1{9} = {rand_prior(i,1),rand_prior(i,2),rand_prior(i,3)};
    temp_run = SFP_parameters(em_Data,param1);
    temp_out =  SFP_classification(temp_run);
    % Calculate error
    temp_error = 0;
    for cnt=1:length(temp_out.SWIFT_type.exp)
        if (temp_out.SWIFT_type.exp(cnt) ~= temp_out.SWIFT_type.type_and_probabilities(cnt,1))
            temp_error = temp_error + 1;
        end
    end
    % Output 
    prior_temp(i,1) = rand_prior(i,1);
    prior_temp(i,2) = rand_prior(i,2);
    prior_temp(i,3) = rand_prior(i,3);
    prior_temp(i,4) = temp_error;    
end
% Reset parameters
param1 = {2,10,5,{1000,30,200},{NaN,0.38,0.26},3000,20,0.2,{0.33,0.33,0.33}};

% Look at prior probabiltiies at error
[err_prior, indx_prior] = sort(prior_temp(:,4));
resort_prior = prior_temp(indx_prior,:);
 

% Extract samples with prior probabilities in correct range
prior_range = NaN(length(resort_prior),4);
count = 1;
for i=1:length(resort_prior)
    if ((resort_prior(i,1) < 0.48)&&(resort_prior(i,1)>0.18))
        if ((resort_prior(i,2) < 0.48)&&(resort_prior(i,2)>0.18))
            if ((resort_prior(i,3) < 0.48)&&(resort_prior(i,3)>0.18))
                prior_range(count,:) = resort_prior(i,:);
                count = count + 1;
            end
        end
    end
end

% Run best trial!!!
param1 = {2,10,5,{1000,30,200},{NaN,0.38,0.26},3000,20,0.2,{0.33,0.33,0.33}};
param1{9} = {0.0417830711057856,0.153022283686667,0.805194645207547};
best_run = SFP_parameters(em_Data,param1);
best_out =  SFP_classification(best_run);
% Calculate error
best_error = 0;
for i=1:length(best_out.SWIFT_type.exp)
    if (best_out.SWIFT_type.exp(i) ~= best_out.SWIFT_type.type_and_probabilities(i,1))
        best_error = best_error + 1;
    end
end

% Run trial with 10% formation water prior for assessment
param1 = {2,10,5,{1000,30,200},{NaN,0.38,0.26},3000,20,0.2,{0.33,0.33,0.33}};
param1{9} = {0.1,0.2,0.7};
best_run2 = SFP_parameters(em_Data,param1);
best_out2 =  SFP_classification(best_run2);
% Calculate error
best_error2 = 0;
for i=1:length(best_out2.SWIFT_type.exp)
    if (best_out2.SWIFT_type.exp(i) ~= best_out2.SWIFT_type.type_and_probabilities(i,1))
        best_error2 = best_error2 + 1;
    end
end

% Look at sensitivity to mixing caps
mix_temp = NaN(1000,4);
rand_mix = rand(1000,3);
for i=1:length(rand_mix)
    % Set parameters
    param1{5} = {rand_mix(i,1),rand_mix(i,2),rand_mix(i,3)};
    temp_run = SFP_parameters(em_Data,param1);
    temp_out =  SFP_classification(temp_run);
    % Calculate error
    temp_error = 0;
    for cnt=1:length(temp_out.SWIFT_type.exp)
        if (temp_out.SWIFT_type.exp(cnt) ~= temp_out.SWIFT_type.type_and_probabilities(cnt,1))
            temp_error = temp_error + 1;
        end
    end
    % Output 
    mix_temp(i,1) = rand_mix(i,1);
    mix_temp(i,2) = rand_mix(i,2);
    mix_temp(i,3) = rand_mix(i,3);
    mix_temp(i,4) = temp_error;    
end
% Reset parameters
param1 = {2,10,5,{1000,30,200},{NaN,0.38,0.26},3000,20,0.2,{0.33,0.33,0.33}};

% Look at mixing percentages error
[err_mix, indx_mix] = sort(mix_temp(:,4));
resort_mix = mix_temp(indx_mix,:);

% Look at no mixing cap
no_mix_run = SFP_parameters(em_Data,param1);
no_mix_out = SFP_classification(no_mix_run);
% Calculate error
no_mix_error = 0;
for i=1:length(test_out.SWIFT_type.exp)
    if (no_mix_out.SWIFT_type.exp(i) ~= no_mix_out.SWIFT_type.type_and_probabilities(i,1))
        no_mix_error = no_mix_error + 1;
    end
end
