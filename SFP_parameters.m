function [ output ] = SFP_parameters( endmember_data, params )

% Function reads in excel data file given as input
% Also sets parameters for 'SFP_classification.m'
% Note that the Excel file must be in the same format as tables in the 
% Supplementary Information
%
% INPUT: 
%   - endmember_data: array of end member data
%   - params: array of parameters in the following order
%       1)  Solute starting index
%       2)  Solute ending index
%       3)  Location of chloride in list of solutes
%       4)  Cell array of covariance matrix iterations
%       5)  Cell array of mixing percentages capped
%       6)  Synthetic data test size
%       7)  Chloride Cutoff
%       8)  Chloride detection limit
%       9)  Prior probability, numeric vector
%       
%   EX) SFP_parameters('File_name',[3,1,7,])
%
% OUTPUT: Data arrays formated for input to SFP_classification.m

% Read end member data in
GW = endmember_data{2};
FW = endmember_data{5};
RS = endmember_data{6};
ORG = endmember_data{7};

% Determines which solutes to pick from entire GW and EM set
% GW = GW(:,1:9) is default/all 9 main solutes from Lautz et al. (2014)
GW = GW(:,params{1}:params{2});
FW = FW(:,params{1}:params{2});
RS = RS(:,params{1}:params{2});
ORG = ORG(:,params{1}:params{2});
endmembers = {FW, RS, ORG};

% Note number and order of the solutes. Need to define where Cl is.
% I = 1; Na = 2; K = 3; Mg = 4; Ca = 5; 
% Cl = 6; Br = 7; Sr = 8; Ba = 9; SO4 = 10
Cl = params{3};
n_sol = params{2} - params{1} + 1;        % Sets # of solutes in the input data

% Salinity parameters
test_size = params{6};
Cl_cutoff = params{7};
Cl_mdl = params{8};

% Covariance parameters (maximum iterations for ecmnmle fxn)
maxI = params{4};

% Mixing cap: Cap the possible percent contamination for different 
% endmembers (these numbers set by trial and error). Target is 20%
mixcap = params{5};

% Set prior probabilities
prior = cell2mat(params{9});

% Outputs
output.test_size = test_size;         % max size of synthetic datasets
output.Cl_cutoff = Cl_cutoff;         % Cl cutoff for high v. low salinity
output.Cl_mdl = Cl_mdl;               % Chloride Detection Limit
output.n_sol = n_sol;                 % number of solutes
output.Cl = Cl;                       % chlorine number     
output.headers =  endmember_data{1};  % headers for solute columns
output.GW = GW;                       % groundwater data only
output.GW_name =  endmember_data{3};  % GW identification
output.GW_exp =  endmember_data{4};   % expected GW values
output.endmembers = endmembers;       % number of endmembers
output.maxI = maxI;
output.mixcap = mixcap;
output.prior = prior;

end


