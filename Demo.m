

%%%% Demo code to execute wonparafac %%

% Note that the code requires tensortoolbox to be on the path environment.
% The tensor toolbox is available at: https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
% You can add path by addpath command.



%% load GDSC 1000 data (1815 genes)  %%
load Demo.mat

% X is GDSC tensor (1815 gene by 935 sample by 5 data type)
% 5 data types = GE(+), GE(-), MT, CN(+), CN(-)
% gene_names contains names of genes in X

%%% select subset of genes (100 genes) %%%
rand('seed', 1111);
gene_select = randsample(size(X,1), 100);

X_use = X(gene_select,:,:);

%%  configure options  %%

options = struct();
options.tol = 1.0e-7;
options.init = 'nvecs';
optoins.printitn = 1000;

% orthogonal constraint
options.orthogonal = [.2, 0, 0]; % orthogonal constraint only on gene mode (strength = 0.2)

% weight setting
var = 1./sum(sum(X_use.^2,1),2);
var = var./sum(var(:));
options.weight_W = ones(size(X_use)).*repmat(var, [size(X_use,1), size(X_use,2)]);



% get outcome across the number of basis / strength of orthogonal constraint
Nbasis = [10:10:200];
Ortho = [0,0.2,0.5,1];


%% Execute WON-PARAFAC %%
Factors = cell(0);
fits = cell(0);
final_fit = [];

for o=1:length(Ortho)
    options.orthogonal = [Ortho(o),0,0];
    for n=1:length(Nbasis)
        [Fs, fit_info] = wonparafac(tensor(X_use), Nbasis(n), options);
        Factors{n,o} = Fs;
        fits{n,o} = fit_info;
        final_fit(n,o) = fits{n,o}.fits(end);
    end
end

%% Plot the outcome

fig = plot(Nbasis, final_fit);
title('Performanc of WON-PARAFAC with varying number of factors/constraint')
xlabel('Number of factors')
ylabel('Explained variation (reconstruction performance)')

xticks(Nbasis)
xticklabels(Nbasis)
legends = cell(0);
for o=1:length(Ortho)
    legends{o} = ['orthogonal constraint = ', num2str(Ortho(o))];
end

legend(legends, 'Location', 'northwest')
saveas(gcf, 'Demo_plot', 'png')
