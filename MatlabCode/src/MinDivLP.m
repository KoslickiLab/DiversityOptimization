%% import the data
% set variables
small_k = 4;  % smaller k-mer size
large_k = 6;  % larger k-mer size

A_k = load(sprintf('../data/97_otus_subset.fasta_A_%d.mat', large_k));
A_k_large = A_k.A_k;

A_k = load(sprintf('../data/97_otus_subset.fasta_A_%d.mat', small_k));
A_k_small = A_k.A_k;

clear('A_k')  % get rid of unneeded variable

%% sub-select the data so things run quickly

cols_vs_rows = 3;  % fix 3-times more columns than rows
num_species = cols_vs_rows*4^small_k;  % reduce the number of columns of 
% the sensing matrix so pictures will be generated in a reasonable amount
% of time.
A_k_small = A_k_small(:, 1:num_species);
A_k_large = A_k_large(:, 1:num_species);

%% Generate the B matrix
B = (A_k_small > 0);

%% Set some variables
q = .1;  % fixed, small q value s.t. 0<q<1
support_size = start:step_size:max_support;  % number of non-zero entries in the simulated ground truth

%% create the simulated ground truth
% create
supp = datasample(1:num_species, support_size, 'Replace', false);  % location of the support
true_x = zeros(num_species,1);  % the true x vector we are trying to reconstruct
true_x(supp) = rand(suppSize,1);  % populate with random data
true_x = true_x./sum(true_x);  % normalize to be a probability vector

y_small_true = A_k_small*true_x;
slop_factor = 2;  % this is for the non-regularized MinDivLP on noisy data
noise_eps = .00001;
y_small_noise = A_k_small*true_x + noise_eps*abs(rand(size(A_k_small,1),1));
y_small_noise = y_small_noise./sum(y_small_noise);

% optimization parameters
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off','Diagnostics','off', 'ConstraintTolerance', .0000001, 'OptimalityTolerance', .0000001);

%% Noisless computations

% Using matlab built in linprog: too slow
%[x_l1, ~, ~, ~, ~] = linprog(ones(1,num_species), [], [], A_k_small, A_k_small*true_x, zeros(1,num_species), ones(1,num_species), options);  % if you do not have a Gurobi license, you can use the built in linprob


for i=1:length(h_sizes)
            y_s{i} = A_hs{i}*true_x;
        end
        
        slop_factor = 1;
        noise_eps = .00001;
        y_noise = A_k_small*true_x + noise_eps*abs(rand(size(A_k_small,1),1));
        y_noise = y_noise./sum(y_noise);
        
        % optimization parameters
        options = optimoptions('linprog','Algorithm','dual-simplex','Display','off','Diagnostics','off', 'ConstraintTolerance', .0000001, 'OptimalityTolerance', .0000001);

        % Unweighted L1 optimization
        % no noise
        %[x_l1, ~, ~, ~, ~] = linprog_gurobi(ones(1,num_species), [], [], A_k_small, A_k_small*true_x, zeros(1,num_species), ones(1, num_species), options);
        
        % with noise
        %[x_l1, ~, ~, ~, ~] = linprog_gurobi(ones(1,num_species), [A_k_small; -A_k_small], [y_noise+slop_factor*noise_eps; -y_noise+slop_factor*noise_eps], [], [], zeros(1,num_species), ones(1, num_species), options);
        
        % matlab linprog too slow
        %[x_l1, ~, ~, ~, ~] = linprog(ones(1,num_species), [], [],A_k_small, A_k_small*true_x, zeros(1,num_species), ones(1,num_species), options);  % if you do not have a Gurobi license, you can use the built in linprob
        
        %x_l1 = x_l1./sum(x_l1);  % normalize
        %temp_l1(rep) = sum(abs(true_x - x_l1));  % store the error
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Temp little test
        % Weighted quikr/regularized Mindiv (store it in x_l1 just so I
        % don't need to recode things
        lambda = 10000;
        i = length(h_sizes);
        f = 1./(B_transpose_hs{i}*y_s{i}).^(1-q);
        x_l1 = lsqnonneg([f'; lambda*A_k_small], [0;lambda*y_noise]);
        x_l1 = x_l1./sum(x_l1);  % normalize
        temp_l1(rep) = sum(abs(true_x - x_l1));  % store the error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Conclusion: This appears to be the best!!!!!!!
        % Go with this for real-world data
        % For a possible python clone (no idea how it performs), see:
        % https://github.com/matthew-brett/diffusion_mri/blob/master/Python/lsqnonneg.py
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % quikr
        lambda = 10000;
        tic;
        % no noise
        %x_q = lsqnonneg([ones(1, num_species); lambda*A_k_small], [0;lambda*A_k_small*true_x]);
        % with noise
        x_q = lsqnonneg([ones(1, num_species); lambda*A_k_small], [0;lambda*y_noise]);
        x_q = x_q/sum(x_q);
        %x_q = ones(num_species, 1);
        temp_q(rep) = sum(abs(true_x - x_q));

        % weighted optimization, over all h sizes
        for i=1:length(h_sizes)
            f = 1./(B_transpose_hs{i}*y_s{i}).^(1-q);
            % no noise
            %[x_star, ~, ~, ~, ~] = linprog_gurobi(f, [], [], A_k_small, A_k_small*true_x, zeros(1,num_species), ones(1, num_species), options);  % again, can substitute linprog
            % with noise
            [x_star, ~, ~, ~, ~] = linprog_gurobi(f, [A_k_small; -A_k_small], [y_noise+slop_factor*noise_eps; -y_noise+slop_factor*noise_eps], [], [], zeros(1,num_species), ones(1, num_species), options);  % again, can substitute linprog
            x_star = x_star./sum(x_star);
            temp_weighteds(i, rep) = sum(abs(true_x - x_star));
            temp_posteriori_test(i, rep) = sum(abs(A_k_large*true_x - A_k_large*x_star));
        end
    end
    unweighted_errors(:, suppSizeInd) = temp_l1;
    quikr_errors(:, suppSizeInd) = temp_q;
    weighted_errors(:, :, suppSizeInd) = temp_weighteds;
    posteriori_test(:, :, suppSizeInd) = temp_posteriori_test;
    %ppm.increment()
end
fprintf('Finished\n')
toc

%% L1 norm error plot: plot of L1 error between true x and reconstructed x 
% (averaged over the number of replicates) as a function of support size ||x||_0
colors = linspecer(2+length(h_sizes));
fs = 15;
set(groot,'defaultAxesColorOrder', [0 0 0], 'DefaultAxesLineStyleOrder','-|--|:|-.|-*')
line_width = 2;
figure();
hold on
plot(support_sizes, mean(unweighted_errors), 'LineWidth', line_width, 'Color', colors(1,:))
plot(support_sizes, mean(quikr_errors), 'LineWidth', line_width, 'Color', colors(2,:))
for i=1:length(h_sizes)
     plot(support_sizes, mean(squeeze(weighted_errors(i,:,:))), 'LineWidth', line_width, 'MarkerSize', 4, 'Color', colors(i+2,:))
end
ylabel('Mean L1 error', 'FontSize', fs)
xlabel('Support size', 'FontSize', fs)
title(sprintf('Reconstruction performance: L1 norm using k = %d', small_k))
legends = {};
legends{1} = 'Unweighted';
legends{2} = 'Quikr';
for i=1:length(h_sizes)
       legends{i+2} = sprintf('Weighted, h = %d', h_sizes(i));
end
lgd = legend(legends{:});
lgd.FontSize = fs;
%set(gca,'DataAspectRatio',[15 1 1])
%x= 1.5;
%set(gcf,'Position',[x*100 x*100 x*500 x*500])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%saveas(gcf, 'Figures/L1NormErrorMultiK.png')

%% Percent recovered plot:
% The percentage (over all replicates) of successful recoveries as a 
% function of the support size
thresh = 1e-1;  % anything L1 error smaller than this is considered "successful"

colors = linspecer(2+length(h_sizes));
fs = 15;
set(groot,'defaultAxesColorOrder', [0 0 0], 'DefaultAxesLineStyleOrder','-|--|:|-.|-*')
line_width = 2;
figure();
hold on
plot(support_sizes, mean(unweighted_errors<thresh), 'LineWidth', line_width, 'Color', colors(1,:))
plot(support_sizes, mean(quikr_errors<thresh), 'LineWidth', line_width, 'Color', colors(2,:))
for i=1:length(h_sizes)
     plot(support_sizes, mean(squeeze(weighted_errors(i,:,:))<thresh), 'LineWidth', line_width, 'MarkerSize', 4, 'Color', colors(i+2,:))
end
ylabel('Percent recovered', 'FontSize', fs)
xlabel('Support size', 'FontSize', fs)
legends = {};
legends{1} = '$\ell_1$';
legends{2} = 'Quikr';
for i=1:length(h_sizes)
       legends{i+2} = sprintf('h = %d', h_sizes(i));
end
lgd = legend(legends{:},'Interpreter','latex');
lgd.FontSize = fs;
set(gca,'DataAspectRatio',[20 1 1])
%x= 1.5;
%set(gcf,'Position',[x*100 x*100 x*500 x*500])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%saveas(gcf, 'Figures/UnweightedVsWeightedMultiK.png')

%% Demonstrate the a posteriori test
thresh = 1e-5;  % anything L1 error smaller than this is considered "successful"

colors = linspecer(1+length(h_sizes));
fs = 15;
set(groot,'defaultAxesColorOrder', [0 0 0], 'DefaultAxesLineStyleOrder','-|--|:|-.|-*')
line_width = 2;
fig = figure();
set(fig,'defaultAxesColorOrder',[colors(1,:); colors(2,:)]);
hold on
i = length(h_sizes);
x = support_sizes;
yyaxis left
plot(x, mean(squeeze(weighted_errors(i,:,:))<thresh), 'LineWidth', line_width, 'MarkerSize', 4, 'Color', colors(1,:))
yyaxis right
plot(x, mean(squeeze(posteriori_test(i,:,:))), 'LineWidth', line_width, 'MarkerSize', 4, 'Color', colors(2,:))
yyaxis left
ylabel('Percent recovered', 'FontSize', fs)
xlabel('Support size', 'FontSize', fs)
yyaxis right
ylabel('$||A^{(h)}x - y^{(h)}||_1$', 'FontSize', fs, 'Interpreter','latex')
%set(fig,'DataAspectRatio',[112.5,1,2.5])
set(gca,'DataAspectRatio',[112.5,1,2.5])
set(gca,'PlotBoxAspectRatio',[1,0.533333333333333,0.533333333333333])
%x = 1.5;
%set(gcf,'Position',[x*100 x*100 x*500 x*500])
ax = gca;
%outerpos = ax.OuterPosition;
%ti = ax.TightInset; 
%left = outerpos(1) + ti(1);
%bottom = outerpos(2) + ti(2);
%ax_width = outerpos(3) - ti(1) - ti(3);
%ax_height = outerpos(4) - ti(2) - ti(4);
%ax.Position = [left bottom ax_width ax_height];
%saveas(gcf, 'Figures/aposterioriMultiK.png')
%% export the data
%if ~isfile('CooccurenceReproduciblesMultiData.mat')
%    save('CooccurenceReproduciblesMultiData.mat')
%end
