% Based on the code by Elizabeth Newman, available at
% https://github.com/elizabethnewman/projected-products

% Edited by: Susana Lopez Moreno, Pusan National University, 2025

clear; clc;

%% Reproducibility seed
rng(2025);

%% Load data

hcube = hypercube('indian_pines.dat');
A     = hcube.DataCube;

% normalize
A = A / fronorm(A);

saveDir  = './results_hyperspectral_injective';
fileName = 'results';
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
if ~exist([saveDir,'orig'], 'dir'), mkdir([saveDir,'/orig']); end

% get sizes
[m,n,p]  = size(A);
nrmA        = fronorm(A);
storeA      = m * n * p;

%% Form matrices

% Injective matrices
% JL
oversamp = 2;
target  = oversamp * p;          % originally 440
N       = 2^nextpow2(target);    % rounds up to 512
D       = diag(sign(randn(N,1)));
Sfull   = (1/sqrt(p))*hadamard(N)*D;   % N×N
rows    = randperm(N, p);
M       = Sfull(rows, :)';            % now N×p
%q       = N;                          % use this N for cost
Mp = (M.'*M)\(M.');

%DCT
q = 2*p;
Cfull = dctmtx(q);             % orthogonal q x q
M_DCT = Cfull(:,1:p);          % q x p
M_DCT_inv = (M_DCT' * M_DCT) \ (M_DCT');

%DST
Sdst = zeros(q,q);
for a = 1:q
    for b = 1:q
        Sdst(a,b) = sqrt(2/(q+1)) * sin(pi*a*b/(q+1));
    end
end

M_DST = Sdst(:,1:p);           % q x p
M_DST_inv = (M_DST' * M_DST) \ (M_DST');

% Surjective/invertible matrices for comparison: U_3 and I
[Z,S,~] = svd(modeUnfold(A,3),'econ');

I = eye(p);

% List of matrices and names
% Q = {'JL', M};
% Q_inv = {'JL_inv', Mp};
Q = {'JL', M, 'U3', Z', 'I', I, 'DCT', M_DCT, 'DST', M_DST};
Q_inv = {'JL_inv', Mp, 'U3_inv', Z, 'I_inv', I', 'DCT_inv', M_DCT_inv, 'DST_inv', M_DST_inv};

% Headers
% headers = {'k','s','rel_err','comp_ratio','rel_comp_ratio'};
headers = {'k','s','rel_err','storage','comp_ratio'};

%% *_M-SVD

for i = 1:2:length(Q)
    
    if ~exist([saveDir,Q{i}], 'dir'), mkdir([saveDir,'/',Q{i}]); end
    if ~exist([saveDir,Q{i},'/img'], 'dir'), mkdir([saveDir,'/',Q{i},'/img']); end
    
    % Push‐forward & facewise SVD
    AHat             = modeProduct(A,Q{i+1});
    [UHat,SHat,VHat] = facewiseSVD(AHat);

    results = [];
    q = size(AHat,3);
    
    % Truncation loop
    for k = 1:min(size(A,1:2))
        % Facewise multiplication of Uhat,Shat,Vhat
        AkHat = facewise(UHat(:,1:k,:),facewise(SHat(1:k,1:k,:), tran(VHat(:,1:k,:))));

        % Pull-back
        Ak  = modeProduct(AkHat,Q_inv{i+1});
        for s = 1:size(A,3)

            % Slice the approximation
            approx = Ak(:,:,1:s);

            % Compute error and storage
            err  = fronorm(A(:,:,1:s) - approx)/fronorm(A(:,:,1:s));
            comp = m * k * s + n * k * s;

            if ismember(Q{i}, {'JL','DCT','DST'})
                comp = comp + q * p;
            elseif strcmp(Q{i},'I')
                comp = comp;
            else
                comp = comp + p * p;
            end
            
            % Store results
            results = cat(1,results,[k,s,err, comp,storeA / comp]);
 
        end
         % Print progress
        fprintf('%s: completed k = %d\n', Q{i}, k);
    end
    
    % save results for fixed Q as csv file
    T = array2table(results, 'VariableNames',headers);
    writetable(T,[saveDir,'/',Q{i},'/',fileName,'.csv'])

end


%% ------------------------------------------------------------------------
%  RGB snapshots + comparison plots
%% ------------------------------------------------------------------------

% Parameters
kpPairs = [1,10; 5,10; 20,50; 20,100];   % [k, p] pairs to snapshot
ps_list = [10,50,100,220];               % for the error‐vs‐compression curves
[R,G,B]  = deal(26,16,8);                % bands used for RGB

% reload the original (normalized) cube and matrices
hcube = hypercube('indian_pines.dat');
A     = hcube.DataCube;      
A     = A / fronorm(A);
[m,n,p] = size(A);

% Injective matrices
% JL
oversamp = 2;
target  = oversamp * p;               % originally 440
N       = 2^nextpow2(target);         % rounds up to 512
D       = diag(sign(randn(N,1)));
Sfull   = (1/sqrt(p))*hadamard(N)*D;  % N×N
rows    = randperm(N, p);
M       = Sfull(rows, :)';            % now N×p
q       = N;                          % use this N for cost
Mp = (M.'*M)\(M.');

%DCT
q = 2*p;
Cfull = dctmtx(q);                    % orthogonal q x q
M_DCT = Cfull(:,1:p);                 % q x p
M_DCT_inv = (M_DCT' * M_DCT) \ (M_DCT');

%DST
Sdst = zeros(q,q);
for a = 1:q
    for b = 1:q
        Sdst(a,b) = sqrt(2/(q+1)) * sin(pi*a*b/(q+1));
    end
end

M_DST = Sdst(:,1:p);                  % q x p
M_DST_inv = (M_DST' * M_DST) \ (M_DST');

% Surjective/invertible matrices for comparison: U_3 and I
[Z,S,~] = svd(modeUnfold(A,3),'econ');


% List of matrices and names
Q = {'JL', M, 'U3', Z', 'DCT', M_DCT, 'DST', M_DST};
Q_inv = {'JL_inv', Mp, 'U3_inv', Z, 'DCT_inv', M_DCT_inv, 'DST_inv', M_DST_inv};

%% 1) RGB snapshots

for i = 1:2:length(Q)
    outImg  = fullfile(saveDir,'/',Q{i},'/img');

    if ~exist(outImg,'dir'), mkdir(outImg); end

    Ah = modeProduct(A,Q{i+1});
    [Uhat,Shat,Vhat] = facewiseSVD(Ah);

    % original RGB
    fig = figure('Visible','off');
    imagesc(rescale(A(:,:, [R G B])));
    axis off; pbaspect([n m 1]);
    exportgraphics(fig, fullfile(outImg,'orig_RGB.png'));
    close(fig);

    % approximations
    for idx = 1:size(kpPairs,1)
      k  = kpPairs(idx,1);
      ps = kpPairs(idx,2);

      % reconstruct rank‐k, pull back all p bands
      Ahat_k     = facewise( ...
                      Uhat(:,1:k,:), ...
                      facewise(Shat(1:k,1:k,:), tran(Vhat(:,1:k,:))) ...
                   );

      approx_full= modeProduct(Ahat_k,Q_inv{i+1});        % m×n×p

      pad_cube = zeros(m, n, p);
      pad_cube(:,:,1:size(approx_full,3)) = approx_full;

      rgb = rescale(pad_cube(:,:, [R G B]));

      fig = figure('Visible','off');
      imagesc(rgb);
      axis off; pbaspect([n m 1]);
      fname = sprintf('img_RGB_k%02d_p%03d.png',k,ps);
      exportgraphics(fig, fullfile(outImg,fname));
      close(fig);
    end
end

%% 2) Plot relative error vs k for all JL and U3

% 1) Read CSV results
projectDir = pwd;  
cd(projectDir);
T_JL  = readtable('./results_hyperspectral_injective/JL/results.csv');
T_DCT  = readtable('./results_hyperspectral_injective/DCT/results.csv');
T_DST  = readtable('./results_hyperspectral_injective/DST/results.csv');
T_U3  = readtable('./results_hyperspectral_injective/U3/results.csv');

% 2) Fix the slice‐counts
s_JL  = 220;  
s_DCT = 220; 
s_DST = 220; 
s_U3 = 220;   

% 3) Select rows
idx_JL = T_JL.s == s_JL;
idx_DCT = T_DCT.s  == s_DCT; 
idx_DST = T_DST.s  == s_DST; 
idx_U3 = T_U3.s  == s_U3; 

% 4) Extract k and relative error
k_JL   = T_JL.k(idx_JL);
err_JL = T_JL.rel_err(idx_JL);

k_DCT   = T_DCT.k(idx_DCT);
err_DCT = T_DCT.rel_err(idx_DCT);

k_DST   = T_DST.k(idx_DST);
err_DST = T_DST.rel_err(idx_DST);

k_U3   = T_U3.k(idx_U3);
err_U3 = T_U3.rel_err(idx_U3);

% 5) Plot
figure; hold on;
plot(k_JL, err_JL, '-o', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','JL');
plot(k_DCT, err_DCT, '-^', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','DCT');
plot(k_DST, err_DST, '-d', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','DST');
plot(k_U3, err_U3, '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','U_3');
xlabel('Truncation k');
ylabel('Relative error');
legend('Location','best','FontSize',14);
grid on;

% Inset
ax2 = axes('Position',[0.55 0.22 0.28 0.28]);
box on; hold on;
plot(ax2, k_JL, err_JL, '-o', 'LineWidth',1.2, 'MarkerSize',4);
plot(ax2, k_DCT, err_DCT, '-^', 'LineWidth',1.2, 'MarkerSize',4);
plot(ax2, k_DST, err_DST, '-d', 'LineWidth',1.2, 'MarkerSize',4);
plot(ax2, k_U3, err_U3, '-s', 'LineWidth',1.2, 'MarkerSize',4);
xlim(ax2,[60 110]);
ylim(ax2,[0.002 0.012]);
grid(ax2,'on');
set(ax2,'FontSize',9);

figure; hold on;

plot(k_JL,  err_JL  - err_U3, '-o', 'LineWidth',1.5, 'MarkerSize',6, ...
    'DisplayName','JL - U_3');
plot(k_DCT, err_DCT - err_U3, '-^', 'LineWidth',1.5, 'MarkerSize',6, ...
    'DisplayName','DCT - U_3');
plot(k_DST, err_DST - err_U3, '-d', 'LineWidth',1.5, 'MarkerSize',6, ...
    'DisplayName','DST - U_3');

yline(0,'k--','LineWidth',1);

xlabel('Truncation k','FontSize',14);
ylabel('Difference in relative error','FontSize',14);
legend('Location','best','FontSize',14);
grid on;
set(gca,'FontSize',14);

%% 3) Relative‐error vs. truncation k for each method & ps
methods = {'JL'};      % match your Q names
ps_list = [10,50,100,220];      % same as above

fig = figure; hold on;
for m = 1:numel(methods)
    method = methods{m};
    T = readtable(fullfile(saveDir, method, 'results.csv'));
    for ps = ps_list
        idx = (T.s == ps);
        plot( T.k(idx), ...
              T.rel_err(idx), ...
              '-o', 'LineWidth',1.5, ...
              'DisplayName', sprintf('%s, s=%d', method, ps) );
    end
end
hold off;

xlabel('Truncation k');
ylabel('Relative error');
title('Error vs. Truncation k, all methods');
legend('Location','best','FontSize',14);
grid on;

exportgraphics(fig, fullfile(saveDir,'error_vs_truncation_k_all_methods.png'));

disp('Injective experiments complete.');
