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
[n1,n2,n3]  = size(A);
nrmA        = fronorm(A);
storeA      = n1 * n2 * n3;

%% Form matrices

% Injective matrix JL

oversamp = 2;
target  = oversamp * n3;          % originally 440
N       = 2^nextpow2(target);    % rounds up to 512
D       = diag(sign(randn(N,1)));
Sfull   = (1/sqrt(n3))*hadamard(N)*D;   % N×N
rows    = randperm(N, n3);
M       = Sfull(rows, :)';            % now N×p
q       = N;                          % use this N for cost
Mp = (M.'*M)\(M.');

% Surjective/invertible matrices for comparison: U_3 and I
[Z,S,~] = svd(modeUnfold(A,3),'econ');

I = eye(n3);

% List of matrices and names
Q = {'JL', M, 'U3', Z', 'I', I};
Q_inv = {'JL_inv', Mp, 'U3_inv', Z, 'I_inv', I'};

% Headers
headers = {'k','p','rel_err','comp_ratio','rel_comp_ratio'};

%% *_M-SVD

for i = 1:2:length(Q)
    
    if ~exist([saveDir,Q{i}], 'dir'), mkdir([saveDir,'/',Q{i}]); end
    if ~exist([saveDir,Q{i},'/img'], 'dir'), mkdir([saveDir,'/',Q{i},'/img']); end
    
    % Push‐forward & facewise SVD
    AHat             = modeProduct(A,Q{i+1});
    [UHat,SHat,VHat] = facewiseSVD(AHat);

    results = [];
    
    % Truncation loop
    for k = 1:min(size(A,1:2))
        % Facewise multiplication of Uhat,Shat,Vhat
        AkHat = facewise(UHat(:,1:k,:),facewise(SHat(1:k,1:k,:), tran(VHat(:,1:k,:))));

        % Pull-back
        Ak  = modeProduct(AkHat,Q_inv{i+1});
        for p = 1:size(A,3)

            % Slice the approximation
            approx = Ak(:,:,1:p);

            % Compute error and storage
            err  = fronorm(A(:,:,1:p) - approx);
            comp = n1 * k * p + n2 * k * p;

            if strcmp(Q{i},'JL')
                comp = comp + q * p;
            else
                comp = comp + n3 * p;
            end
            
            % Store results
            results = cat(1,results,[k,p,err, comp,storeA / comp]);
 
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
[n1,n2,n3] = size(A);

% Injective matrix JL

oversamp = 2;
target  = oversamp * n3;          % originally 440
N       = 2^nextpow2(target);    % rounds up to 512
D       = diag(sign(randn(N,1)));
Sfull   = (1/sqrt(n3))*hadamard(N)*D;   % N×N
rows    = randperm(N, n3);
M       = Sfull(rows, :)';            % now N×p
q       = N;                          % use this N for cost
Mp = (M.'*M)\(M.');

% Surjective/invertible matrices for comparison: U_3 and I
[Z,S,~] = svd(modeUnfold(A,3),'econ');

I = eye(n3);

% List of matrices and names
Q = {'JL', M, 'U3', Z', 'I', I};
Q_inv = {'JL_inv', Mp, 'U3_inv', Z, 'I_inv', I'};


%% 1) RGB snapshots
for i = 1:2:length(Q)
    outImg  = fullfile(saveDir,'/',Q{i},'/img');

    if ~exist(outImg,'dir'), mkdir(outImg); end

    Ah = modeProduct(A,Q{i+1});
    [Uhat,Shat,Vhat] = facewiseSVD(Ah);

    % original RGB
    fig = figure('Visible','off');
    imagesc(rescale(A(:,:, [R G B])));
    axis off; pbaspect([n2 n1 1]);
    exportgraphics(fig, fullfile(outImg,'orig_RGB.png'));
    close(fig);

    % approximations
    for idx = 1:size(kpPairs,1)
      k  = kpPairs(idx,1);
      ps = kpPairs(idx,2);

      % reconstruct rank‐k in transform domain, pull back all n3 bands
      Ahat_k     = facewise( ...
                      Uhat(:,1:k,:), ...
                      facewise(Shat(1:k,1:k,:), tran(Vhat(:,1:k,:))) ...
                   );

      approx_full= modeProduct(Ahat_k,Q_inv{i+1});        % n1×n2×n3
      % rgb         = rescale(approx_full(:,:, [R G B]));

      % 2) pad out to n3 bands
      pad_cube = zeros(n1, n2, n3);
      pad_cube(:,:,1:size(approx_full,3)) = approx_full;

      % 3) now safely index any R,G,B (bands > p are zero)
      rgb = rescale(pad_cube(:,:, [R G B]));

      fig = figure('Visible','off');
      imagesc(rgb);
      axis off; pbaspect([n2 n1 1]);
      fname = sprintf('img_RGB_k%02d_p%03d.png',k,ps);
      exportgraphics(fig, fullfile(outImg,fname));
      close(fig);
    end
end


%% 2) Plot relative error vs k for all JL and U3

% 1) read in your CSV results
projectDir = pwd;   % or hard-code it
cd(projectDir);
T_JL  = readtable('./results_hyperspectral_injective/JL/results.csv');
T_U3  = readtable('./results_hyperspectral_injective/U3/results.csv');

% 2) fix the slice‐counts
p_JL  = 220;  
p_U3 = 220;   

% 3) select rows
idx_JL = T_JL.p == p_JL;
idx_U3 = T_U3.p  == p_U3; 

% 4) extract k and relative error
k_JL   = T_JL.k(idx_JL);
err_JL = T_JL.rel_err(idx_JL);

k_U3   = T_U3.k(idx_U3);
err_U3 = T_U3.rel_err(idx_U3);

% 5) plot
figure; hold on;
plot(k_JL, err_JL, '-o', 'LineWidth',1.5, 'DisplayName','Injective (JL), p=220');
plot(k_U3, err_U3, '-s', 'LineWidth',1.5, 'DisplayName','Surjective (U_3), p=220');
hold off;

xlabel('Truncation k');
ylabel('Relative error');
title('Relative error vs. k  (fixed p = 220)');
legend('Location','best','FontSize',14);
grid on;

%% 3) Relative‐error vs. truncation k for each method & ps
methods = {'JL'};      % match your Q names
ps_list = [10,50,100,220];      % same as above

fig = figure; hold on;
for m = 1:numel(methods)
    method = methods{m};
    T = readtable(fullfile(saveDir, method, 'results.csv'));
    for ps = ps_list
        idx = (T.p == ps);
        plot( T.k(idx), ...
              T.rel_err(idx), ...
              '-o', 'LineWidth',1.5, ...
              'DisplayName', sprintf('%s, p=%d', method, ps) );
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