% Based on the code by Elizabeth Newman, available at
% https://github.com/elizabethnewman/projected-products

% Edited by: Susana Lopez Moreno, Pusan National University, 2026

clear; clc;

%% Reproducibility seed
rng(2025);

%% Load data
video = read(VideoReader('traffic.mj2'));
A = zeros(size(video,1), size(video,2), size(video,4));

for i = 1:size(video,4)
    A(:,:,i) = squeeze(im2gray(video(:,:,:,i)));
end

% normalize
A = A / fronorm(A);

saveDir  = './results_traffic_injective/';
fileName = 'results';
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
if ~exist([saveDir,'orig'], 'dir'), mkdir([saveDir,'orig']); end

% sizes
[m,n,p] = size(A);
nrmA    = fronorm(A);
storeA  = m * n * p;

%% Mode-3 SVD for U3
[Z,S,~] = svd(modeUnfold(A,3),'econ');

T = array2table([(1:size(S,1))', diag(S)], ...
    'VariableNames', {'index','sigma'});
writetable(T, [saveDir,'mode3_singular_values.csv']);

%% Injective matrices

% JL
oversamp = 2;
target   = oversamp * p;
N        = 2^nextpow2(target);
D        = diag(sign(randn(N,1)));
Sfull    = (1/sqrt(p)) * hadamard(N) * D;
rows     = randperm(N, p);
M_JL     = Sfull(rows, :)';              % N x p
M_JL_inv = (M_JL' * M_JL) \ (M_JL');

% DCT
q_dct     = 2*p;
Cfull     = dctmtx(q_dct);
M_DCT     = Cfull(:,1:p);                % q x p
M_DCT_inv = (M_DCT' * M_DCT) \ (M_DCT');

% DST
Sdst = zeros(q_dct, q_dct);
for a = 1:q_dct
    for b = 1:q_dct
        Sdst(a,b) = sqrt(2/(q_dct+1)) * sin(pi*a*b/(q_dct+1));
    end
end
M_DST     = Sdst(:,1:p);                 % q x p
M_DST_inv = (M_DST' * M_DST) \ (M_DST');


Q = {'JL',  M_JL, ...
     'U3',  Z', ...
     'DCT', M_DCT, ...
     'DST', M_DST};

Q_inv = {'JL_inv',  M_JL_inv, ...
         'U3_inv',  Z, ...
         'DCT_inv', M_DCT_inv, ...
         'DST_inv', M_DST_inv};

headers = {'k','s','rel_err','storage','comp_ratio'};

%% *_M-SVD 
for i = 1:2:length(Q)

    methodName = Q{i};
    M          = Q{i+1};
    Minv       = Q_inv{i+1};

    if ~exist([saveDir,methodName], 'dir'), mkdir([saveDir,'/',methodName]); end
    if ~exist([saveDir,methodName,'/img'], 'dir'), mkdir([saveDir,'/',methodName,'/img']); end

    % Push-forward
    AHat = modeProduct(A, M);

    % Facewise SVD
    [UHat, SHat, VHat] = facewiseSVD(AHat);

    results = [];
    q = size(AHat,3);

    for k = 1:min(m,n)

        % rank-k truncation in transformed domain
        AkHat = facewise(UHat(:,1:k,:), ...
                 facewise(SHat(1:k,1:k,:), tran(VHat(:,1:k,:))));

        % Pull back full approximation
        Ak = modeProduct(AkHat, Minv);

        for s = 1:p

            % Slice in original domain
            approx = Ak(:,:,1:s);

            % Relative error on first s slices
            err = fronorm(A(:,:,1:s) - approx) / fronorm(A(:,:,1:s));

            % Storage
            comp = m * k * s + n * k * s;

            if ismember(methodName, {'JL','DCT','DST'})
                comp = comp + q * p;
            elseif strcmp(methodName, 'I')
                comp = comp;
            else
                comp = comp + p * p;
            end

            results = cat(1, results, [k, s, err, comp, storeA/comp]);
        end

        fprintf('%s: completed k = %d\n', methodName, k);
    end

    T = array2table(results, 'VariableNames', headers);
    writetable(T, [saveDir,'/',methodName,'/',fileName,'.csv']);
end

%% ------------------------------------------------------------------------
% Create relative-error heatmaps and snapshot images
%% ------------------------------------------------------------------------

% Preload tensor
K = min(size(A,1), size(A,2));
P = size(A,3);

methodNames = {'JL','U3','DCT','DST'};

REScales = zeros(2, numel(methodNames));

for count = 1:numel(methodNames)
    name = methodNames{count};

    T  = readtable([saveDir,'/',name,'/',fileName,'.csv']);

    % Columns are: {'k','s','rel_err','storage','comp_ratio'}
    % Need rel_err, reshaped as (#s) x (#k)
    RE = reshape(T.rel_err, P, []);

    fig = figure(1);
    set(fig, 'MenuBar', 'none', 'ToolBar', 'none', 'Color', 'w');
    clf(fig);
    imagesc(RE);
    colormap hot;
    axis square
    axis off;
    pbaspect([K P 1])
    exportgraphics(gcf, [saveDir,'/',name,'/RE.png'])

    % store color bounds
    REScales(1,count) = min(RE(:));
    REScales(2,count) = max(RE(:));
end

Tscale = array2table(REScales, ...
    'VariableNames', methodNames, ...
    'RowNames', {'min','max'});
writetable(Tscale, [saveDir,'rel_error_color_scales.csv'], ...
    'WriteRowNames', true);

%% ------------------------------------------------------------------------
% Store original frame and reconstructed frames
%% ------------------------------------------------------------------------

% Choose (k,s) pairs for snapshots
ksPairs = [1,10; 5,10; 20,50; 20,100];
frame   = min(75, size(A,3));   

clf(fig);
imagesc(A(:,:,frame));
cax = clim;
colormap gray
axis off
pbaspect([size(A,2) size(A,1) 1])
exportgraphics(gcf, [saveDir,'/frame_',num2str(frame),'.png'])

for j = 1:2:length(Q)
    methodName = Q{j};
    M          = Q{j+1};
    Minv       = Q_inv{j+1};

    % Transform data
    AHat = modeProduct(A, M);

    % Store transformed-domain feature slices
    featureIdx = unique(round(linspace(1, size(AHat,3), min(6,size(AHat,3)))));
    Tfeat = [];

    for ii = featureIdx
        clf(fig);
        imagesc(AHat(:,:,ii));
        colormap parula
        Tfeat = cat(1, Tfeat, ...
            [ii, min(AHat(:,:,ii),[],'all'), max(AHat(:,:,ii),[],'all')]);
        axis off
        pbaspect([size(A,2) size(A,1) 1])
        exportgraphics(gcf, [saveDir,'/',methodName,'/img/feature_',num2str(ii),'.png'])
    end

    Tfeat = array2table(Tfeat, ...
        'VariableNames', {'slice','clim_min','clim_max'});
    writetable(Tfeat, [saveDir,'/',methodName,'/feature_clim.csv'])

    % Compute facewise SVD
    [UHat, SHat, VHat] = facewiseSVD(AHat);

    % Store reconstructed snapshots for selected (k,s)
    for i = 1:size(ksPairs,1)
        kk = ksPairs(i,1);
        ss = ksPairs(i,2);

        % Rank-k approximation in transformed domain
        AkHat = facewise(UHat(:,1:kk,:), ...
                 facewise(SHat(1:kk,1:kk,:), tran(VHat(:,1:kk,:))));

        % Pull back full approximation
        Ak_full = modeProduct(AkHat, Minv);

        % Keep only first s slices
        Ak = Ak_full(:,:,1:ss);

        % Choose a frame that exists in the truncated tensor
        frame_use = min(frame, ss);

        clf(fig);
        imagesc(Ak(:,:,frame_use));
        colormap gray
        clim(cax)
        axis off
        pbaspect([size(A,2) size(A,1) 1])
        exportgraphics(gcf, ...
            [saveDir,'/',methodName,'/img/frame_',num2str(frame_use), ...
             '_k',num2str(kk),'_s',num2str(ss),'.png'])
    end
end
disp('Traffic experiments complete.');