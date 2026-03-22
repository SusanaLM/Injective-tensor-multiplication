
clear; clc;

%% load data

hcube = hypercube('indian_pines.dat');
A     = hcube.DataCube;

% normalize
A = A / fronorm(A);

saveDir  = './results_hyperspectral_hosvd/';
fileName = 'results';
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
if ~exist([saveDir,'orig'], 'dir'), mkdir([saveDir,'orig']); end

% get sizes
[n1,n2,n3]  = size(A);
nrmA        = fronorm(A);
storeA      = n1 * n2 * n3;

%% hosvd

if ~exist([saveDir,'hosvd/'], 'dir'), mkdir([saveDir,'hosvd/']); end
if ~exist([saveDir,'hosvd/img'], 'dir'), mkdir([saveDir,'hosvd/img']); end

kRange  = [1,10:10:140,size(A,1)];
results = zeros(length(kRange),length(kRange),size(A,3));
storage = zeros(size(results));

for k3 = 1:size(A,3)
    disp(['Starting k3 = ', num2str(k3),'...'])

    % full HOSVD
    [~,U] = hosvd3(A,size(A));
    
    count1 = 1;
    for k1 = kRange

        count2 = 1;
        for k2 = kRange
            % form core
            G = hosvdCore3(A,U,[k1,k2,k3]);
            
            % form approximation
            HApprox = hosvdApprox3(G,U,[k1,k2,k3]);
            
            % store results
            results(count1,count2,k3) = fronorm(A - HApprox);

            % storage costs
            storage(count1,count2,k3) = k1 * k2 * k3 + k1 * size(A,1) + k2 * size(A,2) + k3 * size(A,3);

            % update counter
            count2 = count2 + 1;
        end
        % update counter
        count1 = count1 + 1;
    end

    disp(['Finished k3 = ', num2str(k3),'.'])
end

varNames1 = cellfun(@(x) ['p',num2str(x),'_err'],num2cell(1:size(A,3)),'UniformOutput',false);
varNames2 = cellfun(@(x) ['p',num2str(x),'_store'],num2cell(1:size(A,3)),'UniformOutput',false);
varNames = cat(2,{'k1','k2'},varNames1,varNames2);

kk1 = kron(ones(1,length(kRange)),kRange);
kk2 = kron(kRange,ones(1,length(kRange)));

results = reshape(results / fronorm(A),[],size(results,3));
storage = reshape(numel(A) ./ storage,[],size(storage,3));
results = [kk1(:),kk2(:),results, storage];

T       = array2table(results,'VariableNames',varNames);

writetable(T,[saveDir,'hosvd/',fileName,'.csv'])


return;

%% setup storage

frame   = 75;
[R,G,B] = deal(26,16,8);

figure(1); clf;
imagesc(rescale(A(:,:,[R,G,B])));
cax = clim;
axis('off'); 
pbaspect([size(A,2) size(A,1) 1])
exportgraphics(gcf,[saveDir,'/img_RGB_.png'])

%% images for hosvd

paramList = {[2,5,n1,1],[2,5,9,9],[2,20,n1,1],[2,20,32,32],[100,5,n1,8],[100,5,36,36],[100,20,n1,38],[100,20,74,74]};

[~,U] = hosvd3(A,size(A));

results = [];
for count = 1:length(paramList)
    p = paramList{count}(1);
    k = paramList{count}(2);
    k1 = paramList{count}(3);
    k2 = paramList{count}(4);

    HCore       = hosvdCore3(A,U,[k1,k2,p]);
    HApprox     = hosvdApprox3(HCore,U,[k1,k2,p]);

    figure(1); clf;
    imagesc(rescale(HApprox(:,:,[R,G,B])));
    axis('image'); axis('off');
    exportgraphics(gcf,[saveDir,'hosvd/img/img_RGB_p_',num2str(p),'_k1_',num2str(k1),'_k2_',num2str(k2),'.png']);
    

    rel_err = fronorm(A - HApprox) / fronorm(A);
    rel_store = numel(A) / (k1 * k2 * p + k1 * n1 + k2 * n2 + p * n3);

    results = cat(1,results,[p,k1,k2,rel_err,rel_store]);
end

T = array2table(results,'VariableNames',{'p','k1','k2','rel_err','rel_store'});
writetable(T,[saveDir,'hosvd/',fileName,'_small.csv'])

