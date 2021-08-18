function [allscore] = simtest3d(method,fn)

% sim dataでcell sorting methodの評価をする
% INPUT
% method: cell sorting method, pcaica/plainCNMF/CNMFw3/CNMFE
% fn: file name to save
% OUTPUT
% allscore: temporal correlation (niter x nexp x ncells)

% simdata params
load('bg');
trucells = 100;
snr = [0.4];
scalef = 200;
niter = 10;
nexp = length(snr);

% cell sorting method params
nPCs = 200;
d1=100;
d2=100;
d3=10;
T = 1000;
% CNMFEは名前が衝突するので普段パスを通していない。使うときだけ通し、終わったら削除
if strcmp(method, 'CNMFE')
    pathCNMFE = genpath('E:\work\CNMF_E-master');
    addpath(pathCNMFE);
end

allscore = zeros(niter,nexp,trucells);
for j = 1:niter
    for i = 1:length(snr)
        [Y, fspk] = gensim3dmovRealBG(bg, snr(i),scalef);
        if strcmp(method,'pcaica')
            [ica_sig, ica_filters] = runcellsort3d(Y, nPCs);
            trace = ica_sig;
            filters = ica_filters;
        elseif strcmp(method, 'plainCNMF')
            [A,C] = runplainCNMF(Y,nPCs);
            trace = C;
            filters = permute(reshape(full(A),d1,d2,d3,size(A,2)),[4 1 2 3]);
        elseif strcmp(method, 'CNMFw3')
            [A,C] = runCNMFw3mat(Y, nPCs);
            trace = C;
            filters = permute(reshape(full(A),d1,d2,d3,size(A,2)),[4 1 2 3]);
        elseif strcmp(method, 'CNMFE')
            % これだけ2Dなのでデータを2次元にする
            slno = 5;
            Y = squeeze(Y(:,:,slno,:));
            Ysiz = size(Y)';
            save('Y2d', 'Y', 'Ysiz')
            try
                runCNMFEmat
            catch
                results.C = zeros(trucells,T);
                results.A = zeros(d1*d2,trucells)
                continue
            end
            trace = results.C;
            filters = reshape(full(results.A),d1,d2,[]);
            clearvars -except trace filters fspk allscore i j snr trucells bg method scalef fn niter pathCNMFE
        elseif strcmp(method, 'CNMFE10');
            
        end
        
        [r,p] = corr(fspk, trace');
        maxr = max(r,[],1);
        [maxr, ind] = sort(maxr,'descend');
        % traceがtrucells個に足りない場合がある。そのときはゼロ埋め
        if length(maxr) < trucells
            maxr = [maxr zeros(1,(trucells-length(maxr)))];
        else
            maxr = maxr(1:trucells);% 上位truecells個のみ取ってくる
        end
        maxr(isnan(maxr)) = 0; % nanつぶし
        allscore(j,i,:) = maxr;
        
        %         figure(i)
        %         subplot(1, 2, 1), plotsig(trace,0.5,'descend','sd',1),title(snr(i)),hold on
        %         if strcmp(method, 'CNMFE')
        %             subplot(1, 2, 2),imagesc(std(filters,0,3))
        %         else
        %             subplot(1, 2, 2), ploticafilters3d(filters,0.8,'smooth')
        %         end
        figure(10)
        subplot(length(snr), 1, i), histogram(maxr,[0:0.05:1]), title(snr(i)),hold on
        figure(11)
        cdfplot(maxr),hold on
        
        save(fn, 'allscore', 'snr')
    end
end
if strcmp(method, 'CNMFE')
    rmpath(pathCNMFE);
end




