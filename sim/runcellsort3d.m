function [ica_sig, ica_filters] = runcellsort3d(Y, nPCs);

% cellsortbatch3d MukamelÇÃcellsortÇ3dópÇ…â¸ïœ

%load('simvbg');
%nPCs = 200;
nIC = nPCs;
PCuse = [1:nPCs];

% mu: spatial/temporalÇÃî‰ó¶ÅA0: pure spatial, 1: pure temporal
mu = 0.1;% 0.1->0.2Ç…ïœçX
termtol = 1e-5;
iter = 1000;

% perform 3D PCA 
[mixedsig,mixedfilters,CovEvals,covtrace,movm,movtm]=Cellsort3dPCAmat(Y,nPCs);

% perform 3D ICA
tic
[ica_sig,ica_filters,ica_A,numiter]=Cellsort3dICA(mixedsig,mixedfilters,CovEvals,PCuse,mu,nIC,[],termtol,iter);
toc

% subplot(2,2,1),plotsig(ica_sig,0.5,'none','sd',1);
% subplot(2,2,2),ploticafilters3d(ica_filters,0.8,'smooth')
% subplot(2,2,3),histogram