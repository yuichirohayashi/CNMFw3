function [A,C,b,f,center,Ain,Cin,bin,fin] = runplainCNMF(Y,K)

% plain CNMF‚ð3d movie‚É‘Î‚µŽÀs‚·‚é
% INPUT
% Y: 3d movie XYZT
% K: number of components to be found
% OUTPUT
% A:
% C:
% center: 

% zŽ²‚ð•Ï‚¦‚Ä‚Ý‚é
%Y = permute(Y, [1 3 2 4]);
%tau = [4,2,4];% default [3,3,2]->[4,4,2]‚š‚ð•Ï‚¦‚é‚Æ‚«

if ndims(Y) == 4
    [d1,d2,d3,T] = size(Y);                            % dimensions of dataset
else
    [d1,d2,T] = size(Y);
    d3 = 1;
end
d = d1*d2*d3;                                          % total number of pixels

%% Set parameters

%K = 200;                                          % number of components to be found
tau=[1,1,1];                            % std of gaussian kernel (size of neuron) 
p = 1;                  % order of autoregressive system (p = 0 no dynamics for slow imaging rate)
merge_thr = 0.95;                                 % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,'d3',d3,...                  % dimensions of datasets
    'search_method','dilate',...                 % search locations when updating spatial components
    'maxIter',15,...                             % number of NMF iterations during initialization default 15
    'deconv_method','constrained_foopsi',...     % activity deconvolution method
    'temporal_iter',2,...                        % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau,'nb',2 ...
    );
%% Data pre-processingS

[P,Y] = preprocess_data(Y,p);
%Cn = correlation_image_3D(Y); % for large datasets change with reshape(P.sn,d1,d2,d3), %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize
ff = find(sum(Ain)<1e-3*mean(sum(Ain)));   % remove very small components
Ain(:,ff) = [];
Cin(ff,:) = [];
center(ff,:) = [];

%% update spatial components
Yr = reshape(Y,d,T);
%clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);
