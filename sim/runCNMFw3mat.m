function [A,C,b,f,centerin,Ain,Cin,bin,fin,P,S,YrA,d1,d2,d3] = runCNMFw3mat(Y,K)

% CNMFw3を走らせる
% INPUT
% Y: input array (XYZT)
% K: num of components

% initmode: initialization settings [subtractBG, spatialHPF]
initmode = [1 1];

subtractBG = initmode(1);
preHPF = initmode(2);
[d1,d2,d3,T] = size(Y);                            % dimensions of dataset
d = d1*d2*d3;                                      % total number of pixels

%% Set parameters
% greedyの時はtau=1,1,1でないとうまくいかない。
% greedy_corr3dのときは実際のcell直径くらいにする[3 3 3][5 5 5]あたり
%

p = 1;%本番は1か2　framerate遅いときは１でいい    % order of autoregressive system (p = 0 no dynamics for slow imaging rate)
% merging threshold
% preprocess options
%subtractBG = 0; % estimate bg and subtract it from raw data
%preHPF = 1;     % apply hpf to raw data
% greedyROI_w3 options
% corr*pnrの計算をHPFかけてからやるか　既にHPF掛けた画像を使うときは0にする
if preHPF == 1;
    options.hpf = 0;
else
    options.hpf = 1;
end
options.init_method = 'greedy_corr3d';
options.nb = 0; % num of background
options.gSig = 3; %
options.zSig = 3; % neuron size
options.gSiz = 9; % maximum neuron size in XY axis
options.zSiz = 9; % maximum neuron size in Z axis
options.hpfsigma = [3 3 3]; % SD of gaussian HPF real dataだと222くらいか
options.hpfwin = 10; % gaussian HPF window size
options.min_corr = 0.8;
options.min_pnr = 8;
options.lmaxwin = 3; % lmaxwin^3の範囲でlocal maxをサーチ　奇数でないといけない
options.min_pixel = 9; % cellのpixel numberでの足切りライン
options.bd = 1; % 画像の縁の無視する部分
options.merge_thr = 0.85;
tau = options.gSig;

% update_spatial_components params
options.spatial_method = 'regularized';
options.min_size = 1;
options.max_size = 4;
options.dist = 3;
options.search_method = 'ellipse'
% update_temporal_components params
options.deconv_method = 'constrained_foopsi';
options.fr = 5; % frame rate
options.decay_time = 0.2;
options.lam_pr = 0.99; % false positive probability for determing lambda penalty
options.spk_SNR = 0.5;% spike SNR for min spike value
% Ain = double(full(Ain));
options.temporal_iter = 2;                        % number of block-coordinate descent steps
options.fudge_factor = 0.98;       % bias correction for AR coefficients
% general params
options.d1 = d1;
options.d2 = d2;
options.d3 = d3;

%% background estimation
if subtractBG ==1
    iter = 5;
    Y = reshape(Y, d, T);
    
    [bg_s, bg_t] = rank_1_factorization(Y,iter);
    Y = Y - (bg_s*bg_t);
    Y = reshape(Y, d1,d2,d3,T);
end

if preHPF ==1
    psf = fspecial3('gaussian', options.hpfwin, options.hpfsigma);
    ind_nonzero = (psf(:) >= max(psf(:,1)));
    psf = psf - mean(psf(ind_nonzero));% これでhpfになる
    psf(~ind_nonzero) = 0;
    Y = imfilter(Y,psf,'replicate');
end


%% pre-processing
[P,Y] = preprocess_data(Y,p);

%% initialization of spatial components
[Ain,Cin,bin,fin,centerin,res] = initialize_components_w3(Y,K,options,P);  % initialize
initdebug = 0;
if initdebug
    figure
    plotsig(Cin,0.5,'none','sd',1);
    ncells = size(Ain,2);
    figure
    ploticafilters3d(permute(reshape(full(Ain),d1,d2,d3,ncells),[4 1 2 3]),0.8,'smooth')
    hold on
    plot3(centerin(:,2),centerin(:,1),centerin(:,3),'k*')
    hold off
end

ff = find(sum(Ain)<1e-3*mean(sum(Ain)));   % remove very small components
Ain(:,ff) = [];
Cin(ff,:) = [];
centerin(ff,:) = [];

% Ain = double(full(Ain));
% Cin = double(Cin);
% bin = double(bin);
% fin = double(fin);
if ~isempty(Ain)
    %% update spatial components
    Yr = reshape(Y,d,T);
    %clear Y;
    [A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
    
    %% update temporal components
    P.p
    [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);
else
    A = Ain;
    C = Cin;
end

% figure
% subplot(3,2,1),plotsig(Cin,0.5,'none','sd',1);title Cin
% subplot(3,2,2),ploticafilters3d(permute(reshape(full(Ain),d1,d2,d3,size(Ain,2)),[4 1 2 3]),0.8,'smooth'); title Ain
% subplot(3,2,3),plotsig(C,0.5,'none','sd',1); title C
% subplot(3,2,4),ploticafilters3d(permute(reshape(full(A),d1,d2,d3,size(A,2)),[4 1 2 3]),0.8,'smooth'); title A
%
% [rin,p] = corr(fspk,Cin');
% maxrin = max(rin,[],1);
% [r,p] = corr(fspk,C');
% maxr = max(r,[],1);
%
% subplot(3,2,5),histogram(maxrin,[0:0.05:1])
% subplot(3,2,6),histogram(maxr,[0:0.05:1])
