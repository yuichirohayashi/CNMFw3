function [Ain, Cin,  bin, fin, center, res] = greedyROI_w3(Y, K, options, sn)
%% a greedy method for detecting ROIs and initializing CNMF. in each iteration,
% moviewにspatial HPFをかけ、local corr*PNRを指標にROIを取っていく
%% CNMFEのgreedyROI_endoscopeを参考にしている
%% Input:
%   Y:  d1 x d2 x d3 X T matrx, imaging data
%   K:  scalar, maximum number of neurons to be detected.
%   options: struct data of paramters/options
%       gSig:   neuron diameter
%       zSig:   neuron diameter (Z axis) Z方向は光学系の都合で伸びているので別の値を入れる
%       gSiz:   maximum size of a neuron　ニューロンの探索範囲（箱）
%       nb:     number of background
%       min_corr: minimum threshold of correlation for segementing neurons
%   sn:     d X 1 vector, noise level of each pixel
%   debug_on: options for showing procedure of detecting neurons
%% Output:
%       Ain:  d X K' matrix, estimated spatial component
%       Cin:  K'X T matrix, estimated temporal component
%       bin:  d X nb matrix/vector, spatial components of the background
%       Cin:  nb X T matrix/vector, temporal components of the background
%       center: K' X 2, coordinate of each neuron's center
%       res:  d X T, residual after initializing Ain, Cin, bin, fin

% Written by Yuichiro Hayashi 2021 with modifications from  Pengchen Zhou (github.com/flatironinstitute/CaImAn-MATLAB/blob/master/utilities/greedyROI_corr.m)
% 

debug = false;

%% parameters
[d1,d2,d3,T] = size(Y);
if exist('sn', 'var')
    Y_std = sn;
else
    Y_std = std(Y, 0, ndims(Y));
end

Y_std = Y_std(:);
Y_median = median(Y, ndims(Y));
Y_median = reshape(Y_median, d1*d2*d3, []);
sigma = options.hpfsigma;
win = options.hpfwin;
gSig = options.gSig;
zSig = options.zSig;
gSiz = options.gSiz;
min_corr = options.min_corr;    % minimum local correaltion value to start one neuron
min_pnr = options.min_pnr;
lmaxwin = options.lmaxwin;
min_pixel = options.min_pixel;

sig = 5;    %%%%% これでどれだけ変わるか　HY_thrに使うthresholding noise by sig*std()

min_v_search = min_corr*min_pnr;
bd = options.bd; % 画像の縁の無視する部分pixel
nb = options.nb; % number of the background
%pSiz = 1;       % after selecting one pixel, take the mean of square box

%near the pixel as temporal activity. the box size is (2*pSiz+1)
if options.hpf == 1;
    psf = fspecial3('gaussian', win, sigma);
    ind_nonzero = (psf(:) >= max(psf(:,1)));
    psf = psf - mean(psf(ind_nonzero));% これでhpfになる
    psf(~ind_nonzero) = 0;
    HY = imfilter(Y,psf,'replicate');
else
    HY = Y; % 既にHPF掛けた画像をサーチするならHPF不要
end
HY = reshape(HY, d1*d2*d3, []);% ここでHYはmatrixに変形
Y = reshape(Y, d1*d2*d3, []);% ここでYはmatrixに変形
Yorig = Y;% Yはcell見つけるたびに引き算されていくので元を残しておく
% calc PNR
HY = bsxfun(@minus, HY, median(HY, 2));% medianを引く
HY_max = max(HY, [], 2);
Ysig = get_noise_fft(HY, options);% noise推定
PNR = reshape(HY_max./Ysig, d1, d2, d3);% max/noise=peak to noise ratio
PNR0 = PNR;
PNR(PNR<min_pnr) = 0;
% 小さすぎる（noiselevel=Ysigのsig倍以下）は0にする
HY_thr = HY;
HY_thr(bsxfun(@lt, HY_thr, Ysig*sig)) = 0;

% calc local corr
Cn = correlation_image_3D(HY_thr,6,[d1,d2,d3]); % 3dにおける6近傍でcorrを計算
Cn0 = Cn;   % backup
Cn(isnan(Cn)) = 0;
Cn = Cn + rand(size(Cn))*(1e-6);

% screen seeding pixels as center of the neuron
v_search = Cn.*PNR; % volume data
v_search0 = v_search;
v_search(or(Cn<min_corr, PNR<min_pnr)) = 0;

ind_search = false(d1*d2*d3,1);  % サーチ用フラグ　サーチしたら1を立てる
ind_search(v_search==0) = true; % ignore pixels with small correlations or low peak-noise-ratio

%% start init
Ain = zeros(d1*d2*d3, K);  % spatial components
Cin = zeros(K, T);      % temporal components
center = zeros(K, 3);   % center of the initialized components

% set boudndary to be 0
ind_bd = false(size(v_search));
ind_bd(1:bd,:,:) = true;
ind_bd(end-bd+1:end,:,:) = true;
ind_bd(:,1:bd,:) = true;
ind_bd(:,end-bd+1:end,:) = true;
ind_bd(:,:,1:bd) = true;
ind_bd(:,:,end-bd+1:end) = true;

searching_flag = 1;
k = 0;
while searching_flag
    % find local maxima
    % CNMFのgreedyROI_corrはPNRの高いところから順になめる
    % CNMFEはcorr*PNRの高いところからなめる
    %% なんで毎回medfiltが要るんだ？
    v_search = medfilt3(v_search, lmaxwin*[1 1 1])+randn(size(v_search))*(1e-10);% medfilt3は範囲奇数でないといけない
    v_search(ind_search) = 0;% 既にサーチしたところはゼロにしておく
    %tmp_d = 2*round(gSig)+1;% CNMFEでは、max3x3＝1個となりの範囲で最大値フィルタをかけている
    v_max = ordfilt3(v_search, 'max', lmaxwin); % 最大値フィルタ
    v_search(ind_bd) = 0;% 縁の無視する部分
    
    ind_search(v_search < min_v_search) = true;
    ind_localmax = find(and(v_search(:)==v_max(:), v_max(:)>0)); % 最大値フィルタの効く範囲でmaxと同じ＝localmaxとなる
    if(isempty(ind_localmax)); break; end
    %% localmaxを大きい順に調べていく
    [~, ind_sort] = sort(v_search(ind_localmax), 'descend');
    ind_localmax = ind_localmax(ind_sort);
    [r_peak, c_peak, z_peak] = ind2sub([d1,d2,d3], ind_localmax);
    for mcell = 1:length(ind_localmax);
        ind_p = ind_localmax(mcell);
        max_v = v_search(ind_p);
        ind_search(ind_p) = true;
        if max_v < min_v_search; continue; end % 足切り
        [r, c, z] = ind2sub([d1,d2,d3], ind_p);% ピーク座標
        %% ca traceを見て既処理or弱いものは飛ばす
        y0 = HY(ind_p, :);% ca trace
        y0_std = std(diff(y0));
        y0(y0<median(y0)) = 0;
        if (k>=1) && any(corr(Cin(1:k, :)', y0')>0.9); continue; end % 既に似たtraceを処理していれば飛ばす マジックナンバー注意
        if max(diff(y0)) < 3*y0_std; continue; end % too weak　弱すぎるやつ飛ばす　マジックナンバー注意
        %% ai,ciを推定するための箱 範囲+-gSiz bg推定必要なのでcellより大きい必要あり
        rsub = max(1, -gSiz+r):min(d1, gSiz+r);
        csub = max(1, -gSiz+c):min(d2, gSiz+c);
        zsub = max(1, -gSiz+z):min(d3, gSiz+z);
        [cind, rind, zind] = meshgrid(csub, rsub, zsub); % 順番注意
        [nr, nc, nz] =size(cind); % 箱の大きさ　端だと切れるので場所によって違う
        ind_nhood = sub2ind([d1,d2,d3], rind(:), cind(:), zind(:)); % neighbourhood
        HY_box = HY(ind_nhood, :);
        Y_box = Y(ind_nhood, :);
        ind_ctr = sub2ind([nr, nc, nz], r-rsub(1)+1, c-csub(1)+1, z-zsub(1)+1);% subscripts of the center本当？箱の最初の端では
        
        %% HYのudpate範囲psfでconvするのでaici推定範囲より2倍広くとってある
        rsub = max(1, -2*gSiz+r):min(d1, 2*gSiz+r);
        csub = max(1, -2*gSiz+c):min(d2, 2*gSiz+c);
        zsub = max(1, -2*gSiz+z):min(d3, 2*gSiz+z);
        [cind, rind, zind] = meshgrid(csub, rsub, zsub); % 順番注意
        ind_nhood_HY = sub2ind([d1,d2,d3], rind(:), cind(:), zind(:));% こっちの方はboxが大きい
        [nr2, nc2, nz2] =size(cind);
        
        % extract ai, ci
        sz = [nr, nc, nz];
        [ai, ci_raw, ind_success] = extract_ac3d(HY_box, Y_box, ind_ctr, sz); %箱の中でa,c推定
        if or(any(isnan(ai)), any(isnan(ci_raw))); ind_success=false; end
        if max(ci_raw/get_noise_fft(ci_raw)) < min_pnr; ind_success=false; end
        if sum(ai(:)>0) < min_pixel; ind_success=false; end
        if ind_success
            k = k+1;
            
            ci = ci_raw;
            Ain(ind_nhood, k) = ai;
            Cin(k, :) = ci_raw;
            Cin_raw(k, :) = ci_raw;
            center(k,:) = [r, c, z];
            % avoid buddy　近くに活動が似ているやつがいればそれは仲間なので探索から除外する
            ind_search(ind_nhood(ai>max(ai)*options.merge_thr)) = true;
            Y(ind_nhood, :) = Y_box - ai*ci; % ai,ciの推定値を生データから引く
            if options.hpf == 1;
                Hai = imfilter(reshape(Ain(ind_nhood_HY,k), nr2, nc2, nz2), psf, 'replicate');
            else
                Hai = reshape(Ain(ind_nhood_HY,k), nr2, nc2, nz2);
            end
            HY_box = HY(ind_nhood_HY, :) - Hai(:)*ci; %% HYの引き算
            HY(ind_nhood_HY,:) = HY_box;
            
            % update the maximum projection of HY
            Ysig_box = Ysig(ind_nhood_HY);
            temp = max(HY_box, [], 2);
            tmp_PNR = temp./Ysig_box;
            tmp_PNR(or(isnan(tmp_PNR), tmp_PNR<min_pnr)) = 0;
            PNR(ind_nhood_HY) = tmp_PNR;
            
            HY_box_thr = HY_box;  %thresholded version of HY
            HY_box_thr(bsxfun(@lt, HY_box, Ysig_box*sig)) = 0;
            
            % update correlation image ここあとで検討する
            tmp_Cn = correlation_image_3D(HY_box_thr, 6, [nr2, nc2, nz2]);% 3Dようにする
            tmp_Cn(or(isnan(tmp_Cn), tmp_Cn<min_corr)) = 0;
            Cn(ind_nhood_HY) = tmp_Cn;
            
            % update search value
            v_search(ind_nhood_HY) = Cn(ind_nhood_HY).*PNR(ind_nhood_HY);
            v_search(ind_bd) = 0;
            v_search(ind_search) = 0;
        else
            continue;
        end
        if k==K
            searching_flag = false;
            break;
        end
    end
end
center = center(1:k, :);
Ain = sparse(Ain(:, 1:k));
Cin = Cin(1:k, :);
Cin(Cin<0) = 0;
res = max(Yorig - Ain*Cin, 0); % 0以下を０にする

clear Y Yorig

fin = [mean(res);rand(nb-1,T)];

for nmfiter = 1:100
    bin = max((res*fin')/(fin*fin'),0);
    fin = max((bin'*bin)\(bin'*res),0);
end

if debug
    figure
    Cn0p = squeeze(max(Cn0,[],3));
    subplot(2,2,1),imagesc(Cn0p),title('Corr img')
    PNR0p = squeeze(max(PNR0,[],3));
    subplot(2,2,2),imagesc(PNR0p),title('PNR')
    v_search0p = squeeze(max(v_search0,[],3));
    subplot(2,2,3),imagesc(v_search0p),title('Corr x PNR')
    subplot(2,2,4),imagesc(v_search0p),title('init'),hold on
    plot(center(:,2),center(:,1),'k*')
end
%Cn = Cn0;
%PNR = PNR0;
