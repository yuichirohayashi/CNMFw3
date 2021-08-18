function [Y, fspk, cellpos] = gensim3dmovRealBG(bg,r,scalef,fn);
% カルシウムイメージングシミュレーションデータの作成
% cellは楕円体にする　Zは間隔が大きいので相対的に薄くなる
% それにpsfを掛ける
% cellはギチギチに並べランダムにactive cellを選択
% bgをreal data由来とする
% INPUT
% bg: background data, XYZT
% r: signal to bg ratio
% scalef: #photons
% fn: filename to save
% OUTPUT
% Y: movie data, XYZT
% fspk: calcium traces, t x ncells
% cellpos: cell position, ncells x 3

%%%%% signal parms %%%%%
%r = 0.5; % signal/bg ratio: sPeak/bgPeak
p = 200; % peak height: sPeak+bgPeak
fs = 5; % sampling rate in Hz
tau = 0.5; % decay time constant;
sfr = 0.1; % mean firing rate in Hz

%%%%% image params %%%%%
fovsize = [100, 100, 10]; % xyz
movlen = 1000; % movie length
% set cell pos
ncells = 100; % cell number
xyrad = 2; % radius (xy)
zrad = 1; % radius (z)
% z stepは20um位を想定
margin = [6 6 3];% 端に来ない方がいいのでマージン空ける filter掛けるとおかしくなるので
xyspace = 5; % distance between cells　XY
zspace = 3;% distance between cells Z


% t=tauの時にexp(-1)になっていればいい。
m = exp(-1/(tau*fs));
sheight = r*p/(1+r)
bgheight = p/(1+r)

% make white 3d movie
Y = zeros([fovsize movlen]);

[posX,posY,posZ] = meshgrid(margin(1):xyspace:(fovsize(1)-margin(1)), margin(2):xyspace:(fovsize(2)-margin(2)), margin(3):zspace:fovsize(3));
cellpick = randperm(length(posX(:)),ncells);
posX = posX(cellpick)';
posY = posY(cellpick)';
posZ = posZ(cellpick)';
cellpos= [posX, posY, posZ];
% プラマイ１の範囲で位置に揺らぎを与えるつまり最小距離はxyspace-2,zspace-2になる
randpos = [randi([-1,1],ncells,1), randi([-1,1],ncells,1), randi([-1,1],ncells,1)];
cellpos = cellpos + randpos;

% generate waveforms
spk = poissrnd(sfr/fs,[movlen,ncells]); % poisson distribution
fspk = filter(1,[1 -m],spk); % rise and decay 

for i=1:ncells
    Y(cellpos(i,1),cellpos(i,2),cellpos(i,3),:) = fspk(:,i);
end

% cell型にふくらます
h = fspecial3('ellipsoid',[xyrad xyrad zrad]);
h = h/max(h(:));
Y = imfilter(Y,h); % replicateにすると端（下）の値が大きくなってしまう

% 光学系のPSFを掛ける。
psf = fspecial3('gaussian',10,[2 2 2]);
Y = imfilter(Y, psf);

% スケーリングしてsignalとbackground足し算
Y = Y/max(Y(:))*sheight + bg/max(bg(:))*bgheight;

% add shot noise
% imnoiseは1e-12が1quantumなのでその何倍をsignalのpeakにするか

Y = Y/max(Y(:))*1e-12*scalef;

Y = imnoise(Y, 'poisson');

Y = Y/max(Y(:))*p;

if nargin == 4
    save(fn, 'Y','fspk','cellpos')
end












