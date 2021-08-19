function [newsig]=plotsig(sig,mag,order,normalize,show)

% 細胞毎のCa波形をプロットする
% INPUT
% sig: cellsort's cell_sig(ncells,t)
% mag: 拡大率　大きくすると波形は高さ方向に拡大。しかし重なる。default = 3 SDのとき0.2位
% order: 波形並び順　snr/none/descend　snrにするとSNRの高い順 
% normalize: no/height/sd 'no' ノーマライズしない、'height' 高さそろえる、'sd'SD揃える
% show: 0/1 0だと表示しない。
% OUTPUT
% newsig: プロット用にスケール、間隔を調節したデータ ncells x t

if ~exist('mag', 'var') || isempty(mag)
    mag = 3;
end
if ~exist('order', 'var') || isempty(order)
    order = 'none'
end
if ~exist('normalize', 'var') || isempty(order)
    normalize = 0;
end

% 並び替え
if strcmp(order,'snr')
    sd = rms(sig');
    sd = sd';
    peak = max(sig,[],2);
    snr = peak./sd;
    [B,I] = sort(snr,'ascend');
    sig = sig(I,:);
elseif strcmp(order,'descend')
    sig = sig(flipud([1:size(sig,1)]'),:);
end

% 倍率
ncells = size(sig,1); % total cell number in sig
range = max(max(sig))-min(min(sig));

sigsize = size(sig);
newsig = zeros(sigsize(1),sigsize(2));

for i = 1:ncells
    %geta = ncells - i;
    geta = -i;
    sig1 = sig(i,:);
    sdsig = std(sig1);
    % ノーマライズ
    if strcmp(normalize,'height')
        normsig1 = mag*(sig1-min(sig1))/(max(sig1)-min(sig1));
    elseif strcmp(normalize,'sd')
        normsig1 = mag*sig1/sdsig;
    else
        normsig1 = mag*sig1;
    end
    if show == 1
        plot(normsig1 + geta)
    end
    newsig(i,:) = normsig1 + geta;
    hold on, axis tight
end
hold off