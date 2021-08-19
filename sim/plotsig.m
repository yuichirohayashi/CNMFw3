function [newsig]=plotsig(sig,mag,order,normalize,show)

% �זE����Ca�g�`���v���b�g����
% INPUT
% sig: cellsort's cell_sig(ncells,t)
% mag: �g�嗦�@�傫������Ɣg�`�͍��������Ɋg��B�������d�Ȃ�Bdefault = 3 SD�̂Ƃ�0.2��
% order: �g�`���я��@snr/none/descend�@snr�ɂ����SNR�̍����� 
% normalize: no/height/sd 'no' �m�[�}���C�Y���Ȃ��A'height' �������낦��A'sd'SD������
% show: 0/1 0���ƕ\�����Ȃ��B
% OUTPUT
% newsig: �v���b�g�p�ɃX�P�[���A�Ԋu�𒲐߂����f�[�^ ncells x t

if ~exist('mag', 'var') || isempty(mag)
    mag = 3;
end
if ~exist('order', 'var') || isempty(order)
    order = 'none'
end
if ~exist('normalize', 'var') || isempty(order)
    normalize = 0;
end

% ���ёւ�
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

% �{��
ncells = size(sig,1); % total cell number in sig
range = max(max(sig))-min(min(sig));

sigsize = size(sig);
newsig = zeros(sigsize(1),sigsize(2));

for i = 1:ncells
    %geta = ncells - i;
    geta = -i;
    sig1 = sig(i,:);
    sdsig = std(sig1);
    % �m�[�}���C�Y
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