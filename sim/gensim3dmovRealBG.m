function [Y, fspk, cellpos] = gensim3dmovRealBG(bg,r,scalef,fn);
% �J���V�E���C���[�W���O�V�~�����[�V�����f�[�^�̍쐬
% cell�͑ȉ~�̂ɂ���@Z�͊Ԋu���傫���̂ő��ΓI�ɔ����Ȃ�
% �����psf���|����
% cell�̓M�`�M�`�ɕ��׃����_����active cell��I��
% bg��real data�R���Ƃ���
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
% z step��20um�ʂ�z��
margin = [6 6 3];% �[�ɗ��Ȃ����������̂Ń}�[�W���󂯂� filter�|����Ƃ��������Ȃ�̂�
xyspace = 5; % distance between cells�@XY
zspace = 3;% distance between cells Z


% t=tau�̎���exp(-1)�ɂȂ��Ă���΂����B
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
% �v���}�C�P�͈̔͂ňʒu�ɗh�炬��^����܂�ŏ�������xyspace-2,zspace-2�ɂȂ�
randpos = [randi([-1,1],ncells,1), randi([-1,1],ncells,1), randi([-1,1],ncells,1)];
cellpos = cellpos + randpos;

% generate waveforms
spk = poissrnd(sfr/fs,[movlen,ncells]); % poisson distribution
fspk = filter(1,[1 -m],spk); % rise and decay 

for i=1:ncells
    Y(cellpos(i,1),cellpos(i,2),cellpos(i,3),:) = fspk(:,i);
end

% cell�^�ɂӂ���܂�
h = fspecial3('ellipsoid',[xyrad xyrad zrad]);
h = h/max(h(:));
Y = imfilter(Y,h); % replicate�ɂ���ƒ[�i���j�̒l���傫���Ȃ��Ă��܂�

% ���w�n��PSF���|����B
psf = fspecial3('gaussian',10,[2 2 2]);
Y = imfilter(Y, psf);

% �X�P�[�����O����signal��background�����Z
Y = Y/max(Y(:))*sheight + bg/max(bg(:))*bgheight;

% add shot noise
% imnoise��1e-12��1quantum�Ȃ̂ł��̉��{��signal��peak�ɂ��邩

Y = Y/max(Y(:))*1e-12*scalef;

Y = imnoise(Y, 'poisson');

Y = Y/max(Y(:))*p;

if nargin == 4
    save(fn, 'Y','fspk','cellpos')
end












