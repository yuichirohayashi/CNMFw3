function ploticafilters3d(icafilters,thr,smooth)

% cellsort�ō��ꂽ3D spatial filter��}��
% INPUT
% icafilters: (N x X x Y xZ)
% thr: ���l�ʂ̒l max��thr�{
% smooth: filter���X���[�W���O���邩'smooth' or 'none'�A�f�t�H���g none

sigma = 0.4; %�@�K�E�V�A���t�B���^�̃V�O�}

if nargin == 2 % �f�t�H�̓X���[�W���O�Ȃ�
    smooth = 'none';
end

clmap = hsv;
for i=1:size(icafilters,1)
    data = squeeze(icafilters(i,:,:,:));
    if strcmp(smooth,'smooth')
        data = smooth3(data,'gaussian',3,sigma);
        %data = smooth3(data);
    end
    % patch: �h��Ԃ��ꂽ���p�`
    % isocaps: ���l��
    % isonormals: ���l�ʂ̒��_�̖@��
    % view ���_ view(3) �f�t�H
    patch(isocaps(data,thr*max(data(:))),'FaceColor',clmap(rem(i*11,64)+1,:),'EdgeColor','none');
    p1 = patch(isosurface(data,thr*max(data(:))),'FaceColor',clmap(rem(i*11,64)+1,:),'EdgeColor','none');
    alpha(0.5);
    isonormals(data,p1)
    view(3);
    axis vis3d %tight
    set(gca,'XLim',[1,size(icafilters,3)],'YLim',[1,size(icafilters,2)],'ZLim',[1,size(icafilters,4)]);% �}��yxz�ɂȂ��Ă���̂�
    
    %colormap jet
    lighting gouraud
    title(i); drawnow;
    hold on
end
hold off
daspect([1 1 0.1]);% z����10�{�ɐL�΂��B��̂������񂶂̃A�X�y�N�g��ɂȂ�
for i=1:4% ���邭����
    %camlight left;
end