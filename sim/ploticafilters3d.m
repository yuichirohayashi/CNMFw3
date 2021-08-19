function ploticafilters3d(icafilters,thr,smooth)

% cellsortで作られた3D spatial filterを図示
% INPUT
% icafilters: (N x X x Y xZ)
% thr: 等値面の値 maxのthr倍
% smooth: filterをスムージングするか'smooth' or 'none'、デフォルト none

sigma = 0.4; %　ガウシアンフィルタのシグマ

if nargin == 2 % デフォはスムージングなし
    smooth = 'none';
end

clmap = hsv;
for i=1:size(icafilters,1)
    data = squeeze(icafilters(i,:,:,:));
    if strcmp(smooth,'smooth')
        data = smooth3(data,'gaussian',3,sigma);
        %data = smooth3(data);
    end
    % patch: 塗りつぶされた多角形
    % isocaps: 等値面
    % isonormals: 等値面の頂点の法線
    % view 視点 view(3) デフォ
    patch(isocaps(data,thr*max(data(:))),'FaceColor',clmap(rem(i*11,64)+1,:),'EdgeColor','none');
    p1 = patch(isosurface(data,thr*max(data(:))),'FaceColor',clmap(rem(i*11,64)+1,:),'EdgeColor','none');
    alpha(0.5);
    isonormals(data,p1)
    view(3);
    axis vis3d %tight
    set(gca,'XLim',[1,size(icafilters,3)],'YLim',[1,size(icafilters,2)],'ZLim',[1,size(icafilters,4)]);% 図がyxzになっているので
    
    %colormap jet
    lighting gouraud
    title(i); drawnow;
    hold on
end
hold off
daspect([1 1 0.1]);% z軸を10倍に伸ばす。大体いいかんじのアスペクト比になる
for i=1:4% 明るくする
    %camlight left;
end