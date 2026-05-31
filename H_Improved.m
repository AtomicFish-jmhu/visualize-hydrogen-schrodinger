%% 初始化

close all;
clear all;

a = 1;

theta = 0:pi/60:pi;  % 定义域取样
phi = 0:pi/60:2*pi;
r = 0:0.1:30;

%% 球函数概率幅的可视化（Improved.）

% plotY_lm(3,2,theta,phi) 
plotY_l(5,theta,phi)

%% 径向几率分布可视化

% plotR_nl(3,1,a,r)
% plotR_n(3,a,r)

%% 完整氢原子定态波函数的概率幅可视化

% plotpsi_nlm(3,2,1,a,r,theta,phi)
% plotpsi_n(2,a,r,theta,phi)

%% 绘图函数定义

function plotY_l(l_max,theta,phi) % 绘制从l=0到l=l_max的所有球函数，m取非负值
figure()
for l = 0:l_max
for m = 0:l
    subplot(l_max+1,l_max+1,l*(l_max+1)+m+1)
    plotY_lm(l,m,theta,phi)
end
end
end

function plotY_lm(l,m,theta,phi) % 绘制球函数Y(l,m)，仅surf()，无figure(),subplot()等任何前后缀
    [Th,Ph] = meshgrid(theta,phi); % 犹记Th元素随列标变，Ph元素随行标变
    P = legendre(l,cos(theta),'norm'); % legendre(l,x)返回(l,m)=(l,0)至(l,m)=(l,l)阶连带勒让德多项式构成的矩阵，m为列坐标，x为输入值向量
     
    % 绘制复球函数（可选），以概率密度形式呈现
%     Y_lmMol = abs(P(m+1,:)).^2;
%     Y_lmMol = ones(length(phi),1)*Y_lmMol; % Y_lmMol与theta同型，将其制成与格点网对应的格式。还可用repmat(Y_lmMol,length(phi),1)
    
    % 绘制实球函数（可选）
    Y_lmMol = ones(length(phi),1)*P(m+1,:); % 先转成与格点网对应的格式。还可用repmat(P(m+1,:),length(phi),1)
    Y_lmMol = abs(Y_lmMol.*cos(m.*Ph)); % abs(Y_lmMol.*cos(m.*Ph)).^2 or abs(Y_lmMol.*cos(m.*Ph))
    
    x = Y_lmMol.*sin(Th).*cos(Ph); % 转换到直角坐标，用径长反映球函数的模
    y = Y_lmMol.*sin(Th).*sin(Ph);
    z = Y_lmMol.*cos(Th);
    
    % 用surf()绘图（可选）
    set(0,'defaultfigurecolor',[0.3 0.3 0.3]) % 背景色
    surf(x,y,z,'EdgeColor','none') % x,y,z矩阵同型，且相同位置上的一组(x(a,b),y(a,b),z(a,b))值由同一组theta - phi坐标生成
    title(sprintf('l=%d, m=%d',l,m),'color','w','FontSize',8)
    axis equal off 
    
    % 用scatter3()绘图（可选）
%     x = reshape(x,1,length(theta)*length(phi)); % scatter3只接受向量，故reshape
%     y = reshape(y,1,length(theta)*length(phi));
%     z = reshape(z,1,length(theta)*length(phi));
%     set(0,'defaultfigurecolor',[0.3 0.3 0.3]) % 背景色
%     h = scatter3(x,y,z,0.5); % 最后一个为记号大小
%     h.Marker = '*';
%     h.MarkerEdgeColor = 'g'; % 还有MarkerFaceColor、LineWidth（边缘宽度）属性
%     title(sprintf('l=%d, m=%d',l,m),'color','w')
%     axis equal off
end

function plotR_n(n_max,a,r) % 绘制从n=1到n=n_max的径向几率分布曲线，画在同一幅图中
    figure( )
    for n = 1:n_max
    for l = 0:(n-1)
        hold on;
        plotR_nl(n,l,a,r)
        title(sprintf('n<=%d',n_max)) % 这个title能够覆盖plotR_nl()里的title.
        str{n*(n-1)/2+l+1} = sprintf('n=%d, l=%d',n,l); % 字符型向量
        hold off;
    end
    end
    legend(str)
%     set(0,'defaultfigurecolor',[0.3 0.3 0.3])
%     set(gca, 'XColor', 'w', 'YColor', 'w')
end

function plotR_nl(n,l,a,r) % 绘制(n,l)态的径向几率分布曲线
    N = (2/(a^1.5*n^2*factorial(2*l+1)))*sqrt(factorial(n+l)/factorial(n-l-1));
    R_nl = N*exp(-r/(n*a)).*(2*r/(n*a)).^l.*hypergeom(l+1-n,2*l+2,2*r/(n*a));
    Rprob = r.^2.*R_nl.^2;
    plot(r,Rprob)
    title(sprintf('n=%d, l=%d',n,l))
end

function plotpsi_n(n_max,a,r,theta,phi) % 绘制氢原子定态n=1至n=n_max的可视化
    figure()
    index = 1;
    index_max = 0;
    for n = 1:n_max
    for l = 0:n-1
    for m = 0:l
        index_max = index_max+1;
    end
    end
    end
    for n = 1:n_max
    for l = 0:n-1
    for m = 0:l
        fprintf('Calculating n=%d, l=%d, m=%d, task %d/%d......\n',n,l,m,index,index_max);
        subplot(ceil(sqrt(index_max)),ceil(sqrt(index_max)),index)
        plotpsi_nlm(n,l,m,a,r,theta,phi)
        index = index+1;
    end
    end 
    end
end 

function plotpsi_nlm(n,l,m,a,r,theta,phi) % 绘制氢原子定态(n,l,m)的可视化，使用散点图，一点处的概率幅用记号颜色反映
    [R,Th,Ph] = meshgrid(r,theta,phi); % R,Th,Ph均为立方体矩阵，其中元素分别随第二位、第一位、第三位改变

    N_nl = (2/(a^1.5*n^2*factorial(2*l+1)))*sqrt(factorial(n+l)/factorial(n-l-1));
    R_nl = N_nl*exp(-R/(n*a)).*(2*R/(n*a)).^l.*hypergeom(l+1-n,2*l+2,2*R/(n*a));

    Y_lm = legendre(l,cos(theta),'norm');
    Y_lm = repmat(transpose(Y_lm(m+1,:)),1,length(r),length(phi));
%     Y_lm = Y_lm.*cos(m*Ph); % 绘制实球函数时插入

    psi2 = (abs(R_nl.*Y_lm)).^2;
    psi2 = psi2/max(psi2,[],'all');

    x = R.*sin(Th).*cos(Ph);
    y = R.*sin(Th).*sin(Ph);
    z = R.*cos(Th);

    N = length(r)*length(theta)*length(phi);
    x = reshape(x,1,N);
    y = reshape(y,1,N);
    z = reshape(z,1,N);
    psi2 = reshape(psi2,1,N);

    C = repmat(255*transpose(psi2),1,3); % 颜色映射
    S = 5*psi2+0.000000001; % 大小映射

    h = scatter3(x,y,z,S,C);
    h.Marker = '*';
    h.AlphaData = psi2; % 透明度映射
    h.MarkerEdgeAlpha = 'flat';
    title(sprintf('(%d,%d,%d)',n,l,m),'color','w','FontSize',14)
    set(0,'defaultfigurecolor',[0.3 0.3 0.3])
    axis equal off
    view(85,10)
end
