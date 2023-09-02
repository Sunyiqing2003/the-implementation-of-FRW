%导入导体1,2,主导体为1。
%导体和区域均为正方形，用横纵坐标在导体区域中的极小/极大值组成的四维向量进行表示。
B = [0, 1, 0, 1]; %边界
C1 = [0.2, 0.3, 0.2, 0.3]; %主导体
C2 = [0.7, 0.8, 0.7, 0.8]; %环境导体
G = [0.05, 0.45, 0.05, 0.45];%为主导体生成高斯面G，此处取为与主导体中心重合、边长为0.2的正方形
n = 0;    %采样次数为n
x = [];    %采样值构成的向量为x
C = 0;     %C是需要计算的导体1、2之间的电容
sigma = 1e-5;  %定义阈值,决定了采样的精度
Var = 1;  %Var表示采样点结果的方差
delta = 1e-3; %设置距离临界值
epi = 1e-7
while  (n < Var/sigma) && (n< 5000) %设置最大迭代次数                                                                                                                                                 %根据精度和计算结果的标准差来决定采样次数
    n = n+1; 
r = choose_cubic(G);%高斯面上取样

if (dist(r,B) < delta) || (dist(r, C1) < delta) ||(dist(r, C2) < delta) %边界小于临界值
   S = [r(1,1)-radius(r),r(1,1)+ radius(r),r(1,2)-radius(r),r(1,2)+radius(r) ];  %S是以r中心为半径正方体表面
   t = choose_cubic(S);
   w = weight_cubic(r,t);
   r = t;
else
   S = [r,radius(r)]; %S是以r中心为半径的球或者正方体
   t = choose_sphere(S);
   w = weight_sphere(r,t);
   r = t;
end

 while (onthesurface(r, C2) == 0) && (onthesurface(r, C2) == 0) && (onthesurface(r, B) == 0)%r不在导体j表面 
  if (dist(r,B) < delta) || (dist(r, C1) < delta)||(dist(r, C2) < delta)
   S=[r(1,1)-radius(r),r(1,1)+ radius(r),r(1,2)-radius(r),r(1,2)+radius(r) ]; %S是以r中心为半径正方体表面
   t = choose_cubic(S);
   w = weight_cubic(r,t);
   r = t;
  else
   S=[r,radius(r)]; %S是以r中心为半径的球
   t = choose_sphere(S);
   w = weight_sphere(r,t);
   r = t;
  end
 end
 x = [x,w];
 C = C + w;
Var = var(x);   %计算采样值的方差；
end
C = C/n;    %最后由采样的均值得到导体i和导体j之间的电容


function[r] = choose_cubic(G)
cho1 = randi([1,4]);
cho2 = unifrnd(-1,1,1,1);
cho3 = [0,0,0,0];
center = [(G(1,1)+G(1,2))/2, (G(1,3)+G(1,4))/2, (G(1,2)-G(1,1)+G(1,4)-G(1,3))/4];
%center 前两项表示G的中心的坐标，第三项表示G的半径
if cho1 == 1
    cho3 = [0,1] + [cho2,0];
end
if cho1 == 2
    cho3 = [0,-1] + [cho2,0];
end
if cho1 == 3
    cho3 = [1,0] + [0,cho2];
end
if cho1 == 4
    cho3 = [-1,0] + [0,cho2];
end
r = center(1,1:2) + cho3.*center(1,3)
end

function[r] = choose_sphere(G)
%函数表示在球体的表面随机选择一点。2D情况下为圆,此时输入的G为三维向量，前两项表示圆心，第三项为半径
cho1 = unifrnd(-pi,pi,1,1);
cho2 = [cos(cho1),sin(cho1)];
r = G(1,1:2)+G(1,3).*cho2;
end

function[t]=onthesurface(r,C2)
%判断r在不在导体2的表面上,误差值为1e-7
t = 0;
while (r(1,1) - C2(1) <1e-7 )||(r(1,1) - C2(2) >1e-7 )||(r(1,2)- C2(3) <1e-7)||(r(1,2) - C2(4) >1e-7)
    t = 0;
    break
end
if (abs(r(1,1)-C2(1))<1e-7)||(abs(r(1,1)-C2(2))<1e-7)||(abs(r(1,2)-C2(3))<1e-7)||(abs(r(1,2)-C2(4))<1e-7)
    t = 1;
end
end

function[d] = dist(r,B)
%计算采样点到导体或者边界面的距离，生成转移域
d = 0;
if (r(1) > B(1))&& (r(1) < B(2)) && (r(2) > B(3)) &&(r(2) < B(4))
    d = min([r(1)-B(1), B(2) - r(1), r() - B(3), B(4) - r(2)]);
end
if (r(1) < B(1))&& (r(2) > B(3)) &&(r(2) < B(4))
    d = B(1) - r(1);
end
if (r(1) > B(2))&& (r(2) > B(3)) &&(r(2) < B(4))
    d = r(1) - B(2);
end
if (r(1) > B(1))&& (r(1) < B(2))&& (r(2) < B(3))
    d = B(3) -r(2);
end
if (r(1) > B(1))&& (r(1) < B(2))&& (r(2) > B(4))
    d = r(2) - B(4);
end
if (r(1) < B(1)) && (r(2) < B(3))
    d = max([B(1) - r(1),B(3) -r(2)]);
end
if (r(1) < B(1)) && (r(2) > B(4))
    d = max([B(1) - r(1),r(2) - B(4)]);
end
if (r(1) > B(2)) && (r(2) < B(3))
    d = max([r(1) - B(2),B(3) -r(2)]);
end
if (r(1) > B(2)) && (r(2) > B(4))
    d = max([r(1) - B(2),r(2) - B(4)]);
end
end

function[L] = radius(r)
C1 = [0.2,0.3,0.2,0.3];
C2 = [0.7,0.8,0.7,0.8];
B = [0,1,0,1];
M = [dist(r,C1), dist(r, C2),dist(r,B)];
L=min(M);
%子函数dist，确定r到导体的距离
end

function[w] = weight_sphere(r,t)
%求解得到球面上的权重函数
%考虑2D球面上，由平均值定理，表面格林函数G(r,t)为1/2*pi*(r-t)
%经过计算，可得权重函数为 1/R,其中R为球面的半径。
w = 1 / norm(r - t)
end

function[w] = weight_cubic(r,t)
%求解得到立方体面上的权重函数
w = 0;
radius = max([abs(r(1) - t(1)), abs(r(2) - t(2))]);
x = [1/2, 1/2]+ (t-r)./(2*radius);
ymin = min([x(2),1/2]);
ymax = max([x(2),1/2]);
%对傅里叶级数作截断，取前100项
for n = 1:100
    w = w + 8* sin(n*pi/2)* sin(n*pi*x(1)) * sinh(n*pi*ymin) * sinh(n*pi*(1-ymax))/(n*sinh(n*pi));
end
end

