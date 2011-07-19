

[nCG,t,r] = xlsread('Basepair_Conservation_23S_CG.xls');
[nGC,t,r] = xlsread('Basepair_Conservation_23S_GC.xls');
[nAU,t,r] = xlsread('Basepair_Conservation_23S_AU.xls');
[nUA,t,r] = xlsread('Basepair_Conservation_23S_UA.xls');

AU = [1 0 0];
UA = [-0.5 sqrt(3)/2 0];
CG = [-0.5 -sqrt(3)/2 0];
GC = [0 0 sqrt(2)];

P = [AU; UA; CG; GC];

figure(1)
clf

hold on
text(AU(1),AU(2),AU(3),'AU');
text(UA(1),UA(2),UA(3),'UA');
text(CG(1),CG(2),CG(3),'CG');
text(GC(1),GC(2),GC(3),'GC');

plot3(P([1 2],1),P([1 2],2),P([1 2],3));
plot3(P([1 3],1),P([1 3],2),P([1 3],3));
plot3(P([1 4],1),P([1 4],2),P([1 4],3));
plot3(P([2 3],1),P([2 3],2),P([2 3],3));
plot3(P([2 4],1),P([2 4],2),P([2 4],3));
plot3(P([3 4],1),P([3 4],2),P([3 4],3));

for i = 1:length(nCG(:,1)),
  T = (nCG(i,4)*AU + nCG(i,7)*CG + nCG(i,10)*GC + nCG(i,13)*UA)/sum(nCG(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'y.');
end


for i = 1:length(nGC(:,1)),
  T = (nGC(i,4)*AU + nGC(i,7)*CG + nGC(i,10)*GC + nGC(i,13)*UA)/sum(nGC(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'g.');
end


for i = 1:length(nAU(:,1)),
  T = (nAU(i,4)*AU + nAU(i,7)*CG + nAU(i,10)*GC + nAU(i,13)*UA)/sum(nAU(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'r.');
end


for i = 1:length(nUA(:,1)),
  T = (nUA(i,4)*AU + nUA(i,7)*CG + nUA(i,10)*GC + nUA(i,13)*UA)/sum(nUA(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'b.');
end

T = (0.11*AU + 0.25*CG + 0.25*GC + 0.11*UA)/(0.11*2+0.25*2);
plot3(T(1),T(2),T(3),'k*');

rotate3d on
axis vis3d
axis off

% ------------------------------------ make a planar plot

figure(2)
clf

AU = [0 0 0];
UA = [1 0 0];
CG = [1 1 0];
GC = [0 1 0];

hold on
text(AU(1)-0.05,AU(2)-0.05,AU(3),'AU');
text(UA(1)     ,UA(2)-0.05,UA(3),'UA');
text(CG(1)     ,CG(2)+0.05,CG(3),'CG');
text(GC(1)-0.05,GC(2)+0.05,GC(3),'GC');

for i = 1:length(nCG(:,1)),
  T = (nCG(i,4)*AU + nCG(i,7)*CG + nCG(i,10)*GC + nCG(i,13)*UA)/sum(nCG(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'y.');
end


for i = 1:length(nGC(:,1)),
  T = (nGC(i,4)*AU + nGC(i,7)*CG + nGC(i,10)*GC + nGC(i,13)*UA)/sum(nGC(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'g.');
end


for i = 1:length(nAU(:,1)),
  T = (nAU(i,4)*AU + nAU(i,7)*CG + nAU(i,10)*GC + nAU(i,13)*UA)/sum(nAU(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'r.');
end


for i = 1:length(nUA(:,1)),
  T = (nUA(i,4)*AU + nUA(i,7)*CG + nUA(i,10)*GC + nUA(i,13)*UA)/sum(nUA(i,[4 7 10 13]));
  plot3(T(1),T(2),T(3),'b.');
end

T = (0.11*AU + 0.25*CG + 0.25*GC + 0.11*UA)/(0.11*2+0.25*2);
plot3(T(1),T(2),T(3),'k*');

rotate3d on
axis([-0.1 1.1 -0.1 1.1 0 1]);
axis off
view(2)
