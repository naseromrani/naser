
clc
clear
close all

%% Array factor caculation ============================================================ 

M=20;%input('Enter the number of unitcells in one row: ');
N=M;
k=1;
j=1;  
c=3e8;
f=4.05e12;%input('Enter the operation frequency: ');
landa=c./f;
k=2*pi*f/c;
dx=landa/2.9608;
dy=dx;
A=M*dx*N*dy;
eta=120*pi;

r=200;
phi=linspace(0,2*pi,r);
theta=linspace(0,pi/2,r);
dtheta=theta(2)-theta(1);
dphi=phi(2)-phi(1);

%% phase calculation ================================================================
% desired main lobes directions in degree

tetades1=input('Enter desired first lobe direction in degree: ');
tetades2=input('Enter desired second lobe direction in degree: ');
tetades3=50;%input('Enter desired third lobe direction in degree: ');
tetades4=50;%input('Enter desired fourth lobe direction in degree: ');

%% supercell theta =====================================================================
% tetades11=asin(nnn*sind(tetades1));

kx1=k*sind(tetades1);
kx2=k*sind(tetades2);
kx3=k*sind(tetades3);
kx4=k*sind(tetades4);
x=dx*(1:N);
%% ======================================================================================
phase1=exp(1i*kx1*x);
% phase1=fliplr(phase1);
% cat phase11
phase11=phase1;
for h=1:N-1
    phase11=cat(1,phase11,phase1);
end
%% ======================================================================================
phase2=exp(1i*kx2*x);
% phase2=fliplr(phase2);
% cat phase22 =
phase22=phase2;
for h=1:N-1
    phase22=cat(1,phase22,phase2);
end
%% ======================================================================================
phase3=exp(1i*kx3*x);
% phase2=fliplr(phase2);
% cat phase33 
phase33=phase3;
for h=1:N-1
    phase33=cat(1,phase33,phase3);
end
%% ================================================================================
phase4=exp(1i*kx4*x);
% phase2=fliplr(phase2);
% cat phase44
phase44=phase4;
for h=1:N-1
    phase44=cat(1,phase44,phase4);
end

%% nimsaz phase1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=phase1;
NN=M;
for ke=1:NN
    B1(1,ke)=A(1,ke);
end
for i=1:NN-1
    for j=1:NN
        if j==1
            B1(i+1,j+NN-1)=B1(i,j);
        else
        B1(i+1,j-1)=B1(i,j);
        end
    end
end
 
%% nimsaz phase2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=phase2;
NN=M;
for ke=1:NN
    B2(1,ke)=A(1,ke);
end
for i=1:NN-1
    for j=1:NN
        if j==1
            B2(i+1,j+NN-1)=B2(i,j);
        else
        B2(i+1,j-1)=B2(i,j);
        end
    end
end

%% nimsaz phase3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=phase3;
NN=M;
for ke=1:NN
    B3(1,ke)=A(1,ke);
end
for i=1:NN-1
    for j=1:NN
        if j==1
            B3(i+1,j+NN-1)=B3(i,j);
        else
        B3(i+1,j-1)=B3(i,j);
        end
    end
end
C3=rot90(B3);
V3=rot90(C3);
F3=rot90(V3);

%% nimsaz phase4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=phase4;
NN=M;
for ke=1:NN
    B4(1,ke)=A(1,ke);
end
for i=1:NN-1
    for j=1:NN
        if j==1
            B4(i+1,j+NN-1)=B4(i,j);
        else
        B4(i+1,j-1)=B4(i,j);
        end
    end
end
C4=rot90(B4);
V4=rot90(C4);
F4=rot90(V4);


%% rot90 ==================================================================================

g1=rot90(phase11);
g11=rot90(rot90(phase11));
g2=rot90(phase22);
g22=rot90(rot90(phase22));
g3=rot90(phase33);
g33=rot90(rot90(phase33));
g4=rot90(phase44);
g44=rot90(rot90(phase44));
C1=rot90(B1);
V1=rot90(C1);
F1=rot90(V1);
C2=rot90(B2);
V2=rot90(C2);
F2=rot90(V2);

%% %%%%%%%%%%%%%%%%%%%%%supercell
% nnn=4;
% ppp=N;
% for ki=1:ppp
%     for kj=1:ppp
%         Tphase11(nnn*(ki-1)+1:nnn*(ki-1)+nnn ,nnn*(kj-1)+1:nnn*(kj-1)+nnn)=phase11(ki,kj);
%     end 
% end
% nnn=4;
% ppp=N;
% for ki=1:ppp
%     for kj=1:ppp
%         Tgg1(nnn*(ki-1)+1:nnn*(ki-1)+nnn ,nnn*(kj-1)+1:nnn*(kj-1)+nnn)=gg1(ki,kj);
%     end 
% end


%% Array factor ===========================================================================

a=input('Enter first coefficient between 0 and 1: ');
b=input('Enter second coefficient between 0 and 1: ');
c=0;%input('Enter third coefficient between 0 and 1: ');
d=0;%input('Enter forth coefficient between 0 and 1: ');
% totalph1=(a*phase11);
totalph1=(a*phase11+b*g2+c*g3+d*g4);
% totalph1=(a*exp(1i*c1))+(b*exp(1i*c2));
totalph1=totalph1/max(abs(totalph1(:)));
 AC=abs(totalph1);
% AC=ACq;
 phimn=angle(totalph1);
%  x=ones(M)*inv(AC1);
%  x1=x*a;
%  x2=x*b;
%  totalph=(x1*ones(M).*phase11)+(x2*ones(M).*B2);
%  totalph=totalph/max(abs(totalph(:)));
%   phimn=angle(totalph);
%   AC=abs(totalph);

%% =================2bit-discrete phase========================================================

 for i=1:N
     for j=1:N
         if phimn(i,j)>=0 && phimn(i,j)<pi/2
           phimn(i,j)=0;     
         elseif phimn(i,j)>=pi/2 && phimn(i,j)<pi
           phimn(i,j)=pi/2;
            
         elseif phimn(i,j)>=pi && phimn(i,j)<3*pi/2
           phimn(i,j)=pi;
         
         else
           phimn(i,j)=3*pi/2;
         end
 
     end
 end
 ddd=phimn;
nnn=5;
ppp=N;
for ki=1:ppp
    for kj=1:ppp
        Tphimn(nnn*(ki-1)+1:nnn*(ki-1)+nnn ,nnn*(kj-1)+1:nnn*(kj-1)+nnn)=ddd(ki,kj);
    end 
end

%% ========================= 2bit-discrete amplitute==============================

% 
%  for i=1:N
%   for j=1:N
%           if AC(i,j)>=0 && AC(i,j)<0.25
%            AC(i,j)=0.2;     
%          elseif AC(i,j)>=0.25 && AC(i,j)<0.5
%           AC(i,j)=0.4;
%             
%           elseif AC(i,j)>=0.5 && AC(i,j)<0.75
%            AC(i,j)=0.6;
%          
%            else
%            AC(i,j)=0.8;
%           end
%   end
%  end

%% ========================3bit-discrete phase=====================================

%  for i=1:N
%      for j=1:N
%           if phimn(i,j)>=-pi && phimn(i,j)<-3*pi/4
%            phimn(i,j)=-7*pi/8;     
%          elseif phimn(i,j)>=-3*pi/4 && phimn(i,j)<-pi/2
%          phimn(i,j)=-5*pi/8;
%             
%           elseif phimn(i,j)>=-pi/2 && phimn(i,j)<-pi/4
%            phimn(i,j)=-3*pi/8;
%             elseif phimn(i,j)>=-pi/4 && phimn(i,j)<0
%            phimn(i,j)=-pi/8;
%             elseif phimn(i,j)>=0 && phimn(i,j)<pi/4
%            phimn(i,j)=pi/8;
%             elseif phimn(i,j)>=pi/4 && phimn(i,j)<pi/2
%           phimn(i,j)=3*pi/8;
%             elseif phimn(i,j)>=pi/2 && phimn(i,j)<3*pi/4
%            phimn(i,j)=5*pi/8;
%           else 
%            phimn(i,j)=7*pi/8;
%          
%           end
%      end
%  end

%% ==========================3bit-discrete amplitude===========================================

  for i=1:N
      for j=1:N
  
          if AC(i,j)>=0 && AC(i,j)<0.09375
            AC(i,j)=0.089625;     
          elseif AC(i,j)>=0.09375 && AC(i,j)<0.1875
            AC(i,j)=0.17925;
          elseif AC(i,j)>=0.1875 && AC(i,j)<0.28125
            AC(i,j)=0.268875;
          elseif AC(i,j)>=0.28125 && AC(i,j)<0.375
            AC(i,j)=0.3585;
          elseif AC(i,j)>=0.375 && AC(i,j)<0.46875
            AC(i,j)=0.448125;
          elseif AC(i,j)>=0.46875 && AC(i,j)<0.5625
            AC(i,j)=0.53775;
            
          elseif AC(i,j)>=0.5625 && AC(i,j)<0.65625
            AC(i,j)=0.627375;
         
          else
            AC(i,j)=0.717;
          end
      end
  end

    for i=1:N
              plot(totalph1(i,:),'*')
              hold on
     end
          figure

%% ========================================================================================

for u=1:length(phi)
    for v=1:length(theta)
        E=0;
        for m=1:M
            for n=1:N      
                E=E+AC(m,n)*exp(+1i*k*dx*(n-1)*sin(theta(v))*cos(phi(u))+1i*k*dy*(m-1)*sin(theta(v))*sin(phi(u))+1i*phimn(m,n));
            end
        end
%         AF(u,v)=E*cos(theta(v));
        AF(u,v)=E;
    end
end
%e=FFT2([k*sin(theta).*cos(phi) k*sin(theta).*sin(phi)]);

uu=2*pi/landa*dx*sin(theta).*cos(phi);
vv=2*pi/landa*dx*sin(theta).*sin(phi);
[r1,t,y]=cart2pol(uu,vv,abs(AF)');
% figure(10)
% imagesc(r1,t,y)
%  figure
%%%%%%%%%%%%%%%%%%Directivity%%%%%%%%%%
Prad=0;
for u=1:length(phi)
    for v=1:length(theta)
        Prad=Prad+(abs(AF(u,v)))^2*sin(theta(v));
    end
end
Pradiation=dphi*dtheta*Prad;
for u=1:length(phi)
    for v=1:length(theta)
        dir(u,v)=4*pi*abs(AF(u,v))^2/Pradiation;
    end
end
wat=0;
for u=1:length(phi)
    for v=1:length(theta)
        wat=wat+(dir(u,v))*(sin(theta(v)));
    end
end
wat=wat*dtheta*dphi;
%RCSarray=4*pi*AF(1,1)^2
% RCSLPlate=4*pi*A^2/(landa^2);
% Dplate=4*pi*A/(landa^2);
% % RCSplate=10*log10(RCSLPlate)
% RCSR=10*log10(dir(1,1)/Dplate);
  mdir=max(max(dir));
  mE=max(max(abs(AF)));
% maxplatedir=4*pi*A/(landa^2)
% FOM=mdir/(maxplatedir)
figure(111)
imagesc(theta*180/pi,phi*180/pi,dir)
colormap(jet);

[phi2,theta2] = meshgrid(phi,theta);
[XX,YY,ZZ]=sph2cart(phi2,pi/2-theta2,dir');
% figure(2)
% imagesc(ZZ)
% figure
% subplot(2,2,1)
% surf(XX,YY,ZZ);
% colorbar
% xlabel('x')
% ylabel('y')
% zlabel('Directivity')
% title('Scattering Pattern')
% subplot(2,2,2)
% phiphi=linspace(0,360,length(phi));
% thetatheta=linspace(0,90,length(theta));
% imagesc(phiphi,thetatheta,dir')
% colorbar 
% xlabel('phi')
% ylabel('theta')
% title('Directivity')
% [phim thetam]=find(dir==mdir);
% phimm=phi(phim)*180/pi
% thetamm=theta(thetam)*180/pi
% dirm=dir(phim,thetam)
% phimn=(phimn*180/pi);
figure 
% phimn=flipud(phimn)+157.5; 
imagesc(phimn)
colorbar
title('Surface Phase Distribution')
figure 
AC=flipud(AC); 
imagesc(AC)
colorbar
title('Surface amplitute Distribution')
dirv=reshape(dir,[1 r*r]);

ArrayFactormax=max(max(abs(AF)))
figure
surf(XX,YY,ZZ);
figure
% plot(theta*180/pi,10*log10(dir(180*r/360,:)),'linewidth',2)
% hold on
% plot(theta*180/pi,10*log10(dir(270*r/360,:)),'linewidth',2)
% set(gca,'YLim',[0 40])
% hold on
% plot(theta*180/pi,dir(135*r/360,:),'linewidth',2)
% hold on
plot(theta*180/pi,dir(180*r/360,:),'linewidth',2)
hold on
plot(theta*180/pi,dir(90*r/360,:),'linewidth',2)
% hold on
% plot(theta*180/pi,dir(180*r/360,:),'linewidth',2)

% W11=W1*180/pi;
% UU=180/pi.*x;

save('AC.mat','AC')
save('phimn.mat','phimn')


