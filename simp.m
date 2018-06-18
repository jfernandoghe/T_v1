function ssimp ()
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Coded by: Fernando Gonzalez-Herrera   - CIMAT Zacatecas
%           Carlos Lara-Alvarez         - CIMAT Zacatecas
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clc;
close all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Parameters for Frechet and scanpath comparison algorith
maxFrechet=(3*30)*sqrt(2);
maxA=maxFrechet;
dista=maxFrechet;
dista2=maxFrechet;
qgamas=2;
pung=4;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Random Walk Parameters
% ex, ey = Boundaries
% n = number of points to generate
% npoints = number of scanpath points to evaluate as pattern
    Sal=RandomWalk(50,50,pung);
    qx =Sal(1,:)
    qy =Sal(2,:)
    qpx=Sal(3,:)
    qpy=Sal(4,:)
    hold on;
    grid on
    grid minor
    plot(qx,qy,'.-.r')
    plot(qx,qy,'*r')
    for ii = 1:length(qx)
        text(qx(ii),qy(ii),num2str(ii),'Color','k')
    end
    TWO_Esimp(qx, qy, maxFrechet)
end

function [Sal] = RandomWalk (ex, ey, n)
% % % % % 1920x1080 Screen Resolution
% % % % % First random point within a square on the middle on the screen
% % % % % emulation atention to test
    xmin=0+ex;      xmax=1920-ex;       %ex px offset from top to bottom
    ymin=0+ey;      ymax=1080-ey;       %ey px offset from left to right
    m=n;                                %Number of points to generate
% % Random Coordinates n=m
    x1=round(xmin+rand(1,n)*(xmax-xmin));
    y1=round(ymin+rand(1,m)*(ymax-ymin));
% % Add gaussian noise to each coordenate
    x2=round((randn(1,n)*30)+x1);
    y2=round((randn(1,m)*30)+y1);
% % Add random outlier
    err=abs(rand(1));
    if (err>0.999)
       pla=randi(n,1,1);
       x2(1,pla)=round(xmin+rand(1,1)*(xmax-xmin));       
       y2(1,pla)=round(ymin+rand(1,1)*(ymax-ymin));
    end
% To a Sal vector
    Sal(1,:)=real(x1);
    Sal(2,:)=real(y1);
    Sal(3,:)=real(x2);
    Sal(4,:)=real(y2);
end

function [simp] = TWO_Esimp (scx, scy, maxFrechet)
    if size(scx,2)==size(scy,2)
        i0=2; i1=size(scx,2);
        for i=i0:i1
            fprintf('###############   Run %i      ############### \n',i);
            pi1=[scx(1);scy(1);1];
            pf1=[scx(i);scy(i);1];
            fprintf('###############   Recta entre (%i, %i)  a (%i, %i)\n',scx(1),scy(1),scx(i),scy(i));
            l1=cross(pi1,pf1);
            ll1=l1(1:2,1);
            k1=1/sqrt(dot(ll1,ll1));
            eq1=k1*l1;
            for i=1:i1-i0
                fprintf('###############   Puntos %i',1);
                l2dots(i)=dot(eq1,[scx(i);scy(i);1]);
            end
            l2dots=abs(l2dots)
            l2dot=l2dots;
            l2dot=l2dot(l2dot>0);
            l2dot=l2dot(l2dot<=maxFrechet)
            for i=1:size(l2dots,2)
                for j=1:size(l2dot,2)
                    find(l2dots(i)==l2dot(j))
                end
            end
        end
    end
end

% function [PoI] = poiR2(vx1,vy1,i0,i1,maxleng)
% % vx1       = X coordenate Vector
% % vy1       = Y coordenate Vector
% % i0        = Starting point
% % i1        = Ending point
% % maxleng   = Max distance
% 
% if i1-i0<=1
%     PoI=0;
%     return;
% end
% i=1;j=1; u=1;
% pi1=[vx1(i0);vy1(i0);1];
% pf1=[vx1(i1);vy1(i1);1];
% l1=cross(pi1,pf1);
% ll1=l1(1:2,1);
% k1=1/sqrt(dot(ll1,ll1));
% eq1=k1*l1;                      %Line A - Thru first and last point
% l2dots=zeros(i1-i0+1,1);
% for i=i0:1:i1
%     l2dots(i-i0+1)=dot(eq1,[vx1(i);vy1(i);1]);
% end
% l2dots=abs(l2dots);
% [mL,iL]=max(l2dots);            %Furthest point from line A
% % a=[vx1(iL)-vx1(i0);vy1(iL)-vy1(i0);0];
% % b=[vx1(i1)-vx1(i0);vy1(i1)-vy1(i0);0];
% iLL=iL+i0-1;
% A=[vx1(iLL)-vx1(i0),vy1(iLL)-vy1(i0);
%     vx1(i1)-vx1(i0),vy1(i1)-vy1(i0)];
% areap=abs(0.5*det(A));
% if (areap>maxleng)    
%     PoI=iL+i0-1;
% else
%     PoI=0;
% end
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     for i=1:size(v0,2)
%             if i==size(v0,2)
%                 f = find(w>v0(i),1);
%             else
%                 f = find(w>v0(i)&w<=v1(i),1);
%             end
%     end