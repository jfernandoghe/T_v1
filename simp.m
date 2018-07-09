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
qgamas=3;
pung=3;
ne=3;
npe=pung-2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Random Walk Parameters
% ex, ey = Boundaries
% n = number of points to generate
% npoints = number of scanpath points to evaluate as pattern
    [Sal1, Sal2]=RandomWalk(50, 50, pung, maxFrechet);
%     [Sal1, Sal2]=RandomWalk2(50, 50, pung, ne, npe, maxFrechet);
    qx =Sal1(1,:)
    qy =Sal1(2,:)
    qpx=Sal2(1,:)
    qpy=Sal2(2,:)
    hold on;
    grid on
    grid minor
    plot(qx,qy,'.-.r')
    plot(qx,qy,'*r')
    for ii = 1:length(qx)
        text(qx(ii),qy(ii)+5,num2str(ii),'Color','r')
    end
    hold on;
    grid on
    grid minor
    plot(qpx,qpy,'.-.b')
    plot(qpx,qpy,'*b')
    for ii = 1:length(qpx)
        text(qpx(ii),qpy(ii)+5,num2str(ii),'Color','b')
    end
    
    
%     TWO_Esimp(qx, qy, maxFrechet)
end

function [Sal1, Sal2] = RandomWalk (ex, ey, n, maxFrechet)
% % % % % 1920x1080 Screen Resolution
% % % % % First random point within a square on the middle on the screen
% % % % % emulation atention to test
    xmin=0+ex;      xmax=1920-ex;       %ex px offset from top to bottom
    ymin=0+ey;      ymax=1080-ey;       %ey px offset from left to right
    m=n;                                %Number of points to generate
    Sal1=[]; Sal2=[]; s1=[]; s2=[]; s3=[]; s4=[];
    i=1; k=1;
% % Random Coordinates n=m
%     x1=round(xmin+rand(1,n)*(xmax-xmin))
%     y1=round(ymin+rand(1,m)*(ymax-ymin))
    
    s1=round(xmin+rand(1,1)*(xmax-xmin))
    s2=round(ymin+rand(1,1)*(ymax-ymin))
%     x2=round((randn(1,n)*30)+x1);
%     y2=round((randn(1,m)*30)+y1); 
    while i~=n+1
        pro=abs(rand(1));
        if pro>0.4
            s1=horzcat(s1, (s1(k)+(randn(1,1)*30)))
            s2=horzcat(s2, (s1(k)+(randn(1,1)*30)))
            k=k+1;
        else
            s1=horzcat(s1,round(xmin+rand(1,1)*(xmax-xmin)))
            s2=horzcat(s2,round(xmin+rand(1,1)*(xmax-xmin)))
            s3=horzcat(s3, (s1(k)+(randn(1,1)*30)));
            s4=horzcat(s4, (s2(k)+(randn(1,1)*30)));
            i=i+1;
        end
    end
    s1
    s2
%     Sal1(1,:)=x1;
%     Sal1(2,:)=y1;
    Sal1(1,:)=s1;
    Sal1(2,:)=s2;
    Sal2(1,:)=s3;
    Sal2(2,:)=s4;
%     Sal2(1,:)=real(x2);
%     Sal2(2,:)=real(y2);
end



function [Sal1, Sal2, outl] = RandomWalk2 (ex, ey, n, ne, npe, maxFrechet)
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
    x2=round((randn(1,n)*30)+x1)
    y2=round((randn(1,m)*30)+y1)
    e = randi([1 n-npe],1,1);    
    sx=x2(1:e);
    sy=y2(1:e);
    ex=x2(e:end);
    ey=y2(e:end);
    m=[];
    for i=0:(npe-1)
        er=e+i;
        mx=round(x1(er)+(rand(1,ne)*maxFrechet));
        my=round(y1(er)+(rand(1,ne)*maxFrechet));
%         tx=horzcat(sx, ex);
%         ty=horzcat(sy, ey);
%         m=horzcat(m,mx);
        Sal1(1,:)=horzcat(sx,mx,ex);
        Sal1(2,:)=horzcat(sy,my,ey);
%         Sal1(1,:)=horzcat(x1(1:er), round((x1(er+i)-maxFrechet)+rand(1,ne)*maxFrechet), x1(er:end))
%         Sal1(2,:)=horzcat(y1(1:er), round((y1(er+i)-maxFrechet)+rand(1,ne)*maxFrechet), y1(er:end))
    end
% Add random outlier
    top=fix(n/6)
%     for i=1:top
%        pla=randi(n,1,1);
%        x2(1,pla)=round(xmin+rand(1,1)*(xmax-xmin))       
%        y2(1,pla)=round(ymin+rand(1,1)*(ymax-ymin))
%     end
% To a Sal vector
%     Sal1(1,:)=real(x1);
%     Sal1(2,:)=real(y1);
    Sal2(1,:)=real(x2);
    Sal2(2,:)=real(y2);
    outl=top;
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