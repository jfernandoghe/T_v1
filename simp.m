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
qgamas=2;               %n-Gramas
pung=4;                 %Points to generate
ne=2;                   %Number of errors per point
npe=3;                  %Number of points with error
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Random Walk Parameters
% ex, ey = Boundaries
% n = number of points to generate
% npoints = number of scanpath points to evaluate as pattern
%     Sal=RandomWalk(50,50,pung);
    [Sal1,Sal2]=RandomWalk2(50, 50, pung, npe, ne);
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
        text(qx(ii),qy(ii),num2str(ii),'Color','r')
    end
    hold on;
    grid on
    grid minor
    plot(qpx,qpy,'.-.b')
    plot(qpx,qpy,'*b')
    for ii = 1:length(qpx)
        text(qpx(ii),qpy(ii),num2str(ii),'Color','b')
    end
    
    
%     TWO_Esimp(qx, qy, maxFrechet)
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

function [Sal1, Sal2] = RandomWalk2 (ex, ey, n, npe, ne)
% % % % % 1920x1080 Screen Resolution
% % % % % First random point within a square on the middle on the screen
% % % % % emulation atention to test
    xmin=0+ex;      xmax=1920-ex;       %ex px offset from top to bottom
    ymin=0+ey;      ymax=1080-ey;       %ey px offset from left to right
    m=n;                                %Number of points to generate
    k=1;
% % Random Coordinates n=m
    x1=round(xmin+rand(1,n)*(xmax-xmin))
    y1=round(ymin+rand(1,m)*(ymax-ymin))
% % Add gaussian noise to each coordenate
    x2=round((randn(1,n)*30)+x1);
    y2=round((randn(1,m)*30)+y1); 
    er = randi([1 n-npe],1,1);
    x1p=ones(1,size(x1,2)); y1p=ones(1,size(x1,2));
    for k=1:npe
        if (k==1)
            inix=horzcat(x1(1:er));
            iniy=horzcat(y1(1:er));
        end
        for k=2:npe
            medx=round((x1(er)-10)+rand(1,ne)*(1920/10));
            medy=round((y1(er)-10)+rand(1,ne)*(1920/10));
        end
        if (k==npe)
            finx=horzcat(x1(er+1:end));
            finy=horzcat(y1(er+1:end));
        end 
        a=horzcat(inix,medx,finx);
        b=horzcat(iniy,medy,finy);
        a=horzcat(a);
        b=horzcat(b);
    Sal1(1,:)=a
    Sal1(2,:)=b
%             Sal1(1,:)=horzcat(x1(1:er), round((x1(er)-10)+rand(1,ne)*(1920/10)),    round((x1(er+1)-10)+rand(1,ne+1)*(1920/10)),   round((x1(er+2)-10)+rand(1,ne+2)*(1920/10)),     x1(er+3:end))
%             Sal1(2,:)=horzcat(y1(1:er), round((y1(er)-10)+rand(1,ne)*(1920/10)),    round((y1(er+1)-10)+rand(1,ne+1)*(1920/10)),   round((y1(er+2)-10)+rand(1,ne+2)*(1920/10)),     y1(er+3:end))
     
    end
% Add random outlier
    err=randi([1 n],1,1);
    top=fix(n/5)
    for i=1:top
       pla=randi(n,1,1);
       x1(1,pla)=round(xmin+rand(1,1)*(xmax-xmin))       
       y1(1,pla)=round(ymin+rand(1,1)*(ymax-ymin))
    end
% To a Sal vector
%     Sal1(1,:)=real(x1);
%     Sal1(2,:)=real(y1);
    Sal2(1,:)=real(x2);
    Sal2(2,:)=real(y2);
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