function alg ()
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Coded by: Fernando Gonzalez-Herrera   - CIMAT Zacatecas
%           Carlos Lara-Alvarez         - CIMAT Zacatecas
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clc;
global perror;
perror=0;
close all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Parameters for Frechet and scanpath comparison algorith
maxFrechet=(3*30)*sqrt(2);
maxA=maxFrechet;
dista=maxFrechet;
dista2=maxFrechet;
qgamas=3;
pung=5;
ne=2;
% frechetsal=zeros(1,pung-qgamas+1);
tim=1;
% % Make all the combinations for n length scanpaths in a k-waise maner
% % Pairwaise maner
    allc=nchoosek(1:pung,qgamas);
    frechetsal=zeros(1,length(allc));
    vfrechet=zeros(tim, length(allc));    
    vfrechetV=0;
    hc=zeros(pung-1,qgamas);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for q=1:(pung-1)
    hc(q,1)=q;
    hc(q,2)=q+1;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Number of iterations to evaluate
for z=1:tim
    clc;
    close all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Random Walk Parameters
% ex, ey = Boundaries
% n = number of points to generate
% npoints = number of scanpath points to evaluate as pattern
    [Sal1, Sal2]=RandomWalk2(50,50,pung, ne);
    qx =Sal1(1,:)
    qy =Sal1(2,:)
    qpx=Sal2(1,:)
    qpy=Sal2(2,:)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for qgam=1:size(allc,1)
%     vcom(qgam,1)=qgam;
%     vcom(qgam,2)=qgam+qgamas-1;
% % Selecting the qgamas combinatios till length(qx)
        x=qx(allc(qgam,1):allc(qgam,2));
        y=qy(allc(qgam,1):allc(qgam,2));
        px=qpx(allc(qgam,1):allc(qgam,2));
        py=qpy(allc(qgam,1):allc(qgam,2));
        fprintf('Q GAMAS : \t\t%f \t%f',allc(qgam,1),allc(qgam,2))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Main algoritm to compare Scanpaths
        inits=zeros(1,size(y,2));
        ends=zeros(1,size(y,2));
        inits(1)=1;
        ends(1)=size(y,2);
        pos=1;
        iMiss=1;
        while 1
            ii=poiR2(x,y,inits(pos),ends(pos),maxA);
            if ii>0
                inits(iMiss+1)=inits(pos);
                ends(iMiss+1)=ii;
                inits(iMiss+2)=ii;
                ends(iMiss+2)=ends(pos);
                iMiss=iMiss+2;
            end
            if(pos==iMiss)
                break;
            end
            pos=pos+1;
        end
        POIval=unique(horzcat(inits(1:iMiss),ends(1:iMiss)))
        xx=x(POIval)
        yy=y(POIval)
        distsini=arrayfun(@(x,y) wdistance(x,y,px(1),py(1),1,1), xx,yy)
        v=find(distsini<dista)
        v0=v(1:end)
        v1=v(2:end)
        if v==1
            v0=1;v1=1;
        end
        distsfin=arrayfun(@(x,y) wdistance(x,y,px(end),py(end),1,1), xx,yy)
        w=find(distsfin<dista2)
        figure(1)
        grid on;
        hold on;
        grid minor;
        plot(x,y,'.b');
        plot(x,y,'--b');
        plot(px,py,'m');
        plot(x(POIval),y(POIval),'or');
        lbl_dwn = .1*max(yy);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Add plot numbers to each iteration of q-gamas in a q wize manner
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%     for i = 1:length(POIval)
%         pp=POIval(i);
%         % % Condition for index exceed inner matrix dimmensions
%         if (pp>length(xx))
%             %             fprintf('TEST01');
%             pp=length(xx);
%         end
%         plot(xx(pp),yy(pp),'r+');
%         % % Label the points with the corresponding 'x' value
% %         tex=strcat(num2str(i),'-',num2str(qgam));
%         tex=strcat(num2str(i));
%         text(xx(pp),yy(pp)+0.05,tex);
%     end
%     for i = 1:length(px)
%         plot(px(i),py(i),'b+');
%         % % Label the points with the corresponding 'y' value
% %         tex=strcat(num2str(i),'-',num2str(qgam));
%         tex=strcat(num2str(i));
%         text(px(i),py(i)+0.05,tex);
%     end
        hold on;
        for i=1:size(v0,2)
            if i==size(v0,2)
                f = find(w>v0(i),1);
            else
                f = find(w>v0(i)&w<=v1(i),1);
            end
            if length(f)==1
            start=POIval(v0(i))
            endd=POIval(w(f))
            sx=x(start:endd);
            sy=y(start:endd);
            Pe=(vertcat(px,py));
            Qu=(vertcat(sx,sy));
            figure(1)
            plot(sx,sy,'-hg');
            axis([0 2000 0 1300]);
            fdic=frechet_decide(Pe,Qu,maxFrechet,0,0);
            if fdic==1
                figure(1);
                grid minor;
                    plot(sx,sy,'-k');
                else
                    frechetsal(1,qgam)=frechetsal(1,qgam)+1;
                    fprintf('EEEEEEEEEERRRROOOR\n');
                end
            else
                frechetsal(1,qgam)=frechetsal(1,qgam)+1;
            end
        end
end
vfrechet(z,:)=frechetsal
vfrechet;
%     jx=1-pdist(vertcat(x,px),'jaccard')
%     jy=1-pdist(vertcat(y,py),'jaccard')
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % End of iterative statements for evaluation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for ii = 1:length(qx)
    text(qx(ii),qy(ii),num2str(ii),'color','k')
end
for ii = 1:length(qpx)
    text(qpx(ii),qpy(ii),num2str(ii),'color','k')
end
allc'
vfrechet
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for band=1:size(vfrechet,1)
    if (length(find(vfrechet==1))>ceil(0.2*length(allc)))
        if (vfrechet(1,1)+vfrechet(1,5)+vfrechet(1,8)+vfrechet(1,10)) >2
            vfrechetV=vfrechetV+1;
        end
%         vfrechet(band,size(vfrechet,2)+1)=vfrechet(band,size(vfrechet,2)+1)+1;
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if vfrechetV==0
    fprintf('\nEl scanpath patron cumple con el scanpath generado\t\n\n');
    fprintf('%f errores de %f comparaciones, %i \t\n',size(vfrechet,2)*size(vfrechet,1)-sum(vfrechet(:)==0),size(vfrechet,2)*size(vfrechet,1),(((sum(vfrechet(:)==0))/(size(vfrechet,2)*size(vfrechet,1))))*100);
else
    fprintf('\nEl scanpath patron NO cumple con el scanpath generado\t\n\n');
    fprintf('%f errores de %f comparaciones, %i \t\n',size(vfrechet,2)*size(vfrechet,1)-sum(vfrechet(:)==0),size(vfrechet,2)*size(vfrechet,1),(((sum(vfrechet(:)==0))/(size(vfrechet,2)*size(vfrechet,1))))*100);
end
figure (2)
plot(qx,qy,'r');
hold on;
plot(qpx,qpy,'b');
grid on;
grid minor;
axis([0 2000 0 1300]);
for ii = 1:length(qx)
    text(qx(ii),qy(ii),num2str(ii),'color','r')
end
for ii = 1:length(qpx)
    text(qpx(ii),qpy(ii),num2str(ii),'color','b')
end
hold on;
end

function dis = wdistance(x0,y0,x1,y1,alx,aly)
% % Get the distance between two points with compensation alx,aly
% % 1 & COMPY, only compensation in the Y axis
    dx=(x1-x0)*alx;
    dy=(y1-y0)*aly;
    dis=sqrt(dx*dx+dy*dy);
end

function [PoI] = poiR2(vx1,vy1,i0,i1,maxleng)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Coded by: Jose Fernando Gonzalez Herrera - CIMAT Zacatecas
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% vx1       = X coordenate Vector
% vy1       = Y coordenate Vector
% i0        = Starting point
% i1        = Ending point
% maxleng   = Max distance
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% NOTE:     vx1 and vy1 must be equal size
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if i1-i0<=1
    PoI=0;
    return;
end
i=1;j=1; u=1;
pi1=[vx1(i0);vy1(i0);1];
pf1=[vx1(i1);vy1(i1);1];
l1=cross(pi1,pf1);
ll1=l1(1:2,1);
k1=1/sqrt(dot(ll1,ll1));
eq1=k1*l1;                      %Line A - Thru first and last point
l2dots=zeros(i1-i0+1,1);
for i=i0:1:i1
    l2dots(i-i0+1)=dot(eq1,[vx1(i);vy1(i);1]);
end
l2dots=abs(l2dots);
[mL,iL]=max(l2dots);            %Furthest point from line A
% a=[vx1(iL)-vx1(i0);vy1(iL)-vy1(i0);0];
% b=[vx1(i1)-vx1(i0);vy1(i1)-vy1(i0);0];
iLL=iL+i0-1;
A=[vx1(iLL)-vx1(i0),vy1(iLL)-vy1(i0);
    vx1(i1)-vx1(i0),vy1(i1)-vy1(i0)];
areap=abs(0.5*det(A));
if (areap>maxleng)    
    PoI=iL+i0-1;
else
    PoI=0;
end
end

function [Sal1, Sal2] = RandomWalk (ex, ey, n, ne)
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
%     er = randi([2 n-1],1,1);
%     x2(1:er);
%     y2(1:er);
%     x2(er+1:end);
%     y2(er+1:end);
%     Sal1(1,:)=horzcat(x1(1:er), round((x1(er)-5)+rand(1,ne)*(x1(er)+5)), x1(er+1:end));
%     Sal1(2,:)=horzcat(y1(1:er), round((y1(er)-5)+rand(1,ne)*(y1(er)+5)), y1(er+1:end));
% Add random outlier
    err=randi([1 n],1,1);
    top=fix(n/5)
    for i=1:top
       pla=randi(n,1,1);
       x2(1,pla)=round(xmin+rand(1,1)*(xmax-xmin))       
       y2(1,pla)=round(ymin+rand(1,1)*(ymax-ymin))
    end
% To a Sal vector
    Sal1(1,:)=real(x1);
    Sal1(2,:)=real(y1);
    Sal2(1,:)=real(x2);
    Sal2(2,:)=real(y2);
end

function [Sal1, Sal2] = RandomWalk2 (ex, ey, n, ne)
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
    er = randi([2 n-1],1,1);
    x2(1:er);
    y2(1:er);
    x2(er+1:end);
    y2(er+1:end);
    Sal1(1,:)=horzcat(x1(1:er), round((x1(er)-5)+rand(1,ne)*(1920/10)), x1(er+1:end));
    Sal1(2,:)=horzcat(y1(1:er), round((y1(er)-5)+rand(1,ne)*(1080/10)), y1(er+1:end));
% Add random outlier
    err=randi([1 n],1,1);
    top=fix(n/5)
    for i=1:top
       pla=randi(n,1,1);
       x2(1,pla)=round(xmin+rand(1,1)*(xmax-xmin))       
       y2(1,pla)=round(ymin+rand(1,1)*(ymax-ymin))
    end
% To a Sal vector
%     Sal1(1,:)=real(x1);
%     Sal1(2,:)=real(y1);
    Sal2(1,:)=real(x2);
    Sal2(2,:)=real(y2);
end

function [decide] = frechet_decide(P,Q,len,plotFSD,printFSD)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Frechet Distance Demo
% Script to run the Decision Problem (DP) for Frechet distance
% The DP is based on Algorithm 1 in "Computing the Frechet Distance 
% Between Two Polygonal Curves", Alt and Godau, International Journal of
% Computational Geometry and Applications, 1995.
% Author:     Rich Kenefic
% Created:    14-Dec-2011
% Modified:   18-Oct-2012 To include the Frechet computation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
global perror;
%--solves the decision problem for Frechet distance
[M,N]=size(P);
if M~=2, error('P must be a 2 by I array'); end
I = N;
[M,N]=size(Q);
if M~=2, error('Q must be a 2 by J array'); end
J = N;
%--compute the free space in each cell
for i=1:I-1
    xp = P(1,i); yp = P(2,i);
    xp1 = P(1,i+1); yp1 = P(2,i+1);
    lP(i) = (xp1-xp)^2 + (yp1-yp)^2;  %--length^2 of the ith segment of P
end
for j=1:J-1
    xq = Q(1,j); yq = Q(2,j);
    xq1 = Q(1,j+1); yq1 = Q(2,j+1);
    lQ(j) = (xq1-xq)^2 + (yq1-yq)^2;  %--length^2 of the jth segment of Q
end
for i=1:I
    xp = P(1,i); yp = P(2,i);
    for j=1:J
        xq = Q(1,j); yq = Q(2,j);
        lPQ(i,j) = (xq-xp)^2 + (yq-yp)^2;  %--length^2 of P->Q link
    end
end
for i=1:I-1
    xp = P(1,i); yp = P(2,i);
    xp1 = P(1,i+1); yp1 = P(2,i+1);    
    for j=1:J-1
        xq = Q(1,j); yq = Q(2,j);
        xq1 = Q(1,j+1); yq1 = Q(2,j+1);
        %--solve for the line segment circle intersections
        ap = lP(i);  aq = lQ(j);
        bp = 2*((xp-xq)*(xp1-xp)+(yp-yq)*(yp1-yp));
        bq = 2*((xq-xp)*(xq1-xq)+(yq-yp)*(yq1-yq));
        c = lPQ(i,j) - len^2;
        dp = bp*bp - 4*ap*c;
        dq = bq*bq - 4*aq*c;
        if dp<0.0,  %--  a_ij, b_ij, LF_ij
            %--line and circle do not intersect
            A(i,j) = NaN; B(i,j) = NaN; 
            LF(i,j,1) = NaN; LF(i,j,2) = NaN;
        else
            up = ((-bp)+sqrt(dp))/(2*ap); um = ((-bp)-sqrt(dp))/(2*ap);
            if (((up<0)&&(um<0))||((up>1)&&(um>1)))
                %--line segment outside circle
                A(i,j) = NaN; B(i,j) = NaN;
                LF(i,j,1) = NaN; LF(i,j,2) = NaN;
            elseif ((min([um up])<0)&&(max([um up])>1))
                %--line segment is interior to circle
                A(i,j) = 0; B(i,j) = 1;
                LF(i,j,1) = 0; LF(i,j,2) = 1;
            elseif ((min([um up])<=0)&&(max([um up])<=1))
                %--one intersection (b_i,j)
                A(i,j) = 0; B(i,j) = max([um up]);
                LF(i,j,1) = 0; LF(i,j,2) = B(i,j);
            elseif ((min([um up])>=0)&&(max([um up])>1))
                %--one intersection (a_i,j)
                A(i,j) = min([um up]); B(i,j) = 1;
                LF(i,j,1) = A(i,j); LF(i,j,2) = 1;
            elseif ((min([um up])>=0)&&(max([um up])<=1))
                %--two intersections
                A(i,j) = min([um up]); B(i,j) = max([um up]);
                LF(i,j,1) = A(i,j); LF(i,j,2) = B(i,j);
            else
                error('Unexpected case in frechet_decide at LF');
            end
        end
        if dq<0.0,  %--  c_ij, d_ij, BF_ij
            %--line and circle do not intersect
            C(i,j) = NaN; D(i,j) = NaN; 
            BF(i,j,1) = NaN; BF(i,j,2) = NaN;
        else
            up = (-bq+sqrt(dq))/(2*aq); um = (-bq-sqrt(dq))/(2*aq);
            if (((up<0)&&(um<0))||((up>1)&&(um>1)))
                %--line segment outside circle
                C(i,j) = NaN; D(i,j) = NaN;
                BF(i,j,1) = NaN; BF(i,j,2) = NaN;
            elseif ((min([um up])<0)&&(max([um up])>1))
                %--line segment is interior to circle
                C(i,j) = 0; D(i,j) = 1;
                BF(i,j,1) = 0; BF(i,j,2) = 1;
            elseif ((min([um up])<=0)&&(max([um up])<=1))
                %--one intersection (d_i,j)
                C(i,j) = 0; D(i,j) = max([um up]);
                BF(i,j,1) = 0; BF(i,j,2) = D(i,j);
            elseif ((min([um up])>=0)&&(max([um up])>1))
                %--one intersection (a_i,j)
                C(i,j) = min([um up]); D(i,j) = 1;
                BF(i,j,1) = C(i,j); BF(i,j,2) = 1;
            elseif ((min([um up])>=0)&&(max([um up])<=1))
                %--two intersections
                C(i,j) = min([um up]); D(i,j) = max([um up]);
                BF(i,j,1) = C(i,j); BF(i,j,2) = D(i,j);
            else
                error('Unexpected case in frechet_decide at BF');
            end
        end
    end
end
%--Top row
xp = P(1,I); yp = P(2,I);
for j=1:J-1
    xq = Q(1,j); yq = Q(2,j);
    xq1 = Q(1,j+1); yq1 = Q(2,j+1);
    %--solve for the line segment circle intersections
    aq = lQ(j);
    bq = 2*((xq-xp)*(xq1-xq)+(yq-yp)*(yq1-yq));
    c = lPQ(I,j) - len^2;
    dq = bq*bq - 4*aq*c;
    if dq<0.0,  %--  c_ij, d_ij, BF_ij
        %--line and circle do not intersect
        C(I,j) = NaN; D(I,j) = NaN;
        BF(I,j,1) = NaN; BF(I,j,2) = NaN;
    else
        up = (-bq+sqrt(dq))/(2*aq); um = (-bq-sqrt(dq))/(2*aq);
        if (((up<0)&&(um<0))||((up>1)&&(um>1)))
            %--line segment outside circle
            C(I,j) = NaN; D(I,j) = NaN;
            BF(I,j,1) = NaN; BF(I,j,2) = NaN;
        elseif ((min([um up])<0)&&(max([um up])>1))
            %--line segment is interior to circle
            C(I,j) = 0; D(I,j) = 1;
            BF(I,j,1) = 0; BF(I,j,2) = 1;
        elseif ((min([um up])<=0)&&(max([um up])<=1))
            %--one intersection (d_i,j)
            C(I,j) = 0; D(I,j) = max([um up]);
            BF(I,j,1) = 0; BF(I,j,2) = D(I,j);
        elseif ((min([um up])>=0)&&(max([um up])>1))
            %--one intersection (a_i,j)
            C(I,j) = min([um up]); D(I,j) = 1;
            BF(I,j,1) = C(I,j); BF(I,j,2) = 1;
        elseif ((min([um up])>=0)&&(max([um up])<=1))
            %--two intersections
            C(I,j) = min([um up]); D(I,j) = max([um up]);
            BF(I,j,1) = C(I,j); BF(I,j,2) = D(I,j);
        else
            error('Unexpected case in frechet_decide at BF (top row)');
        end
    end
end
%--Right column
xq = Q(1,J); yq = Q(2,J);
for i=1:I-1
    xp = P(1,i); yp = P(2,i);
    xp1 = P(1,i+1); yp1 = P(2,i+1);
    %--solve for the line segment circle intersections
    ap = lP(i); 
    bp = 2*((xp-xq)*(xp1-xp)+(yp-yq)*(yp1-yp));
    c = lPQ(i,J) - len^2;
    dp = bp*bp - 4*ap*c;
    if dp<0.0,  %--  a_ij, b_ij, LF_ij
        %--line and circle do not intersect
        A(i,J) = NaN; B(i,J) = NaN;
        LF(i,J,1) = NaN; LF(i,J,2) = NaN;
    else
        up = (-bp+sqrt(dp))/(2*ap); um = (-bp-sqrt(dp))/(2*ap);
        if (((up<0)&&(um<0))||((up>1)&&(um>1)))
            %--line segment outside circle
            A(i,J) = NaN; B(i,J) = NaN;
            LF(i,J,1) = NaN; LF(i,J,2) = NaN;
        elseif ((min([um up])<0)&&(max([um up])>1))
            %--line segment is interior to circle
            A(i,J) = 0; B(i,J) = 1;
            LF(i,J,1) = 0; LF(i,J,2) = 1;
        elseif ((min([um up])<=0)&&(max([um up])<=1))
            %--one intersection (b_i,j)
            A(i,J) = 0; B(i,J) = max([um up]);
            LF(i,J,1) = 0; LF(i,J,2) = B(i,J);
        elseif ((min([um up])>=0)&&(max([um up])>1))
            %--one intersection (a_i,j)
            A(i,J) = min([um up]); B(i,J) = 1;
            LF(i,J,1) = A(i,J); LF(i,J,2) = 1;
        elseif ((min([um up])>=0)&&(max([um up])<=1))
            %--two intersections
            A(i,J) = min([um up]); B(i,J) = max([um up]);
            LF(i,J,1) = A(i,J); LF(i,J,2) = B(i,J);
        else
            error('Unexpected case in frechet_decide at LF (right col)');
        end
    end
end
%--Top Right cell
xp = P(1,I-1); yp = P(2,I-1);
xp1 = P(1,I); yp1 = P(2,I);
xq = Q(1,J-1); yq = Q(2,J-1);
xq1 = Q(1,J); yq1 = Q(2,J);
ap = lP(I-1);  aq = lQ(J-1);
bp = 2*((xp-xq1)*(xp1-xp)+(yp-yq1)*(yp1-yp));
bq = 2*((xq-xp1)*(xq1-xq)+(yq-yp1)*(yq1-yq));
cp = lPQ(I-1,J) - len^2;
cq = lPQ(I,J-1) - len^2;
dp = bp*bp - 4*ap*cp;
dq = bq*bq - 4*aq*cq;
if dp<0.0,  %--  a_ij, b_ij, LF_ij
    %--line and circle do not intersect
    A(I,J) = NaN; B(I,J) = NaN;
    LF(I,J,1) = NaN; LF(I,J,2) = NaN;
else
    up = (-bp+sqrt(dp))/(2*ap); um = (-bp-sqrt(dp))/(2*ap);
    if (((up<0)&&(um<0))||((up>1)&&(um>1)))
        %--line segment outside circle
        A(I,J) = NaN; B(I,J) = NaN;
        LF(I,J,1) = NaN; LF(I,J,2) = NaN;
    elseif ((min([um up])<0)&&(max([um up])>1))
        %--line segment is interior to circle
        A(I,J) = 0; B(I,J) = 1;
        LF(I,J,1) = 0; LF(I,J,2) = 1;
    elseif ((min([um up])<=0)&&(max([um up])<=1))
        %--one intersection (b_i,j)
        A(I,J) = 0; B(I,J) = max([um up]);
        LF(I,J,1) = 0; LF(I,J,2) = B(i,j);
    elseif ((min([um up])>=0)&&(max([um up])>1))
        %--one intersection (a_i,j)
        A(I,J) = min([um up]); B(I,J) = 1;
        LF(I,J,1) = A(i,j); LF(I,J,2) = 1;
    elseif ((min([um up])>=0)&&(max([um up])<=1))
        %--two intersections
        A(I,J) = min([um up]); B(I,J) = max([um up]);
        LF(I,J,1) = A(I,J); LF(I,J,2) = B(I,J);
    else
        error('Unexpected case in frechet_decide at LF (top right cell)');
    end
end
if dq<0.0,  %--  c_ij, d_ij, BF_ij
    %--line and circle do not intersect
    C(I,J) = NaN; D(I,J) = NaN;
    BF(I,J,1) = NaN; BF(I,J,2) = NaN;
else
    up = (-bq+sqrt(dq))/(2*aq); um = (-bq-sqrt(dq))/(2*aq);
    if (((up<0)&&(um<0))||((up>1)&&(um>1)))
        %--line segment outside circle
        C(I,J) = NaN; D(I,J) = NaN;
        BF(I,J,1) = NaN; BF(I,J,2) = NaN;
    elseif ((min([um up])<0)&&(max([um up])>1))
        %--line segment is interior to circle
        C(I,J) = 0; D(I,J) = 1;
        BF(I,J,1) = 0; BF(I,J,2) = 1;
    elseif ((min([um up])<=0)&&(max([um up])<=1))
        %--one intersection (d_i,j)
        C(I,J) = 0; D(I,J) = max([um up]);
        BF(I,J,1) = 0; BF(I,J,2) = D(I,J);
    elseif ((min([um up])>=0)&&(max([um up])>1))
        %--one intersection (a_i,j)
        C(I,J) = min([um up]); D(I,J) = 1;
        BF(I,J,1) = C(I,J); BF(I,J,2) = 1;
    elseif ((min([um up])>=0)&&(max([um up])<=1))
        %--two intersections
        C(I,J) = min([um up]); D(I,J) = max([um up]);
        BF(I,J,1) = C(I,J); BF(I,J,2) = D(I,J);
    else
        error('Unexpected case in frechet_decide at BF (top right cell)');
    end
end
%--plot the free space diagram
if plotFSD==1  
    figure(2)
%     figure('Name','FSD','NumberTitle','off')
    if exist('h')==1
        close(2)
        figure(2)
%         figure('Name','Graph','NumberTitle','off')
        h = [];
    end
    ih = 1;
    for i=1:I-1
        for j=1:J-1
            x = [j-1 j j j-1];  y = [i-1 i-1 i i];
            h(ih)=patch(x',y','k'); ih = ih + 1;
            x = []; y = [];
            if ~isnan(C(i,j))
                x = [x j-1+C(i,j)];  y = [y i-1];
            end
            if ~isnan(D(i,j))
                x = [x j-1+D(i,j)];  y = [y i-1];
            end
            if ~isnan(A(i,j+1))
                x = [x j];  y = [y i-1+A(i,j+1)];
            end
            if ~isnan(B(i,j+1))
                x = [x j];  y = [y i-1+B(i,j+1)];
            end
            if ~isnan(D(i+1,j))
                x = [x j-1+D(i+1,j)];  y = [y i];
            end
            if ~isnan(C(i+1,j))
                x = [x j-1+C(i+1,j)];  y = [y i];
            end
            if ~isnan(B(i,j))
                x = [x j-1];  y = [y i-1+B(i,j)];
            end
            if ~isnan(A(i,j))
                x = [x j-1];  y = [y i-1+A(i,j)];
            end
            if length(x)>2
                h(ih) = patch(x',y','w'); ih = ih + 1;
            elseif length(x)>0
                x
                y
                fprintf('patch error\n');
            end
        end
    end
end
%--Compute the reachable sets for each cell
for i=2:I    %--fill in column 1
    LR(i,1,1) = NaN; LR(i,1,2) = NaN;
end
for j=2:J    %--fill in row 1
    BR(1,j,1) = NaN; BR(1,j,2) = NaN;
end
for i=1:I-1
    for j=1:J-1
        if (i==1)&&(j==1)
            if (LF(i,j,1)==0)&&(BF(i,j,1)==0)   %--start at the origin
                LR(i,j+1,1) = LF(i,j+1,1); LR(i,j+1,2) = LF(i,j+1,2);
                BR(j+1,i,1) = BF(j+1,i,1); BR(j+1,i,2) = BF(j+1,i,2);
            else
                LR(i,j+1,1) = NaN; LR(i,j+1,2) = NaN;
                BR(j+1,i,1) = NaN; BR(j+1,i,2) = NaN;
            end
        else
            if isnan(LR(i,j,1))&&isnan(BR(i,j,1))
                LR(i,j+1,1) = NaN; LR(i,j+1,2) = NaN;
                BR(i+1,j,1) = NaN; BR(i+1,j,2) = NaN;
            elseif (~isnan(LR(i,j,1)))&&isnan(BR(i,j,1))
                BR(i+1,j,1) = BF(i+1,j,1); BR(i+1,j,2) = BF(i+1,j,2);
                if (LF(i,j+1,2)<LR(i,j,1))||isnan(LF(i,j+1,2))
                    LR(i,j+1,1) = NaN; LR(i,j+1,2) = NaN;
                elseif LF(i,j+1,1)>LR(i,j,1)
                    LR(i,j+1,1) = LF(i,j+1,1);
                    LR(i,j+1,2) = LF(i,j+1,2);
                else
                    LR(i,j+1,1) = LR(i,j,1);
                    LR(i,j+1,2) = LF(i,j+1,2);
                end
            elseif isnan(LR(i,j,1))&&(~isnan(BR(i,j,1)))
                LR(i,j+1,1) = LF(i,j+1,1); LR(i,j+1,2) = LF(i,j+1,2);
                if (BF(i+1,j,2)<BR(i,j,1))||isnan(BF(i+1,j,2))
                    BR(i+1,j,1) = NaN; BR(i+1,j,2) = NaN;
                elseif BF(i+1,j,1)>BR(i,j,1)
                    BR(i+1,j,1) = BF(i+1,j,1);
                    BR(i+1,j,2) = BF(i+1,j,2);
                else
                    BR(i+1,j,1) = BR(i,j,1);
                    BR(i+1,j,2) = BF(i+1,j,2);
                end
            else
                LR(i,j+1,1) = LF(i,j+1,1); LR(i,j+1,2) = LF(i,j+1,2);
                BR(i+1,j,1) = BF(i+1,j,1); BR(i+1,j,2) = BF(i+1,j,2);
            end
        end
        if printFSD==1
            fprintf('cell(%d,%d):\n',i,j);
            fprintf('\tLF(i,j)=[%f,%f] ',LF(i,j,1),LF(i,j,2));
            fprintf('\tBF(i,j)=[%f,%f] \n',BF(i,j,1),BF(i,j,2));
            fprintf('\tLR(i,j)=[%f,%f] ',LR(i,j,1),LR(i,j,2));
            fprintf('\tBR(i,j)=[%f,%f] \n',BR(i,j,1),BR(i,j,2));
            fprintf('\tLR(i,j+1)=[%f,%f] ',LR(i,j+1,1),LR(i,j+1,2));
            fprintf('\tBR(i+1,j)=[%f,%f] \n',BR(i+1,j,1),BR(i+1,j,2));
        end
    end
end
%--decide
if (BR(I,J-1,2)==1)||(LR(I-1,J,2)==1)
    decide = 1;
else
    decide = 0;
end
title(['Free Space for Leash Length = ' num2str(len)]);
return
end