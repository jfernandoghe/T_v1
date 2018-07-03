function master_comp ()
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
pung=5;
ne=3;
npe=pung-2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Random Walk Parameters
% pung      = number of points to generate
% ne        = number of errors per micro-saccade
% npe       = number of points with micro-saccades
    [Sal1, Sal2]=RandomWalk2(50, 50, pung, ne, npe, maxFrechet);
    qx =Sal1(1,:)
    qy =Sal1(2,:)
    qpx=Sal2(1,:);
    qpy=Sal2(2,:);   
    [Vsal1, strs1,strs2]=edit_dist(qx, qy, 120); 
    [Vsal2, strp1,strp2]=edit_dist(qpx, qpy, 120);
    Vsal1
    Vsal2
    [v1]=EditDistance(Vsal1,Vsal2)
    lev(Vsal1,Vsal2)
end

function [vsal,str1, str2] = edit_dist (rawx, rawy, inte)
str1=''; str2=''; vsal=''; band=1;
bx=(0:inte:1920);
by=(0:inte:1080);
valx=('a':'z');
valy=('A':'Z');
if size(rawx,2)==size(rawy,2)
    for i=1:size(rawx,2)
        for j=1:size(bx,2)
            if rawx(i)<bx(j)
                str1=horzcat(str1,valx(j));
                for k=1:size(by,2)
                    if rawy(i)<by(k)
                        str2=horzcat(str2,valy(k));
%                         vsal(band)=(i*j);
                        vsal(band)=char(33+(i*j));
                        band=band+1;
                        break;
                    end
                end
                break;
            end
        end
    end
end
% for m=1:size(rawy,2)
%     for n=1:size(by,2)
%         if rawy(m)<by(n)
%             str2=horzcat(str2,valy(n));
%             break;
%         end
%     end
% end
% for o=0:size(by,2)
%     for p=size(bx,2)
%         
%     end
% end
end

function [V,v] = EditDistance(string1,string2)
% Edit Distance is a standard Dynamic Programming problem. Given two 
% strings s1 and s2, the edit distance between s1 and s2 is the minimum 
% number of operations required to convert string s1 to s2. The following 
% operations are typically used:
% 
% Replacing one character of string by another character.
% Deleting a character from string
% Adding a character to string
% Example:
% s1='article'
% s2='ardipo'
% EditDistance(s1,s2)
% > 4
% you need to do 4 actions to convert s1 to s2
% replace(t,d) , replace(c,p) , replace(l,o) , delete(e)
% using the other output, you can see the matrix solution to this problem
%
% by : 
% Reza Ahmadzadeh (seyedreza_ahmadzadeh@yahoo.com - reza.ahmadzadeh@iit.it)
% 14-11-2012

m=length(string1);
n=length(string2);
v=zeros(m+1,n+1);
for i=1:1:m
    v(i+1,1)=i;
end
for j=1:1:n
    v(1,j+1)=j;
end
for i=1:m
    for j=1:n
        if (string1(i) == string2(j))
            v(i+1,j+1)=v(i,j);
        else
            v(i+1,j+1)=1+min(min(v(i+1,j),v(i,j+1)),v(i,j));
        end
    end
end
V=v(m+1,n+1);
end

function d = lev(s,t)
% Levenshtein distance between strings or char arrays.
% lev(s,t) is the number of deletions, insertions,
% or substitutions required to transform s to t.
% https://en.wikipedia.org/wiki/Levenshtein_distance

    s = char(s);
    t = char(t);
    m = length(s);
    n = length(t);
    x = 0:n;
    y = zeros(1,n+1);   
    for i = 1:m
        y(1) = i;
        for j = 1:n
            c = (s(i) ~= t(j)); % c = 0 if chars match, 1 if not.
            y(j+1) = min([y(j) + 1
                          x(j+1) + 1
                          x(j) + c]);
        end
        % swap
        [x,y] = deal(y,x);
    end
    d = x(n+1);
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