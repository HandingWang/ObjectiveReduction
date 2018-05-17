%%%%    Authors:    Handing Wang,  Xin Yao
%%%%    Xidian University, China, and University of Birmingham, UK
%%%%    EMAIL:      wanghanding.patch@gmail.com, X.Yao@cs.bham.ac.uk
%%%%    WEBSITE:    http://www.cs.bham.ac.uk/~xin/
%%%%    DATE:       March 2015
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang,  Xin Yao, Objective Reduction Based on Nonlinear Correlation Information Entropy, Soft Computing, Accepted.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
%You can run Two_Arch2 by calling function [ Sr ] = Select_Ojbectives( obj)
%to obtain reducted objective subset. You can use this function in every
%generation of Pareto-based MOEAs to reduce the number of objectives.


function [ Sr ] = Select_Ojbectives( obj )
% Usage: [ Sr ] = Select_Ojbectives( obj )
%
% Input:
% obj           - objective values of the non-dominated set
%
% Output: 
% Sr            - reducted objective subset(index)
% 
 c=size(obj,2);%number of objectives
 St=1:c;
 Sr=[];
 R=RNCC(obj);%NCIE matrix
 i=0;
 while size(St,2)~=0
     t=R<0;
     s=sum(t);
     [sx,index]=max(s);
     if sx==0
         t=sum(R);
         [sx,index]=max(t);
     end
         i=i+1;
         Sr(i)=St(index);
         t=R([1:index-1,index+1:size(R,2)],index);
         St(index)=[];
         R(:,index)=[];  
         R(index,:)=[];

        if max(t)>0
            I1=find(t>=0);
            n=size(I1,1);
            tt=t(I1);
            [ts,I2]=sort(tt);
            ts2=[0;ts(1:n-1)];
% objective classfication
            [tmax,It]=max(ts-ts2);
            index2=I1(I2(It:n));
        else
            index2=[];
        end
         St(index2)=[];
         R(:,index2)=[];  
         R(index2,:)=[];  
 end
end

function [ R ] = RNCC(r)
% Usage: [ R ] = RNCC(r)
%
% Input:
% r             - dataset
%
% Output: 
% R             - NCIE matrix
% 
c=size(r,2);
R=eye(c);
for i=1:c
    for j=i+1:c
        R(i,j)=NCC1(r(:,i),r(:,j));
        R(j,i)=R(i,j);
    end
end
end

function [ ncc ] = NCC1( x,y)
% Usage: [ ncc ] = NCC1( x,y)
%
% Input:
% x              - variable
% y              - variable
%
% Output: 
% ncc            - NCIE of x and y
% 
n=size(x,1);
b=fix(n^0.5);
if max(x)~=min(x) & max(y)~=min(y)
    detax=(max(x)-min(x)+0.00001*(max(x)-min(x)))/b;
    detay=(max(y)-min(y)+0.00001*(max(y)-min(y)))/b;
    if detax~=0& detay~=0
    p=zeros(b,b);
    x1=ceil((x-min(x)+0.000005*(max(x)-min(x)))/detax);
    y1=ceil((y-min(y)+0.000005*(max(y)-min(y)))/detay);
    x1(find(x1<=0))=1;
    y1(find(y1<=0))=1;
    x1(find(x1>b))=b;
    y1(find(y1>b))=b;
    for i=1:n 
        p(x1(i),y1(i))=p(x1(i),y1(i))+1/n;  
    end
    ncc=0;
    for i=1:b
        for j=1:b
            if p(i,j)~=0
            ncc=ncc+(p(i,j))*log(p(i,j))/log(b);
            end
        end
    end
    for i=1:b
        if sum(p(i,:))~=0
            ncc=ncc-(sum(p(i,:)))*log(sum(p(i,:)))/log(b);
        end
    end
    for i=1:b
        if sum(p(:,i))~=0
            ncc=ncc-(sum(p(:,i)))*log(sum(p(:,i)))/log(b);
        end
    end
    else
        ncc=0;
    end
else
    ncc=0;
end
cor=cov(x,y);
if cor(1,2)<0
    ncc=0-ncc;
end
end


