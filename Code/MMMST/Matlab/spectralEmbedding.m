function [m1,m2]=spectralEmbedding(A1,A2,k,beta)
%%
% A1 is one domain
% A2 is another domain
% beta is the parameter
% k is the number of dimentions wanted

[r1,c1]=size(A1);
[r2,c2]=size(A2);

if size(A2,1)<size(A1,1)
    number=size(A1,1)-size(A2,1);
    while number>0
        A2=[A2;A2(1,:)];
        number=number-1;
    end
end

if size(A2,1)>size(A1,1)
    number=size(A2,1)-size(A1,1);
    while number>0
        A1=[A1;A1(1,:)];
        number=number-1;
    end
end

A1T=A1*A1';
A2T=A2*A2';


A1=2*A1T+beta^2/2*A2T;
A2=beta*(A1T+A2T);
A3=A2';
A4=beta^2/2*A1T+2*A2T;


A=[A1,A2;A3,A4];
%I=lambda*[eye(size(A1T,1)),zeros(size(A1T,1),size(A1T,1));zeros(size(A1T,1),size(A1T,1)),-eye(size(A1T,1))];
ATA=A;
[V,D,V2] = svd(ATA);

RV=V(:,1:k);

[row1,col1]=size(A1);
[row2,col2]=size(A2);

m1=RV(1:r1,1:k);
m2=RV(max(r1,r2)+1:max(r1,r2)+r2,1:k);

%% Normalization
for j=1:size(m1,2)
    maxN=0;
    for i=1:size(m1,1)
        if abs(m1(i,j))>maxN
            maxN=abs(m1(i,j));
        end
        %sum=sum+m1(i,j);
    end
    for i=1:size(m1,1)
        m1(i,j)=m1(i,j)/maxN;
    end
end

for j=1:size(m2,2)
    maxN=0;
    for i=1:size(m2,1)
        if abs(m2(i,j))>maxN
            maxN=abs(m2(i,j));
        end
        %sum=sum+m2(i,j);
    end
    for i=1:size(m2,1)
        m2(i,j)=m2(i,j)/maxN;
    end
end
        

