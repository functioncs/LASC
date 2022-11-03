clear all;close all;
addpath(genpath(cd)) ;
name = 'YaleB_32x32';
selected_class = 30;
load(name)
fea = double(fea);
nnClass = length(unique(gnd));     % The number of classes
select_sample = [];
select_gnd    = [];
for i = 1:selected_class
    idx = find(gnd == i);
    idx_sample    = fea(idx,:);
    select_sample = [select_sample;idx_sample];
    select_gnd    = [select_gnd;gnd(idx)];
end
    
fea = select_sample';
fea = fea./repmat(sqrt(sum(fea.^2)),[size(fea,1) 1]);
gnd = select_gnd;
c   = selected_class;


X = fea;
s=gnd;
% clear fea select_gnd select_sample idx
% imshow(X,[]);



para.lamda1=0.01;
para.lamda2=1;

[d,n]=size(X);
nCluster = length(unique(s));
k=nCluster;

% initial subspace learning
C=inv(X'*X+para.lamda1*eye(n))*X'*X;
C = C- diag(diag(C));

% Afa=[0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10];aaa=length(Afa);% 
% Beta=[0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10];bbb=length(Beta);
%  AAAA1=zeros(aaa,bbb);BBBB1=zeros(aaa,bbb);
% ccc1=zeros(aaa,bbb);
% ddd1=zeros(aaa,bbb);
% 
% 
% for i=1:aaa
%     for j=1:bbb
%            para.lamda1=Afa(i);
%          para.lamda2=Beta(j);

iter=3;
for t=1:iter
    
%update S
% D=[];
% for i=1:n
%     for j=1:n
%         A=[C(:,i);C(:,j)];
%         d=pdist(A,'correlation');
%         D(i,j)=d;
%     end
% end
D=1-corr(C,'type','Pearson');
S=-D/(2*para.lamda2);
 S=S-diag(diag(S));
for ic = 1:n
        idx    = 1:n;
        idx(ic) = [];
        S(ic,idx) = EProjSimplex_new(S(ic,idx));          % 
end
 %update C,  
   C=inv(X'*X+para.lamda1*eye(n))*X'*X;
end   
 
 
   S1=abs(S)+abs((S)');
   [nmi1,ACC1,f1,RI1]=clustering(S1, nCluster, s);
 
%    AAAA1(i,j)=nmi1;
%    BBBB1(i,j)=ACC1;
%    ccc1(i,j)=f1;
%    ddd1(i,j)=RI1;
%  
%      end
%    end       
         
% addpath('Ncut_9');
% Z_out = S;
% A = Z_out;
% A = A - diag(diag(A));
% A = (A+A')/2;  
% [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(A,k);
% [value,result_label] = max(NcutDiscrete,[],2);
% result = ClusteringMeasure1(gnd, result_label)
         
         
              