function [Q] = lLRLSC(X,para,k)

[d,n]=size(X);
[l,n]=size(para.H);

 Y1=zeros(d,n);  
 Y2=zeros(n,n);   
 Y3=zeros(l,n);

Q = para.Z;%快对角正则引出来的变量


one = ones(n,1);
L = diag(para.Z*one)-para.Z;

iter = 0;
maxIter = 1000; 

mu = 1e-8;
pho = 1.1;
max_mu = 1e6;
while iter < maxIter
    iter = iter + 1;  
    
         % update G
    A = mu*((para.P)'*para.P); B = (eye(n)+mu*eye(n)-(para.S)'-para.S+para.S*(para.S)')+eye(n)*1e-4;
    C = mu*(para.P)'*X+(para.P)'*Y1+mu*para.H+Y3;
    para.G = lyap(A,B,C); 
    
          %update S
          para.S=inv((para.G)'*para.G+mu*eye(n))*((para.G)'*para.G+mu*para.Z+Y2);
     
    
           %update H
    para.H  = solve_l1l2((para.G-Y3/mu)',para.lamda1/mu);
    para.H = (para.H)';
%     para.H = softth((para.G-Y3/mu),para.lamda1/mu);
    
       %update Z   
    para.Z = para.S-Y2/mu-para.lamda2/mu*(repmat(diag(Q),1,n)-Q);
    para.Z  = max(0,(para.Z +(para.Z)')/2);
    para.Z  = para.Z -diag(diag(para.Z));
    L = diag(para.Z*one)-para.Z;    
    
      % update Q
    [V, D] = eig(L);
    D = diag(D);
    [~, ind] = sort(D);    
    Q = V(:,ind(1:k))*V(:,ind(1:k))';
    
         % update P
       PP=para.G*(X+Y1/mu)';
       [UU,SS,VV] = svd (PP,'econ'); 
       PT = UU*VV';
       para.P = PT';
       
  
      %更新Ei
%     GG1 = X-para.W1*para.X1+Y11/mu;
%     [para.E1] = solve_l1l2(GG1,para.afa/mu); 
    

%        
     % updata multipliers
    Y1 = Y1+ mu*(X-para.P*para.H); 
    Y2 = Y2+ mu*(para.Z-para.S);
    Y3 = Y3+ mu*(para.H-para.G);
     

%     
     mu = min(pho*mu, max_mu);
    
   thrsh = 0.00000001;
    if(norm(X-para.P*para.H,inf)<thrsh && norm(para.Z-para.S,inf)<thrsh && norm(para.H-para.G,inf)<thrsh)
        break;
    end
    

end


end