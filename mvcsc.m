function [r]=mvcsc(G,C,gamma,beta,rho,Y)
maxiter=10;
M = size(G,2);%G is graph
N = size(G{1,1},2); % N is the number of nodes
k=length(unique(Y));% K is the number of classes,Y is the groundtruth
truth=Y;
Lk=generatelk(G,M);%L is the Lapican of G
R=cell(1,15);% save results
for j=1:10 % run times
     %初始化
    n1=size(C,1);
    lambda=zeros(n1,k);
    uk=1/M*ones(M,1);
    %z=zeros(n1,k);
    z=rand(n1,k);
    %f=-1+2*rand(N,k);  
    f=rand(N,k);
    f = MGramSchmidt(f);
for i=1:maxiter
    v=zeros(M,1);
    for a=1:M %生成v
        v(a,1)=trace(f'*cell2mat(Lk(1,a))*f);  
    end
    loss1=uk'*v;
    loss2=gamma*norm(z,1);
    loss3=beta/2*norm(uk,2);
    loss4=loss1+loss2+loss3;
    Loss(i,:)=[loss1, loss2, loss3,loss4];
    f_old=f;
    %update f
    uL=generateulk(Lk,N,M,uk);
    [f,out,tsolve]=ocmvkway(f,uL,C,z,lambda,uk,Lk,rho);
    %update z
    z=1/rho*soft(lambda+rho*C*f,gamma);
    for a=1:M %生成v
        v(a,1)=trace(f'*cell2mat(Lk(1,a))*f);  
    end
    uk=outuk1(v,beta,M);
    lambda=lambda+rho*(C*f-z);  
    if norm(f-f_old)<1e-16
       break;
    end  
end 

for n=1:40
    y2=kmeans(f,k);
    [AR,RI,~,~]=RandIndex(y2,truth);
    NMI=nmi(y2, truth);
    [ACC,MIhat,Purity] = ClusteringMeasure(y2,truth);
    %[f1,p,r] = compute_f(y2,truth);
    result1(n,:,j)=[AR NMI ACC];
end

%max(result1)

 R{1,1}=result1;
end

b=zeros(1,3);
a=max(R{1,1});
for i=1:j
    b=[b;a(:,:,i)];
end
b(1,:)=[];
r=mean(b)
end


% for i=1:10
%     figure(i);spy(G{1,i});
% end
% F=cell(1,10);
% for i=1:10
% f=-1+2*rand(N,k); 
% F{1,i}=f;
% end