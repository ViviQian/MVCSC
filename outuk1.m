function uk=outuk1(v,gama,M)%optimal w
[v,index]=sort(v);
v(M+1,1)=v(M,1);
for p=1:M  %Çó³öp
    theta=(gama+sum(v(1:p,1)))/min(p,M);
    for i=1:M+1
        if theta-v(i,1)<=0
            break
        end
    end
    if p==i-1 
        break
    end
end
W=zeros(M,1);%Wm
for i=1:M
    if i<=p
        W(i,1)=(theta-v(i,1))/(gama);
    end
end
uk=W(index,1);