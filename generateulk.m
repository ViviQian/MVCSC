function ulk=generateulk(Lk,N,M,uk)
ulk=zeros(N,N);

for j=1:M
L=cell2mat(Lk(1,j));
ulk=ulk+uk(j,1)*L;
end

end