
clear;
load('data/synthetic/synthetic.mat')
load('data/synthetic/constraintssyn.mat')
G=W;
C=cell2mat(CM(1,10));
rho=1;
beta=9;
gamma=0.01;
result=mvcsc(G,C,gamma,beta,rho,Y);

