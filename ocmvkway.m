function [X,out,tsolve]=ocmvkway(f,uL,C,z,lambda,uk,Lk,rho)

X=f;

opts.record = 0; 
opts.mxitr  = 1000;
opts.xtol = 1e-5;
opts.gtol = 1e-5;
opts.ftol = 1e-8;

% X0 = randn(n,k); %≥ı ºªØ
% X0 = orth(X0);
tic; [X, out]= OptStiefelGBB(X, @fun, opts, uL,C,z,lambda,rho ); tsolve = toc;
out.fval = -2*out.fval; % convert the function value to the sum of eigenvalues
%fprintf('\nOptM: obj: %7.6e, itr: %d, nfe: %d, cpu: %f, norm(XT*X-I): %3.2e \n', ...
            %out.fval, out.itr, out.nfe, tsolve, norm(X'*X - eye(k), 'fro') );
end

function [F, G] = fun(X,uL,C,z,lambda,rho)
  G = uL*X+C'*lambda+rho*C'*(C*X-z);
  %F = 0.5*trace(uk(1,1)*X'*Lk{1,1}*X+uk(2,1)*X'*Lk{1,2}*X) +trace(lambda'*(C*X-z))+0.5*rho*norm(C*X-z);
  %F = 0.5*trace(uk(1,1)*X'*Lk{1,1}*X+uk(2,1)*X'*Lk{1,2}*X+uk(3,1)*X'*Lk{1,3}*X+uk(4,1)*X'*Lk{1,4}*X+uk(5,1)*X'*Lk{1,5}*X) +trace(lambda'*(C*X-z))+0.5*rho*norm(C*X-z);
  F=0.5*trace(X'*uL*X)+trace(lambda'*(C*X-z))+0.5*rho*norm(C*X-z);
end
%+uk(7,1)*X'*Lk{1,7}*X+uk(8,1)*X'*Lk{1,8}*X+uk(9,1)*X'*Lk{1,9}*X+uk(10,1)*X'*Lk{1,10}*X+uk(11,1)*X'*Lk{1,11}*X+uk(12,1)*X'*Lk{1,12}*X