function [f,J] = exact_lin_est(pars,Y,T,obs)

A = reshape(pars(1:4),2,2)';
c = pars(5:6);

E = zeros(2,2,length(T));

for i = 1:length(T)
    E(:,:,i) = expm(A*T(i));
end

X = [squeeze(E(1,:,obs{1}))'; squeeze(E(2,:,obs{2}))'];

f = Y - X*c;

if nargout > 1
   J = [ (squeeze(E(1,:,obs{1}))'*c).*T(obs{1}) (squeeze(E(2,:,obs{1}))'*c).*T(obs{1})  zeros(length(obs{1}),2);
          zeros(length(obs{2}),2) (squeeze(E(1,:,obs{2}))'*c).*T(obs{2}) (squeeze(E(2,:,obs{2}))'*c).*T(obs{2})];

   J = -[J X];
        
end