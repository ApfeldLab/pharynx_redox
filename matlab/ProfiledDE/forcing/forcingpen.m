function r = forcingpen(coefs,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingpen
%
% penalizes a forcing funciton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

quads = getquadvals(more.basis{1});

forces = Make_fdcell(coefs,more.basis);

pens = cell(length(forces),length(more.deg));

for i = 1:length(more.deg)
    pens(:,i) = eval_fdcell(quads(:,1),forces,more.deg(i));
end

lambda = more.lambda;

if numel(lambda) == 1 
    lambda = lambda*ones(length(forces),length(more.deg));
elseif numel(lambda) == length(more.deg)
    lambda = reshape(lambda,1,length(more.deg));
    lambda = ones(length(more.basis),1)*lambda;
elseif numel(lambda) == length(more.basis)
    lambda = reshape(lambda,length(more.basis),1);
    lambda = lambda*ones(1,length(more.deg));    
else
    error('lambda penalty of unintelligible dimension');
end


for i = 1:length(forces)
    for j = 1:length(more.deg)  
       pens{i,j} = diag(sqrt(lambda(i,j)*quads(:,2)))*pens{i,j};
    end
end

r = [];

for i = 1:length(more.deg)
    r = [r; cell2mat(pens(:,i))];
end