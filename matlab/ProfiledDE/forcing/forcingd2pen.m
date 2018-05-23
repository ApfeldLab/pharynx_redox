function r = forcingd2pen(coefs,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingd2pen
%
% second derivative of forcing penalty wrt coefs
%
% Note that unlike forcingpen and forcingdpen this does not output values
% at quadrature points, but rather provides a full Hessian. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


quads = getquadvals(more.basis{1});

pens = cell(length(more.basis),length(more.deg));

for i = 1:length(more.deg)
    pens(:,i) = eval_basis_cell(quads(:,1),more.basis,more.deg(i));
end

lambda = more.lambda;

if numel(lambda) == 1 
    lambda = lambda*ones(length(more.basis),length(more.deg)); 
 elseif numel(lambda) == length(more.deg)
    lambda = reshape(lambda,1,length(more.deg));
    lambda = ones(length(more.basis),1)*lambda;
elseif numel(lambda) == length(more.basis)
    lambda = reshape(lambda,length(more.basis),1);
    lambda = lambda*ones(1,length(more.deg));
else
    error('lambda penalty of unintelligible dimension');
end


for i = 1:length(more.basis)
    for j = 1:length(more.deg)
        pens{i,j} = lambda(i,j)*pens{i,j}'*diag(quads(:,2))*pens{i,j};
    end
end

r = 0;

for i = 1:length(more.deg)
    r = r + mattdiag_cell(pens,0);
end