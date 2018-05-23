function r = genlinfunode(t,y,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlinfunode
%
% forced linear differential equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    more = [];
end

n = size(y,1);

more = checkmore(more,n);

pmat = more.mat;

pmat(sub2ind([n n],more.sub(:,1),more.sub(:,2))) = p(1:size(more.sub,1));

r = pmat*y;

if ~isempty(more.force)

    fs = zeros(length(t),length(more.force));
        
    for i = 1:length(more.force)
        if isa_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force_input);
        end
    end
    
    b = more.force_mat;
    b(sub2ind(size(b),more.force_sub(:,1),more.force_sub(:,2))) = ...
              p((size(more.sub,1)+1):(size(more.sub,1)+size(more.force_sub,1)));
  
    r = r + b*fs';        
end
    
end