function r = genlinfun(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlinfun
%
% linear differential equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    more = [];
end

n = size(fd_cell,2);

more = checkmore(more,n);

y = cell2mat(eval_fdcell(t,fd_cell,0));
pmat = more.mat;


pmat(sub2ind([n n],more.sub(:,1),more.sub(:,2))) = p(1:size(more.sub,1));

r1 = y*pmat';


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
  
    r1 = r1 + fs*b';    
end
    
r = cell(1,n);

for i = 1:n
    r{i} = r1(:,i);
end

end