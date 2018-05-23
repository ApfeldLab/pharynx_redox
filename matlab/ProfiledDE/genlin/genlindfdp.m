function r = genlindfdp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlindfdp
%
% The derivative of linear differential equations with respect to p 
% at times t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    more = [];
end

n = size(fd_cell,2);

more = checkmore(more,n);

y = eval_fdcell(t,fd_cell,0);

r = cell(n,length(p));
r(:) = {0};

for i = 1:size(more.sub,1)
    r{more.sub(i,1),i} = y{more.sub(i,2)};    
end

if ~isempty(more.force)  
    m = size(more.sub,1);
    fs = zeros(length(t),length(more.force));
    for i = 1:length(more.force)
        if isa_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force_input);
        end
        
%        r{more.force_sub(i),i+m} = r{more.force_sub(i),i+m} + fs;
    end

    for i = 1:size(more.force_sub,1)
       r{more.force_sub(i,1),i+m} =  r{more.force_sub(i,1),i+m} + ...
           fs(:,more.force_sub(i,2));
    end
    
%     if ~isempty(more.dforcedp)) %%% BALLS
%         
%         b = more.force_mat;
%         b(more.force_sub) = p((size(more.sub,1)+1):(size(more.sub,1)+length(more.force_sub)));
%         
%         for i = 1:size(more.dforcedp,1))             
%             dfs = dforcedp{i}(t,p,more.force_input);
%             r{more.dforcedp_sub(i,:)} = r{more.dforcedp_sub(i,:)} + ...
%                 b(more.dforcedp_sub(i,1))*dfs;
%     end
end


end