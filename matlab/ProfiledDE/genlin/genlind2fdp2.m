function r = genlind2fdp2(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlindfdp
%
% The derivative of a forced linear differential equations with respect to
% p at times t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(fd_cell,2);
m = length(p);

r = cell(n,m,m);
r(:) = {0};

% if(~isempty(more.force))
%     if(~isempty(more.dforcedp))
%         for(i = 1:length(more.dforcedp_sub)
%             dfs = dforcedp{i}(t,p,more.force_input);
%             r{more.force_sub(more.dforcedp_sub(i,1)),...
%                 more.force_sub(more.dforcedp_sub(i,1)),more.dforcedp_sub(i,2)} = dfs;
%                      r{more.force_sub(more.dforcedp_sub(i,1)),...
%                 more.dforcedp_sub(i,2),more.force_sub(more.dforcedp_sub(i,1))} = dfs;
%         end
%        
%         if(~isempty(more.d2forcedp2))
%             b = more.force_vec;
%             b(more.force_sub) = p((size(more.sub)+1):(size(more.sub)+1):length(more.force_sub));
%             
%             for(i = 1:length(more.d2forcedp2))
%                 ddfs = d2forcedp2{i}(t,p,more.force_input);
%                 r{more.force_sub(more.d2forcedp2_sub(i,1)),...
%                     more.d2forcedp2_sub(i,2),more.d2forcedp2_sub(i,3)} =...
%                     b(more.d2forcedp2_sub(i,1))*dfs;
%             end
% 
%         end
%     end
% end


end