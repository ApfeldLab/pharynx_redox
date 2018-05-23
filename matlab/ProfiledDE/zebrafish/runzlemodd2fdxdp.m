function y = runzlemodd2fdxdp(t,DEfd,p)

r = eval_fdcell(DEfd);
y = cell(2,2,18);

y(:,1,1) = {0 0};
y(:,1,2) = {1 0};
y(:,1,3) = {2 .* r{1} 0};
y(:,1,4) = {3 .* r{1} .^ 2 0};
y(:,1,5) = {0 0};
y(:,1,6) = {r{2} 0};

y(:,1,7) = {0 0};
y(:,1,8) = {0 1};
y(:,1,9) = {0 0};

y(:,1,10) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) 0};
y(:,1,11) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* r{1} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) 0};
y(:,1,12) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* r{1} .^ 2 + 0.2e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* r{1} 0};
y(:,1,13) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* r{1} .^ 3 + 0.3e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* r{1} .^ 2 0};
y(:,1,14) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* r{2} 0};
y(:,1,15) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* r{1} .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* r{2} 0};

y(:,1,16) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(18) .^ 2 - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .^ 2 .* r{1} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .^ 2 .* r{1} .^ 2 + 0.2e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .* r{1} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .^ 2 .* r{1} .^ 3 + 0.3e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .* r{1} .^ 2 - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) .^ 2 .* r{2} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .^ 2 .* r{1} .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .* r{2} 0};
y(:,1,17) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(18) .^ 2 .* r{2} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .^ 2 .* r{1} .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .* r{2} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .^ 2 .* r{1} .^ 2 .* r{2} + 0.2e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .* r{2} .* r{1} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .^ 2 .* r{1} .^ 3 .* r{2} + 0.3e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .* r{2} .* r{1} .^ 2 - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) .^ 2 .* r{2} .^ 2 - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .^ 2 .* r{1} .* r{2} .^ 2 + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .* r{2} .^ 2 0};
y(:,1,18) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(18) .* (r{1} - p(16) - p(17) .* r{2}) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .* r{1} .* (r{1} - p(16) - p(17) .* r{2}) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* r{1} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* (r{1} - p(16) - p(17) .* r{2}) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .* r{1} .^ 2 .* (r{1} - p(16) - p(17) .* r{2}) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* r{1} .^ 2 - 0.2e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* (r{1} - p(16) - p(17) .* r{2}) .* r{1} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .* r{1} .^ 3 .* (r{1} - p(16) - p(17) .* r{2}) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* r{1} .^ 3 - 0.3e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* (r{1} - p(16) - p(17) .* r{2}) .* r{1} .^ 2 + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) .* r{2} .* (r{1} - p(16) - p(17) .* r{2}) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .* r{1} .* r{2} .* (r{1} - p(16) - p(17) .* r{2}) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* r{1} .* r{2} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* (r{1} - p(16) - p(17) .* r{2}) .* r{2} 0};


y(:,2,1) = {0 0};
y(:,2,2) = {0 0};
y(:,2,3) = {0 0};
y(:,2,4) = {0 0};
y(:,2,5) = {1 0};
y(:,2,6) = {r{1} 0};

y(:,2,7) = {0 0};
y(:,2,8) = {0 0};
y(:,2,9) = {0 1};

y(:,2,10) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* p(17) 0};
y(:,2,11) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* p(17) .* r{1} 0};
y(:,2,12) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* p(17) .* r{1} .^ 2 0};
y(:,2,13) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* p(17) .* r{1} .^ 3 0};
y(:,2,14) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* p(17) .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) 0};
y(:,2,15) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(18) .* p(17) .* r{1} .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* r{1} 0};

y(:,2,16) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(18) .^ 2 .* p(17) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .^ 2 .* r{1} .* p(17) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .^ 2 .* r{1} .^ 2 .* p(17) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .^ 2 .* r{1} .^ 3 .* p(17) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) .^ 2 .* r{2} .* p(17) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .^ 2 .* r{1} .* r{2} .* p(17) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .* r{1} 0};
y(:,2,17) = {0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(18) .^ 2 .* p(17) .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(18) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .^ 2 .* p(17) .* r{1} .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .* r{1} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .^ 2 .* p(17) .* r{1} .^ 2 .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .* r{1} .^ 2 + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .^ 2 .* p(17) .* r{1} .^ 3 .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .* r{1} .^ 3 + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) .^ 2 .* p(17) .* r{2} .^ 2 + 0.2e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) .* r{2} + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .^ 2 .* p(17) .* r{1} .* r{2} .^ 2 + 0.2e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .* r{1} .* r{2} 0};
y(:,2,18) = {-0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(18) .* p(17) .* (r{1} - p(16) - p(17) .* r{2}) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(10) .* p(17) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(18) .* p(17) .* r{1} .* (r{1} - p(16) - p(17) .* r{2}) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(11) .* p(17) .* r{1} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(18) .* p(17) .* r{1} .^ 2 .* (r{1} - p(16) - p(17) .* r{2}) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(12) .* p(17) .* r{1} .^ 2 - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(18) .* p(17) .* r{1} .^ 3 .* (r{1} - p(16) - p(17) .* r{2}) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(13) .* p(17) .* r{1} .^ 3 - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(18) .* p(17) .* r{2} .* (r{1} - p(16) - p(17) .* r{2}) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* p(17) .* r{2} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(14) .* (r{1} - p(16) - p(17) .* r{2}) - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(18) .* p(17) .* r{1} .* r{2} .* (r{1} - p(16) - p(17) .* r{2}) + 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* p(17) .* r{1} .* r{2} - 0.1e1 ./ exp(p(18) .* (r{1} - p(16) - p(17) .* r{2})) .* p(15) .* (r{1} - p(16) - p(17) .* r{2}) .* r{1} 0};

