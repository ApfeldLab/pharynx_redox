function [logrank, chi_square, df, v, sum_V] = j_logrank2(strain,test_type)
%function [logrank, chi_square, df, v, sum_V] = j_logrank(strain)
%stain is a cell array of strains/conditions;
%test type: 1 logrank (no weight)
%           2 Gehan’s Wilcoxon (weight = n at risk at i'th time)
%conditions that are empty are removed (thus not affecting df)
%Each cell in strain contains a table: col1 is death time, col2 is censoring variable (0 = not censored)
%Output is [logrank, chi_square, degrees of freedom, v, sum_V]
%v = sum(D) - sum(E), the sum over all timesteps of the observed deaths minus the expected deaths
%V is the Variance-Covariance Matrix, sum_V is the sum of V over all timesteps


%remove empties
clean{1} = [];
counter = 0;
for n = 1:size(strain,2)
    if isempty(strain{n}) ~= 1
        counter = counter + 1;
        clean{counter} = strain{n};
    end
end

strain = clean;
strain_number = size(strain,2);

% add censoring col = 0 if it is missing
for n = 1:strain_number
    if size(strain{n},2) < 2
        strain{n}(:,2) = 0;
    end
end

all_strains = [];
for n =1:strain_number
    all_strains = cat(1,all_strains,strain{n});
end

all_event = j_count(all_strains(:,1));
time_range = all_event(:,1);

% dead|censored{n} are tables: 'day' 'number of events'
death = [];
censored = [];

for n = 1:strain_number
    death{n}= j_count(strain{n}(find(strain{n}(:,2) == 0),1));
    
    if isempty(find(strain{n}(:,2) ~= 0)) == 1
        censored{n} = [];
    else
        censored{n}= j_count(strain{n}(find(strain{n}(:,2) ~= 0),1));
    end
end

% Build Censoring and Death matrices
C = zeros(size(time_range,1),strain_number);
D = zeros(size(time_range,1),strain_number);


for n = 1:strain_number
    if isempty(censored{n}) == 0
        for m = 1:size(censored{n}(:,1),1)
            index = find(time_range == censored{n}(m,1));
            C(index,n) = censored{n}(m,2);
        end
    end
    for m = 1:size(death{n}(:,1),1)
        index = find(time_range == death{n}(m,1));
        D(index,n) = death{n}(m,2);
    end
end

% Build At risk matrix and Strain_start_size matrix
A = zeros(size(time_range,1),strain_number);
Strain_start_size = zeros(1,strain_number);
for n = 1:strain_number
    for m = 1:size(time_range,1)
        if m == 1
            Strain_start_size(n) = size(strain{n},1);
            A(m,n) =  Strain_start_size(n);
        else
            A(m,n) = (A(m-1,n) - C(m-1,n) - D(m-1,n));
        end
    end
end




% Build Expected deaths matrix
E = zeros(size(time_range,1),strain_number);
for m = 1:size(time_range,1)
    E(m,:) = sum(D(m,:))/sum(A(m,:)) .* A(m,:);
end

switch test_type
    case 1
        W = ones(size(time_range,1),1);
    case 2
        W = sum(A,2);
        % size(W)
        % length(W)
    case 3
        W = 
end
%{
% Approximate calculation of logrank
chi_square1 = sum((sum(E)-sum(D)).^2./sum(E))
logrank1 = 1 - chi2cdf(chi_square,strain_number-1)
%}

% The following is required for the exact calculation of the logrank test

% Build Variance-Covariance  matrix
V = zeros(strain_number, strain_number, size(time_range,1));
for m = 1:size(time_range,1)
    Am = sum(A(m,:)); % Am is the total At risk at timestep m
    if Am == 1 % Correction for the variance when sample size is 1
        Amm = 1;
    else
        Amm = Am -1;
    end
    Dm = sum(D(m,:)); % Dm is the total Dead at timestep m
    for i = 1:strain_number
        for l = 1:strain_number
            if i == l
                V(i,l,m) = A(m,i)*(Am-A(m,i))*Dm*(Am-Dm)/(Am^2*Amm);
            else
                V(i,l,m) = -A(m,i)*A(m,l)*Dm*(Am-Dm)/(Am^2*Amm);
            end
        end
    end
end


%Strain_start_size
%cat(2,time_range,A,D,C,E)
%sum(cat(2,D,E))

% Calculate outputs
% if test_type==1
%     v = sum(D)-sum(E);
%     sum_V = sum(V,3);
%     df = strain_number-1;
%     if df > 0 & trace(sum_V) > 0
%         v_range = 1:df;
%         chi_square = v(v_range) * inv(sum(V(v_range,v_range,:),3)) * v(v_range)';
%         logrank = 1 - chi2cdf(chi_square,df);
%     else
%         %disp('else')
%         chi_square = 0;
%         logrank = 1;
%         
%     end
%     
% else
    v = sum((W*ones(1,strain_number)).*(D-E)) ;

     V= (repmat(reshape(W,1,1,length(W)),[strain_number,strain_number,1]).^2).*V;
    
    sum_V = sum(V,3);
    df = strain_number-1;
    if df > 0 & trace(sum_V) > 0
        v_range = 1:df;
        chi_square = v(v_range) * inv(sum(V(v_range,v_range,:),3)) * v(v_range)';
        logrank = 1 - chi2cdf(chi_square,df);
    else
        %disp('else')
        chi_square = 0;
        logrank = 1;
        
%     end
end





























