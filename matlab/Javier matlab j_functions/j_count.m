function a = j_count(series)
%output is a list of categories total-counts each series counts

%series = [1 2 3 0; 2 0 0 0]'

%disp('start j_count')
cat_count=zeros(1,size(series,2)+2);
cat_count(1,1)=series(1,1);

for n = 1:size(series,1)
    for m = 1:size(series,2)
        index = find(cat_count(:,1) == series(n,m),1);
        %disp([series(n) index]);
        if isempty(index) == 1
            %disp('empty')
            new_ncat = size(cat_count,1)+1;
            cat_count(new_ncat,1) = series(n,m);
            cat_count(new_ncat,2) = 1;
            cat_count(new_ncat,2+m) = 1;
        else
            %disp('add')
            cat_count(index,2)=cat_count(index,2)+1;
            cat_count(index,2+m)=cat_count(index,2+m)+1;
        end
    end
end

%cat_count
a = sortrows(cat_count,1);


        