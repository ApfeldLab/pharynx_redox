%% Find cellular regions and take the mean (builds the stacked per Avg_perCell matrix) 
%The input of this script are Data_Lev or Data_TB matrices (output of Linearized Vectors script)

%% PARAMETERS from user

% ENTER VALUE

%Data=[Data_HD233_2do_120430; Data_HD240_2do_120430; Data_HD244_2do_120430; Data_HD236_2do_120430]; % by now do manual concatenation of the different genotypes

nframes= size(Data,1)/100;

start_medax=3;         %start of cell boundary
end_medax =97;         %end of cell boundary

start_pm7=5;          %start of cell boundary
end_pm7 =10;         %end of cell boundary

start_pm6=15;        %start of cell boundary
end_pm6 =20; 

start_pm5=30;        %start of cell boundary
end_pm5 =50;         %end of cell boundary

start_pm4=59;        %start of cell boundary
end_pm4 =65;         %end of cell boundary

start_pm3=70;        %start of cell boundary
end_pm3 =90;         %end of cell boundary

%% Find pixels that belong to pm7 and take the mean

px_pm7= end_pm7 - start_pm7 +1;     %range of pixels that belongs to pm7
    
           
    px=Data(:,5);
    ind_pm7= find(px>=start_pm7 & px<=end_pm7); % find the linear indices of the pixels that belongs to pm7
    pm7=Data(ind_pm7,:);        %select those rows and all columns  from data matrices
    avg_pm7=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm7;                       %and for each row of the avg matrix
            avg_pm7(j,1:5)= pm7((i-(px_pm7-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm7(j,6)= mean(pm7(((i-(px_pm7-1)):i),6));      %find average 410 per block of pixels
            avg_pm7(j,7)= mean(pm7(((i-(px_pm7-1)):i),7));      %find average 470 per block of pixels
            avg_pm7(j,8)= mean(pm7(((i-(px_pm7-1)):i),8));      %find average R
            avg_pm7(j,9)= mean(pm7(((i-(px_pm7-1)):i),9));      %find average OxD
            avg_pm7(j,10)= mean(pm7(((i-(px_pm7-1)):i),10));      %find average E  
    end    
      

%% pm6

px_pm6= end_pm6 - start_pm6 +1;     %range of pixels that belongs to pm6
    
           
    px=Data(:,5);
    ind_pm6= find(px>=start_pm6 & px<=end_pm6); % find the linear indices of the pixels that belongs to pm6
    pm6=Data(ind_pm6,:);        %select those rows and all columns  from data matrices
    avg_pm6=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm6;                       %and for each row of the avg matrix
            avg_pm6(j,1:5)= pm6((i-(px_pm6-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm6(j,6)= mean(pm6(((i-(px_pm6-1)):i),6));      %find average 410 per block of pixels
            avg_pm6(j,7)= mean(pm6(((i-(px_pm6-1)):i),7));      %find average 470 per block of pixels
            avg_pm6(j,8)= mean(pm6(((i-(px_pm6-1)):i),8));      %find average R
            avg_pm6(j,9)= mean(pm6(((i-(px_pm6-1)):i),9));      %find average OxD
            avg_pm6(j,10)= mean(pm6(((i-(px_pm6-1)):i),10));      %find average E  
    end    
      

%% pm5
px_pm5= end_pm5 - start_pm5 +1;     %range of pixels that belongs to pm5
    
           
    px=Data(:,5);
    ind_pm5= find(px>=start_pm5 & px<=end_pm5); % find the linear indices of the pixels that belongs to pm5
    pm5=Data(ind_pm5,:);        %select those rows and all columns  from data matrices
    avg_pm5=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm5;                       %and for each row of the avg matrix
            avg_pm5(j,1:5)= pm5((i-(px_pm5-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm5(j,6)= mean(pm5(((i-(px_pm5-1)):i),6));      %find average 410 per block of pixels
            avg_pm5(j,7)= mean(pm5(((i-(px_pm5-1)):i),7));      %find average 470 per block of pixels
            avg_pm5(j,8)= mean(pm5(((i-(px_pm5-1)):i),8));
            avg_pm5(j,9)= mean(pm5(((i-(px_pm5-1)):i),9));
            avg_pm5(j,10)= mean(pm5(((i-(px_pm5-1)):i),10));
    end    
       

%% pm4
px_pm4= end_pm4 - start_pm4 +1;     %range of pixels that belongs to pm4
    
           
    px=Data(:,5);
    ind_pm4= find(px>=start_pm4 & px<=end_pm4); % find the linear indices of the pixels that belongs to pm4
    pm4=Data(ind_pm4,:);        %select those rows and all columns  from data matrices
    avg_pm4=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm4;                       %and for each row of the avg matrix
            avg_pm4(j,1:5)= pm4((i-(px_pm4-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm4(j,6)= mean(pm4(((i-(px_pm4-1)):i),6));      %find average 410 per block of pixels
            avg_pm4(j,7)= mean(pm4(((i-(px_pm4-1)):i),7));      %find average 470 per block of pixels
            avg_pm4(j,8)= mean(pm4(((i-(px_pm4-1)):i),8));
            avg_pm4(j,9)= mean(pm4(((i-(px_pm4-1)):i),9));
            avg_pm4(j,10)= mean(pm4(((i-(px_pm4-1)):i),10));
    end    




%% pm3
px_pm3= end_pm3 - start_pm3 +1;     %range of pixels that belongs to pm3
    
           
    px=Data(:,5);
    ind_pm3= find(px>=start_pm3 & px<=end_pm3); % find the linear indices of the pixels that belongs to pm3
    pm3=Data(ind_pm3,:);        %select those rows and all columns  from data matrices
    avg_pm3=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm3;                       %and for each row of the avg matrix
            avg_pm3(j,1:5)= pm3((i-(px_pm3-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm3(j,6)= mean(pm3(((i-(px_pm3-1)):i),6));      %find average 410 per block of pixels
            avg_pm3(j,7)= mean(pm3(((i-(px_pm3-1)):i),7));      %find average 470 per block of pixels
            avg_pm3(j,8)= mean(pm3(((i-(px_pm3-1)):i),8));
            avg_pm3(j,9)= mean(pm3(((i-(px_pm3-1)):i),9));
            avg_pm3(j,10)= mean(pm3(((i-(px_pm3-1)):i),10));
    end    



%% Medial Axis
% take the average of the medial axis, excluding the ends

px_medax= end_medax - start_medax +1;     %range of pixels that belongs to medax
    
           
    px=Data(:,5);
    ind_medax= find(px>=start_medax & px<=end_medax); % find the linear indices of the pixels that belongs to medax
    medax=Data(ind_medax,:);        %select those rows and all columns  from data matrices
    avg_medax=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_medax;                       %and for each row of the avg matrix
            avg_medax(j,1:5)= medax((i-(px_medax-1)),1:5); %fill in metadata (col 1 to 7)
            avg_medax(j,6)= mean(medax(((i-(px_medax-1)):i),6));      %find average 410 per block of pixels
            avg_medax(j,7)= mean(medax(((i-(px_medax-1)):i),7));      %find average 470 per block of pixels
            avg_medax(j,8)= mean(medax(((i-(px_medax-1)):i),8));
            avg_medax(j,9)= mean(medax(((i-(px_medax-1)):i),9));
            avg_medax(j,10)= mean(medax(((i-(px_medax-1)):i),10));
    end    


%% Concatenate Per cell matrices

  %  Avg_perCell= [avg_pm7; avg_pm5; avg_pm3; avg_medax];
  Avg_perCell= [avg_medax; avg_pm7; avg_pm6; avg_pm5; avg_pm4; avg_pm3];