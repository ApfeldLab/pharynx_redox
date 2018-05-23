%% Find cellular regions and take the mean (builds the stacked per Avg_perCell matrix) 
%The input of this script are Data_Lev or Data_TB matrices (output of Linearized Vectors script)

%% PARAMETERS from user

% ENTER VALUE

%Data=[Data_HD233_2do_120430; Data_HD240_2do_120430; Data_HD244_2do_120430; Data_HD236_2do_120430]; % by now do manual concatenation of the different genotypes

nframes= size(Data,1)/100;

start_medax=2;        %start of cell boundary
end_medax =93;         %end of cell boundary

start_pm7=3;        %start of cell boundary
end_pm7 =10;         %end of cell boundary

start_pm5=21;        %start of cell boundary
end_pm5 =40;         %end of cell boundary

start_pm4=53;        %start of cell boundary
end_pm4 =60;         %end of cell boundary

start_pm3=71;        %start of cell boundary
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
    end    
        avg_pm7(:,8)= 1000*(avg_pm7(:,6)./ avg_pm7(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm7(:,9)= log(avg_pm7(:,6))- log(avg_pm7(:,7)); %col 11 is ln(R)


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
    end    
        avg_pm5(:,8)= 1000*(avg_pm5(:,6)./ avg_pm5(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm5(:,9)= log(avg_pm5(:,6))- log(avg_pm5(:,7)); %col 11 is ln(R)

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
    end    
        avg_pm4(:,8)= 1000*(avg_pm4(:,6)./ avg_pm4(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm4(:,9)= log(avg_pm4(:,6))- log(avg_pm4(:,7)); %col 11 is ln(R)



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
    end    
        avg_pm3(:,8)= 1000*(avg_pm3(:,6)./ avg_pm3(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm3(:,9)= log(avg_pm3(:,6))- log(avg_pm3(:,7)); %col 11 is ln(R)


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
    end    
        avg_medax(:,8)= 1000*(avg_medax(:,6)./ avg_medax(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_medax(:,9)= log(avg_medax(:,6))- log(avg_medax(:,7)); %col 11 is ln(R)

%% Concatenate Per cell matrices

  %  Avg_perCell= [avg_pm7; avg_pm5; avg_pm3; avg_medax];
  Avg_perCell= [avg_pm7; avg_pm5; avg_pm4; avg_pm3; avg_medax];