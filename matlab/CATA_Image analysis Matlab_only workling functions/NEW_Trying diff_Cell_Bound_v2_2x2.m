%% Find cellular regions and take the mean (builds the stacked per Avg_perCell matrix) 
%The input of this script are Data_Lev or Data_TB matrices (output of Linearized Vectors script)

%% PARAMETERS from user

% ENTER VALUE

%Data=[Data_HD233_2do_120430; Data_HD240_2do_120430; Data_HD244_2do_120430; Data_HD236_2do_120430]; % by now do manual concatenation of the different genotypes

nframes= size(Data,1)/100;

% 14_02_24 Trying new cell boundaries

start_pm8=3;        %start of cell boundary //NEW
end_pm8 =6;         %end of cell boundary

start_pm7=7;        %start of cell boundary
end_pm7 =11;         %end of cell boundary 

start_pm6=15;        %start of cell boundary
end_pm6 =18;         %end of cell boundary

start_pm5=35;        %start of cell boundary
end_pm5 =45;         %end of cell boundary

start_pm4=58;        %start of cell boundary
end_pm4 =63;         %end of cell boundary

start_pm3=73;        %start of cell boundary
end_pm3 =92;         %end of cell boundary

start_pm3p=73;        %start of cell boundary
end_pm3p =80; 

start_pm3a=85;        %start of cell boundary
end_pm3a =92;         %end of cell boundary

start_medax=2;        %start of cell boundary
end_medax =97;         %end of cell boundary  //NEW before end_medax =95;   

  

%% pm8
px_pm8= end_pm8 - start_pm8 +1;     %range of pixels that belongs to pm8
    
           
    px=Data(:,5);
    ind_pm8= find(px>=start_pm8 & px<=end_pm8); % find the linear indices of the pixels that belongs to pm5
    pm8=Data(ind_pm8,:);        %select those rows and all columns  from data matrices
    avg_pm8=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm8;                       %and for each row of the avg matrix
            avg_pm8(j,1:5)= pm8((i-(px_pm8-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm8(j,6)= mean(pm8(((i-(px_pm8-1)):i),6));      %find average 410 per block of pixels
            avg_pm8(j,7)= mean(pm8(((i-(px_pm8-1)):i),7));      %find average 470 per block of pixels
    end    
        avg_pm8(:,8)= 1000*(avg_pm8(:,6)./ avg_pm8(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm8(:,9)= log(avg_pm8(:,6))- log(avg_pm8(:,7)); %col 11 is ln(R)

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


%% Find pixels that belong to pm7 and take the mean

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
    end    
        avg_pm6(:,8)= 1000*(avg_pm6(:,6)./ avg_pm6(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm6(:,9)= log(avg_pm6(:,6))- log(avg_pm6(:,7)); %col 11 is ln(R)


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





%% pm3p = pm3p posterior end
px_pm3p= end_pm3p - start_pm3p +1;     %range of pixels that belongs to pm3p
    
           
    px=Data(:,5);
    ind_pm3p= find(px>=start_pm3p & px<=end_pm3p); % find the linear indices of the pixels that belongs to pm3p
    pm3p=Data(ind_pm3p,:);        %select those rows and all columns  from data matrices
    avg_pm3p=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm3p;                       %and for each row of the avg matrix
            avg_pm3p(j,1:5)= pm3p((i-(px_pm3p-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm3p(j,6)= mean(pm3p(((i-(px_pm3p-1)):i),6));      %find average 410 per block of pixels
            avg_pm3p(j,7)= mean(pm3p(((i-(px_pm3p-1)):i),7));      %find average 470 per block of pixels
    end    
        avg_pm3p(:,8)= 1000*(avg_pm3p(:,6)./ avg_pm3p(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm3p(:,9)= log(avg_pm3p(:,6))- log(avg_pm3p(:,7)); %col 11 is ln(R)

%% pm3a=pm3 anterior
px_pm3a= end_pm3a - start_pm3a +1;     %range of pixels that belongs to pm3a
    
           
    px=Data(:,5);
    ind_pm3a= find(px>=start_pm3a & px<=end_pm3a); % find the linear indices of the pixels that belongs to pm3a
    pm3a=Data(ind_pm3a,:);        %select those rows and all columns  from data matrices
    avg_pm3a=zeros(nframes, size(Data,2));   %create empty matrix to store aveg per cell
    for j=1:nframes           %for each block of pixels that defines a cell boundary
            i=j*px_pm3a;                       %and for each row of the avg matrix
            avg_pm3a(j,1:5)= pm3a((i-(px_pm3a-1)),1:5); %fill in metadata (col 1 to 7)
            avg_pm3a(j,6)= mean(pm3a(((i-(px_pm3a-1)):i),6));      %find average 410 per block of pixels
            avg_pm3a(j,7)= mean(pm3a(((i-(px_pm3a-1)):i),7));      %find average 470 per block of pixels
    end    
        avg_pm3a(:,8)= 1000*(avg_pm3a(:,6)./ avg_pm3a(:,7)); %col 10 is R1000 colm = 1000*410/470
        avg_pm3a(:,9)= log(avg_pm3a(:,6))- log(avg_pm3a(:,7)); %col 11 is ln(R)
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

    %Avg_perCell= [avg_pm7; avg_pm5; avg_pm3; avg_medax];
   Avg_perCell= [avg_pm8; avg_pm7; avg_pm6; avg_pm5; avg_pm4; avg_pm3; avg_pm3p; avg_pm3a; avg_medax];