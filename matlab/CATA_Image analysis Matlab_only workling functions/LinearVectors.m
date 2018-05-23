%% PARAMETERS from user

% string metadata
expID='111031_HD233_Control_5mM';


% numeric metadata
date = 111031;
condition= 5;   % or '5' (if 1mM needed -> change second case in time vector generation

tp=61;
nworms=54;
n_plt_rep=3;
nframes=tp*nworms;

t_res=1;                %in mins or sec

nQ=3;
nworms_perQ=6;

delay_Q_1= 10;
delay_Q_2= 8;
delay_Q_3= 6;
delay_Q_4= 4;

Gen_Q_1= 233;     % Genotype or treatment per quadrant 
Gen_Q_2= 233; 
Gen_Q_3= 233; 
Gen_Q_4= 233; 





%Exception 1
if nworms_perQ * nQ * n_plt_rep~=nworms
 'total nworm is NOT equal to #wms per quadrant per plate'
end
 
%Exception 2
if condition== 0
        if kym_PreLev_410(1,1)==kym_PreLev_470(1,1)
        'IJ ERROR= Kymographs in the two wavelengths are identical'
        end
    elseif condition== 5

        if kym_5mM_410(1,1)==kym_5mM_470(1,1)
        'IJ ERROR= Kymographs in the two wavelengths are identical'
        end
end



%% 'string' metadata
 % col_expID= repmat(expID, 100*nframes, 1); 
col_date= repmat(date, 100*nframes, 1);
col_conc= repmat(condition, 100*nframes, 1);



plt_rep_mat=zeros(100,nframes);
 for i=1:n_plt_rep 
     plt_rep_block = repmat(i, 100, tp*(nworms/n_plt_rep)); %create block matrix ONLY works when each replicate has the same number of worms
     plt_rep_mat(: , ((i*tp*(nworms/n_plt_rep))-(tp*(nworms/n_plt_rep)-1)): i*tp*(nworms/n_plt_rep))= plt_rep_block;
 end
col_plt_rep= reshape(plt_rep_mat,(100*nframes),1);

strain_mat=zeros(100,nframes);
 for j=1:nQ 
     strain_perQ = repmat(eval(['Gen_Q_' int2str(j)]),100, tp*nworms_perQ);     %create block matrix
     strain_mat(: , ((j*tp*nworms_perQ)-(tp*nworms_perQ-1)): j*tp*nworms_perQ)= strain_perQ; %populate empty matrix with block matrix
 end
 col_strain=reshape(strain_mat,(100*nframes),1);


    
%% Generate linearized worm vector 
% worm number is NOT the worm ID. Each experiment counts worms independently and is determined by their position 
% in the imaging plate. The the right worm ID should include
% ExpName_date_genotype_plt_wormNumber
% Also generates the original frame number = columns

col_frameNum=[1:nframes]';

worm_mat=zeros(100,nframes);                %empty matrix to be filled up with 'worm block'
worm_n=ones(100,tp);                        %worm block matrix to be multiplied through loop
       
    for i=1:nworms       
            worm_mat(: , ((i*tp)-(tp-1)) : (i*tp))= worm_n*i; 
    end
col_worm=reshape(worm_mat,(100*nframes),1); %reshape matrix into a linearized vector

    
%% Generate time vector for either baseline OR response 

time_mat=zeros(100,nframes);
time_n=ones(100,1);

if condition== 0
    %% PreLev doesn't deals with quadrants delays and counts time in negative integers
    t_start=(0-((tp-1)*t_res));
    t_final=0;
    time_row= t_start: t_res : t_final ;
    
        for i=1:nworms
            time_mat(: , ((i*tp)-(tp-1)) : (i*tp))= time_n*time_row; 
        end
    col_time=reshape(time_mat,(100*nframes),1);
    %%
elseif condition== 5 %OJO!!! change to the TB concentration needed
    %% Response to oxidation considers time delays between quadrants
    for j=1:nQ 
        t_start= (0+ eval(['delay_Q_' int2str(j)]));
        t_final= ((tp*t_res - 1) + eval(['delay_Q_' int2str(j)]));
        time_row= t_start: t_res : t_final;
        
        time_block_perQ=repmat(time_row,100,nworms_perQ);  % concatenate horz "time_row" for nworms_perQ times
        
        time_mat(: , ((j*tp*nworms_perQ)-(tp*nworms_perQ-1)): j*tp*nworms_perQ)= time_block_perQ; %populate empty matrix with block matrix
    end
    col_time=reshape(time_mat,(100*nframes),1);
end




%% Generate linearized pixel coordinates vector

pixels=[1:100]';
px_mat=zeros(100,nframes);

    for i=1:nframes
       px_mat(:,i)=pixels;
    end
col_px=reshape(px_mat,(100*nframes),1);

%% Intensity Values in each channel and build the final Data_condition matrices

if condition== 0
    col_Lev_410=reshape(kym_PreLev_410,(100*nframes),1);
    col_Lev_470=reshape(kym_PreLev_470,(100*nframes),1);
    col_R1000= 1000*(col_Lev_410./col_Lev_470);
    col_ln= log(col_Lev_410) - log(col_Lev_470);

    Data_Lev=[col_date, col_plt_rep, col_conc, col_strain, col_worm, col_time, col_px, col_Lev_410, col_Lev_470, col_R1000, col_ln];

elseif condition== 5 %OJO!!! change to the TB concentration needed
    col_TB_410=reshape(kym_5mM_410,(100*nframes),1);
    col_TB_470=reshape(kym_5mM_470,(100*nframes),1);
    col_R1000 = 1000*(col_TB_410 ./ col_TB_470);
    col_ln= log(col_TB_410) - log(col_TB_470);

    Data_TB=[col_date, col_plt_rep, col_conc, col_strain, col_worm, col_time, col_px, col_TB_410, col_TB_470, col_R1000, col_ln];
    
end

% %% Find cellular regions and take the mean
% %col_cellRegion %use ifelse 4 ranges to assign label, how to do this from one column to another?
% start_pm7=3;        %start of cell boundary
% end_pm7=10;         %end of cell boundary
% px_pm7= end_pm7 - start_pm7 +1;     %range of pixels that belongs to pm7
%     % find the linear indices of the pixels that belongs to pm7
% 
% if condition== 0            
%     px=Data_Lev(:,7);
%     ind_pm7= find(px>=start_pm7 & px<=end_pm7); 
%     pm7_Lev=Data_Lev(ind_pm7,:);        %select those rows and all columns  from data matrices
%     avg_pm7_Lev=zeros(nframes, size(Data_Lev,2));   %create empty matrix to store aveg per cell
%     for j=1:nframes           %for each block of pixels that defines a cell boundary
%             i=j*px_pm7;                       %and for each row of the avg matrix
%             avg_pm7_Lev(j,1:7)= pm7_Lev((i-(px_pm7-1)),1:7); %fill in metadata (col 1 to 7)
%             avg_pm7_Lev(j,8)= mean(pm7_Lev(((i-(px_pm7-1)):i),8));      %find average 410 per block of pixels
%             avg_pm7_Lev(j,9)= mean(pm7_Lev(((i-(px_pm7-1)):i),9));      %find average 470 per block of pixels
%     end    
%         avg_pm7_Lev(:,10)= 1000*(avg_pm7_Lev(:,8)./ avg_pm7_Lev(:,9)); %col 10 is R1000 colm = 1000*410/470
%         avg_pm7_Lev(:,11)= log(avg_pm7_Lev(:,8))- log(avg_pm7_Lev(:,9)); %col 11 is ln(R)
% 
% elseif condition== 5
%     px=Data_TB(:,7);
%     ind_pm7= find(px>=start_pm7 & px<=end_pm7); 
%     pm7_TB=Data_TB(ind_pm7,:);        
%     avg_pm7_TB=zeros(nframes, size(Data_TB,2));
%     for j=1:nframes           %for each block of pixels that defines a cell boundary
%             i=j*px_pm7;                       %and for each row of the avg matrix
%             avg_pm7_TB(j,1:7)= pm7_TB((i-(px_pm7-1)),1:7); %fill in metadata (col 1 to 7)
%             avg_pm7_TB(j,8)= mean(pm7_TB(((i-(px_pm7-1)):i),8));      %find average 410 per block of pixels
%             avg_pm7_TB(j,9)= mean(pm7_TB(((i-(px_pm7-1)):i),9));      %find average 470 per block of pixels
%     end
%         avg_pm7_TB(:,10)= 1000*(avg_pm7_TB(:,8)./ avg_pm7_TB(:,9)); %col 10 is R1000 colm = 1000*410/470
%         avg_pm7_TB(:,11)= log(avg_pm7_TB(:,8))- log(avg_pm7_TB(:,9)); %col 11 is ln(R)    
% end


%% Use Cell arrays to store string and numeric data?
% mycell = { 'a' 1 2 3 ; 'b' 4 5 6 };
% [nrows,ncols]= size(mycell);
% 
% filename = 'celldata.dat';
% fid = fopen(filename, 'w');
% 
% for row=1:nrows
%     fprintf(fid, '%s, %10.10f, %10.10f, %10.10f,\n', mycell{row,:}); % '10.10f' is a float precision with 10 digits to the left and 10 to the right
% end
% 
% fclose(fid);
% type celldata.dat %to view
