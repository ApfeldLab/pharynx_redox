%% PARAMETERS from user

% string metadata
expID='120717_18_HD240_daf2';

% numeric metadata
Date = 120718.19;
Gen= 240;
age= 8;

%% Linerrized Vectors

%Exception 1
    if kym_410(1,1)==kym_470(1,1)
        'IJ ERROR= Kymographs in the two wavelengths are identical'
    end

nworms= size(kym_410,2);
col_Gen= repmat(Gen, 100*nworms, 1);
col_Date= repmat(Date, 100*nworms, 1);
col_age= repmat(age, 100*nworms, 1);

worm_mat=zeros(100, nworms);     % Generate linearized worm number vector
    for i=1:nworms
            worm_mat(: , i)=i;
    end
col_Worm=reshape(worm_mat,(100*nworms),1); 


pixels=[1:100]';    % Generate linearized pixel coordinates vector
px_mat=zeros(100,nworms);

    for i=1:nworms
       px_mat(:,i)=pixels;
    end
col_px=reshape(px_mat,(100*nworms),1);


% Intensity Values in each channel and build the final Data_Gen matrices

    col_Lev_410=reshape(kym_410,(100*nworms),1);
    col_Lev_470=reshape(kym_470,(100*nworms),1);
    col_R1000= 1000*(col_Lev_410./col_Lev_470);
    col_ln= log(col_Lev_410) - log(col_Lev_470);

    Data=[col_Date, col_Gen, col_age, col_Worm, col_px, col_Lev_410, col_Lev_470, col_R1000, col_ln];

%     LATER ON name kym_'Gen'_410 and Data_'Gen'

%     Data_(eval(['Gen']))= 240;
%     Data_'int2str(eval(['Gen']))'=2;
%       
%     Data_(disp(Gen)