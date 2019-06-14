%% PARAMETERS from user

% string metadata
expID='140710_11_HD233_3doF_Rep1';

% numeric metadata
Date = 140710.11;
Gen= 233;
age= 3;

% var
%R=(col_Lev_410./col_Lev_470);

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
    col_R=(col_Lev_410./col_Lev_470);
    col_OxD= (col_R - 0.667)./((col_R - 0.667)+(5.207 - col_R)*0.171);
   % col_E= (-265)+((8.315 * (273.15 + 22))./ (2 * 96.48)) * log(col_OxD./
   % (1 - col_OxD)); the extra digits don't explain the difference between
   % JMP result and ML result
    col_E= (-265)+((8.314462 * (273.15 + 22))./ (2 * 96.4853415)) * log(col_OxD./ (1 - col_OxD)); % from ja code
        % E0 = -265;
        % Temp = 22;
        % E = E0 - (8314.462*(273.15+Temp)/(2*96485.3415))*log((1-OxD)./OxD);
        % E = real(E);

    Data=[col_Date, col_Gen, col_age, col_Worm, col_px, col_Lev_410, col_Lev_470, col_R, col_OxD,col_E];
    Table=array2table(Data, 'VariableNames',{'Date' 'Strain' 'Age' 'worm' 'px' 'Int410' 'Int470' 'R' 'OxD' 'E'});
