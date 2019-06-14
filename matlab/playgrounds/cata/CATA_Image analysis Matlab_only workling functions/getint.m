function kym = getint(lnorm)
%   GETINT extracts the intensity values (lny vectors)from the lenght
%   normalized matrix (lnorm) wich results from the SQUARE fuction.
%   This lnorm matrix has the pixel coordinates vectors (lnx) intercalated
%   with the intensity vectors


sq=lnorm; % I am duplicating this var to avoid overwriting by now
nrows=size(sq,1);
ncols=size(sq,2); 
New=zeros(nrows,ncols/2);
for k=2:2:ncols  %Counter for the yi columns of initial matrix
    j=k/2;    %Counter for the j-th column of new (intensity) matrix
    New(:,j)=sq(:,k);
    
end

kym=New;
%plot(kym)