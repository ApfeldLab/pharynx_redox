function FinalImage = loadTiffImageStack(filepath)
    % Read a TIFF image stack into a matrix.
    % 
    % Dimensions: width x height x worm
    %
    %
    % TODO: do the pictures come back as ints or doubles?
    
    InfoImage=imfinfo(filepath);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    FinalImage=zeros(nImage,mImage,NumberImages,'double');

    TifLink = Tiff(filepath, 'r');
    for i=1:NumberImages
       TifLink.setDirectory(i);
       FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
end