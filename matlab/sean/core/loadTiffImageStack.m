function FinalImage = loadTiffImageStack(filepath)
    % Read a TIFF image stack into a matrix (width x height x nFrames) 
    
    InfoImage=imfinfo(filepath);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    FinalImage=zeros(nImage,mImage,NumberImages,'uint16');

    TifLink = Tiff(filepath, 'r');
    for i=1:NumberImages
       TifLink.setDirectory(i);
       FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
end