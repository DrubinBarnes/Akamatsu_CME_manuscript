
function im = openTiffStack(tsfilename)

% based on suggestion at http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/

% you can access the image by FinalImage(rows,cols, frame)
% example: tsfilename = '/Users/Nathaniel/Desktop/TestTiffStack.tif';

InfoImage=imfinfo(tsfilename);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
im=zeros(nImage,mImage,NumberImages,'uint16');

TifLink = Tiff(tsfilename, 'r');
for i=1:NumberImages
	TifLink.setDirectory(i);
	im(:,:,i)=TifLink.read();
end
TifLink.close();