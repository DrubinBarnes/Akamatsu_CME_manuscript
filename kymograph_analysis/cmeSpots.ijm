// List of buttons for analyzing CCPs in TIRF

// Macro to make and save kymographs from line regions
// First select one or more line regions using the Line Tool and save each to the ROI manager using "t"


var kymographWidth = 6 //in pixels. The width of the region from which to make kymograph.


var measureInitiations = 0; // Set this to 1 if you want to measure initiation densities.
var timeinterval = 2; // in seconds. 
var micronsperpixel = 0.108; // default for IXACT camera 60x objective

//var micronsperpixel = 0.108*2; // IXACT camera 60x objective with 2x2 binning
var LUTnb=1; // 1: Grays; 2: inverted grays; 3: fire
var useSubset = 1; // Set this to 1 if you want to measure intitation densities for a (continuous) region smaller than the size of the kymograph.

macro "Make kymographs (shift: batch mode) Action Tool - C059T3e16K" {
	
	// This tool will make kymographs from line regions you create.
	
	if(isKeyDown("shift"))
	{fastmethod=1;
		
		print("This is the batch method of making kymographs. All open images have been closed. Your kymographs will be saved in the same folder as your movie.");
	}
	
	else
	{fastmethod = 0;
	}
	
	if(fastmethod>0){
		
		// Close all windows without saving
		
		while (nImages>0) { 
			selectImage(floor(nImages*random)+1); 
			close(); 
		}
		
		
		imageFile=File.openDialog("Select your movie");
		setBatchMode(true);	
		open(imageFile);
		
		dir=getInfo("image.directory");
		curID=getImageID();
		imageTitle=getTitle();
		shortImageTitle=replace(imageTitle,".tif","_");
		
		roiNb=roiManager("count");
		
		// save ROIs
		print(dir+shortImageTitle+"ROIs.zip");
		
		if(roiNb == 0 ){
			
			kymROIs = File.openDialog("Select ROI zip file to make kymographs");
			roiManager("Open", kymROIs); 
			roiNb=roiManager("count");
			
		}
		else {
			roiManager("deselect");
			roiManager("Remove Slice Info");
			roiManager("Save", dir+shortImageTitle+"w"+kymographWidth+"_r1-"+roiNb+"ROIs.zip");
		}
		
		// Iterate through ROIs
		for(i=0;i<roiNb;i++){		
			
			ii = i+1; // so the naming starts at 1 instead of 0
			
			selectImage(curID);
			roiManager("select",i);
			// straighten image - width defined above
			run("Straighten...", "line="+kymographWidth+" process");
			id1 = getImageID();
			// reslice
			
			run("Reslice [/]...", "output=1.000 start=Top avoid");
			id2 = getImageID();
			
			// maximum project image
			
			run("Z Project...", "projection=[Max Intensity]");
			id3 = getImageID();
			curKym = getTitle();
			
			// rotate 90 degrees left
			
			run("Rotate 90 Degrees Left");
			
			// Change color
			
			//run("Green Fire Blue"); // If you don't have this particular lookup table or want b & w, comment this out.
			
			// save
			
			//saveAs("TIF", dir+shortImageTitle+"r"+ii+"w"+kymographWidth+"_kymograph.tif");
			//saveAs("png", dir+shortImageTitle+"r"+ii+"w"+kymographWidth+"_kymograph.png");
			
			// close intermediate images
			
			selectImage(id1);
			run("Close");
			
			selectImage(id2);
			run("Close");
			
			// FOR BATCH MODE
			// Try to make a stack from both images, turn into composite, and change look up tables. to circumvent protblem from merging
			
			
			//wait(500);
			//mergedID = getImageID();
			
			selectImage(id3); 
			
			saveAs("TIF", dir+shortImageTitle+"_w"+kymographWidth+"r"+ii+"_kymograph.tif");
			saveAs("PNG", dir+shortImageTitle+"w"+kymographWidth+"r"+ii+"_kymograph.png");
			
			selectImage(shortImageTitle+"w"+kymographWidth+"r"+ii+"_kymograph.png");
			run("Close");
			
			
			
			
			
			
		}
		
		setBatchMode("exit and diplay");
	}
	
	else
	{
		
		dir=getInfo("image.directory");
		curID=getImageID();
		imageTitle=getTitle();
		shortImageTitle=replace(imageTitle,".tif","_");
		
		roiNb=roiManager("count");
		// save ROIs
		print(dir+shortImageTitle+"ROIs.zip");
		
		if(roiNb == 0 ){
			
			kymROIs = File.openDialog("Select ROI zip file to make kymographs");
			roiManager("Open", kymROIs); 
			roiNb=roiManager("count");
			
		}
		else {
			
			roiManager("deselect");
			roiManager("Remove Slice Info");
			roiManager("Save", dir+shortImageTitle+"w"+kymographWidth+"_r1-"+roiNb+"ROIs.zip");
		}
		
		// Iterate through ROIs
		
		for(i=0;i<roiNb;i++){
			
			ii = i+1; // so the naming starts at 1 instead of 0
			
			selectImage(curID);
			roiManager("select",i);
			// straighten image - width defined above
			run("Straighten...", "line="+kymographWidth+" process");
			id1 = getImageID();
			// reslice
			
			run("Reslice [/]...", "output=1.000 start=Top avoid");
			id2 = getImageID();
			
			// maximum project image. If you prefer sum projection: comment this out and uncomment "Sum slices" line
			
			run("Z Project...", "projection=[Max Intensity]");
			//		run("Z Project...", "projection=[Sum Slices]");
			
			id3 = getImageID();;
			
			// rotate 90 degrees left
			
			run("Rotate 90 Degrees Left");
			
			// Change color
			
			// run("Green Fire Blue"); 
			// If you prefer green/ blue kymographs, remove the "//" at the beginning of the line above.
			
			// save
			
			saveAs("TIF", dir+shortImageTitle+"_w"+kymographWidth+"r"+ii+"_kymograph.tif");
			saveAs("png", dir+shortImageTitle+"w"+kymographWidth+"r"+ii+"_kymograph.png");
			
			// close intermediate images
			
			selectImage(id1);
			run("Close");
			
			selectImage(id2);
			run("Close");
			
		}
	}
	// Save image of regions + original movie file
	
	selectImage(curID);
	roiManager("deselect");
	run("Duplicate...", "title="+shortImageTitle+"w"+kymographWidth+"r1-"+roiNb+"_regionOverlay.tif");
	run("Enhance Contrast", "saturated=0.5");
	roiManager("Set Line Width", kymographWidth);
	roiManager("Remove Slice Info");
	roiManager("Set Color", "#4dffff00");
	roiManager("Show All with labels");
	run("From ROI Manager");
	saveAs("png", dir+shortImageTitle+"w"+kymographWidth+"r1-"+roiNb+"_regionOverlay.png");
	
	//run("Close");
	
	print("kymograph width is "+kymographWidth+" pixels");
	
}
macro 'OptimizeContrast [o]'{

	// With this toolset installed, you can optimize the brightness/contrast by pressing the "o" key on your keyboard.

	if (LUTnb==1){
		LUTnb=2;
		run("Grays");
		run("Invert LUT"); //run("Rainbow RGB");
		run("Enhance Contrast", "saturated=0.5");
	}
	else if (LUTnb==2){
		LUTnb=3;
		run("Fire");
		//run("Enhance Contrast", "saturated=0.5");
	}
	else if (LUTnb==3){
		LUTnb=1;
		run("Grays");
		//run("Invert LUT"); //run("Rainbow RGB");
		run("Enhance Contrast", "saturated=0.5");
	}
}

macro "Make color composite Action Tool - Cff0D63D64D65D71D72D73D74D75D81D82D83D84D85D93D94D95C0f0D90D91D92Da0Da1Da2Da3Da4Da5Da6Db0Db1Db2Db3Db4Db5Db6Dc0Dc1Dc2Dc3Dc4Dc5Dc6Dc7Dc8Dd1Dd2Dd3Dd4Dd5Dd6Dd7Dd8De1De2De3De4De5De6De7De8Df3Df4Df5Df6Cf0fD39D47D48D49D57D58D59D67D68D69Cf00D03D04D05D06D11D12D13D14D15D16D17D18D21D22D23D24D25D26D27D28D30D31D32D33D34D35D36D37D38D40D41D42D43D44D45D46D50D51D52D53D54D55D56D60D61D62C00fD3aD3bD3cD4aD4bD4cD4dD4eD5aD5bD5cD5dD5eD6aD6bD6cD6dD6eD6fD79D7aD7bD7cD7dD7eD7fD89D8aD8bD8cD8dD8eD8fD9aD9bD9cD9dD9eD9fDaaDabDacDadDaeDbaDbbDbcDbdDbeDcaDcbDccCfffD66D76D77D78D86D87D88D96C0ffD97D98D99Da7Da8Da9Db7Db8Db9Dc9"
// This button makes two-color composite images or movies from two single-channel movies.

{
	
	if(isKeyDown("shift"))
	{greenRed=1;
		print("You're making a green/red composite");
	}
	
	else
	{greenRed = 0;
	}
	
	greenFile=File.openDialog("Select GFP file");
	redFile  =File.openDialog("Select RFP file");
	
	print(greenFile);
	
	open(greenFile);
	gfpImage = getTitle();
	
	open(redFile);
	rfpImage = getTitle();
	
	print(gfpImage);
	print(rfpImage);
	
	dir=getInfo("image.directory");
	curID=getImageID();
	imageTitle=getTitle();
	shortImageTitle=replace(imageTitle,".tif","_");
	
	if(greenRed>0)
	{run("Merge Channels...", "c1=["+rfpImage+"] c2=["+gfpImage+"] create ignore");
	}
	
	else
	{run("Merge Channels...", "c6=["+rfpImage+"] c2=["+gfpImage+"] create ignore");
	}
	
	// This sets the contrast of the two colors such that the red channel is 3x brighter than the green channel.
	
	setSlice(1);		
	
	run("Enhance Contrast", "saturated=0.25");
	
	
	setSlice(2);		
	run("Enhance Contrast", "saturated=0.75");
	
	
	
	saveAs("TIF", dir+shortImageTitle+"composite.tif");
	
}

macro "Make color kymographs Action Tool - C059T3e16C"
// This tool will make two-color kymographs from line regions you create.
// You need to have two movies, one for each channel. 
// You can make these movies with the “Image Registration” application made by Sun Hong.

{
	
	if(isKeyDown("shift"))
	{greenRed=1;
		print("You're making a green/red composite");
	}
	
	else
	{greenRed = 0;
	}
	
	
	// Close all windows without saving
	
	while (nImages>0) { 
		selectImage(floor(nImages*random)+1); 
		close(); 
	}
	
	
	/*
	list = getList("window.titles");
	if (list.length>0) for (i=0; i<list.length; i++) {
	selectWindow(list[i]);
	run ("Close");
	}
	
	ids=getIDs();
	if (ids.length>0) for (i=0;i<ids.length;i++){
	selectImage(ids[i]);
	run ("Close");
	}
	*/
	
	greenFile=File.openDialog("Select GFP file");
	redFile  =File.openDialog("Select RFP file");
	
	setBatchMode(true);	
	
	
	print(greenFile);
	
	open(greenFile);
	gfpImage = getTitle();
	gfpID = getImageID();
	
	open(redFile);
	rfpImage = getTitle();
	rfpID = getImageID();
	
	dir=getInfo("image.directory");
	curID=getImageID();
	imageTitle=getTitle();
	shortImageTitle=replace(imageTitle,".tif","_");
	
	roiNb=roiManager("count");
	// save ROIs
	print(dir+shortImageTitle+"ROIs.zip");
	
	if(roiNb == 0 ){
		
		kymROIs = File.openDialog("Select ROI zip file to make kymographs");
		roiManager("Open", kymROIs); 
		roiNb=roiManager("count");
		
	}
	else {
		roiManager("deselect");
		roiManager("Remove Slice Info");
		//	roiManager("Save", dir+shortImageTitle+"ROIs.zip");
		roiManager("Save", dir+shortImageTitle+"w"+kymographWidth+"_r1-"+roiNb+"ROIs.zip");
		// Iterate through ROIs
	}
	
	
	
	for(i=0;i<roiNb;i++){
		//	for(i=0;i<2;i++){
		
		
		ii = i+1; // so the naming starts at 1 instead of 0
		
		selectImage(gfpID);
		roiManager("select",i);
		// straighten image - width defined above
		run("Straighten...", "line="+kymographWidth+" process");
		id1g = getImageID();
		// reslice
		
		run("Reslice [/]...", "output=1.000 start=Top avoid");
		id2g = getImageID();
		
		// maximum project image
		
		run("Z Project...", "projection=[Max Intensity]");
		id3g = getImageID();
		greenKym = getTitle();
		
		// rotate 90 degrees left
		
		run("Rotate 90 Degrees Left");
		
		// Change color
		
		//run("Green Fire Blue"); // If you don't have this particular lookup table or want b & w, comment this out.
		
		// save
		
		//saveAs("TIF", dir+shortImageTitle+"r"+ii+"w"+kymographWidth+"_kymograph.tif");
		//saveAs("png", dir+shortImageTitle+"r"+ii+"w"+kymographWidth+"_kymograph.png");
		
		// close intermediate images
		
		selectImage(id1g);
		run("Close");
		
		selectImage(id2g);
		run("Close");
		
		// For RFP image
		
		selectImage(rfpID);
		roiManager("select",i);
		// straighten image - width defined above
		run("Straighten...", "line="+kymographWidth+" process");
		id1r = getImageID();
		// reslice
		
		run("Reslice [/]...", "output=1.000 start=Top avoid");
		id2r = getImageID();
		
		// maximum project image
		
		run("Z Project...", "projection=[Max Intensity]");
		id3r = getImageID();
		redKym = getTitle();
		
		// rotate 90 degrees left
		
		run("Rotate 90 Degrees Left");
		
		// close intermediate images
		
		selectImage(id1r);
		run("Close");
		
		selectImage(id2r);
		run("Close");
		
		// Merge green and red kymographs
		
		// As magenta/green if shift not pressed; as red/green if shift pressed
		
		//setBatchMode(false);
		
		//wait(500);
		
		// FOR BATCH MODE
		// Try to make a stack from both images, turn into composite, and change look up tables. to circumvent protblem from merging
		
		
		if(greenRed>0)
			
			//{run("Merge Channels...", "c1=["+redKym+"] c2=["+greenKym+"] create ignore");
		
		{run("Images to Stack", "name=StackTwoColor use");
			run("Make Composite", "display=Composite");
			
			setSlice(1);		
			run("Green");			
			
			setSlice(2);		
			run("Red");		
		}
		
		else
			//{run("Merge Channels...", "c6=["+redKym+"] c2=["+greenKym+"] create ignore");
		
		{run("Images to Stack", "name=StackTwoColor use");
			run("Make Composite", "display=Composite");
			
			setSlice(1);		
			run("Green");			
			
			// set contrast of GFP channel severely
	//		setMinAndMax(255, 386); 
			
			
			setSlice(2);		
			run("Magenta");			
			
		}
		
		
		
		
		//wait(500);
		mergedID = getImageID();
		
		//setBatchMode(true);
		
		//saveAs("TIF", dir+shortImageTitle+"composite.tif");
		//selectImage("Composite"); // THIS LINE IS NOT SO STABLE: it requires the merged image to appear with the title "Composite"
		//		selectImage("StackTwoColor"); 
		selectImage(mergedID); 
		
		if(greenRed>0)
			
		{saveAs("TIF", dir+shortImageTitle+"_w"+kymographWidth+"rg_r"+ii+"_composite_kymograph.tif");
			saveAs("PNG", dir+shortImageTitle+"w"+kymographWidth+"rg_r"+ii+"_composite_kymograph.png");
			selectImage(shortImageTitle+"w"+kymographWidth+"rg_r"+ii+"_composite_kymograph.png");
			run("Close");
		}
		
		else
			
		{saveAs("TIF", dir+shortImageTitle+"_w"+kymographWidth+"r"+ii+"_composite_kymograph.tif");
			saveAs("PNG", dir+shortImageTitle+"w"+kymographWidth+"r"+ii+"_composite_kymograph.png");
			
			selectImage(shortImageTitle+"w"+kymographWidth+"r"+ii+"_composite_kymograph.png");
			run("Close");
		}
		
		
		
		
		
	}
	
	setBatchMode(false);		
	
	
	// Save image of regions + original movie file
	
	selectImage(rfpID);	
	roiManager("deselect");
	run("Duplicate...", "title="+shortImageTitle+"w"+kymographWidth+"r1-"+roiNb+"_regionOverlay.tif");
	
	run("Enhance Contrast", "saturated=0.5");
	roiManager("Set Line Width", kymographWidth);
	roiManager("Remove Slice Info");
	roiManager("Set Color", "#4dffff00");
	roiManager("Show All with labels");
	//roiManager("Show None");
	//roiManager("Show All with labels");		
	run("From ROI Manager");
	saveAs("png", dir+shortImageTitle+"w"+kymographWidth+"_r1-"+roiNb+"_regionOverlay.png");
		
}	

macro "Maximum (+shift for Sum) Intensity Projection Action Tool - C902T3f18Z"
{
	if (Stack.isHyperStack)
	{
		getDimensions(w, h, channels, slices, frames);
		if(isKeyDown("shift"))
			{run("Z Project...", "start=1 stop="+slices+" projection=[Sum Slices] all");}
		else 
			{run("Z Project...", "start=1 stop="+slices+" projection=[Max Intensity] all");}
	}
	else
		
		if(isKeyDown("shift"))
		{run("Z Project...", "start=1 stop="+nSlices+" projection=[Sum Slices]");}
	else 
	{
		run("Z Project...", "start=1 stop="+nSlices+" projection=[Max Intensity]");
	}
	
	
}

macro "Cut out box from composite Action Tool - C059T3e16B"
{
	dir=getInfo("image.directory");
	curID=getImageID();
	imageTitle=getTitle();
	shortImageTitle=replace(imageTitle,".tif","_");
	count=roiManager("count");

	roiManager("deselect");
	roiManager("Remove Slice Info");
	roiManager("Save", dir+shortImageTitle+"boxROIs.zip");
	
	for(i=0;i<count;i++){
		
		j = i+1;
		
		selectImage(curID);
		roiManager("select",i);
	
		run("Duplicate...", "duplicate title=["+shortImageTitle+"box"+j+".tif]");
		saveAs("TIF", dir+shortImageTitle+"box"+j+".tif");

	}
	
	// Save image of regions + original movie file
	
	selectImage(curID);	
	roiManager("deselect");
	run("Select All");
	run("Duplicate...", "title="+shortImageTitle+"box1-"+count+"_regionOverlay.tif");
	curDupID=getImageID();
	
	//run("Enhance Contrast", "saturated=0.5");
	run("Stack to RGB");
	roiManager("Set Line Width", 5);
	roiManager("Remove Slice Info");
	roiManager("Set Color", "#4dffff00");
	roiManager("deselect");
	roiManager("Show None");
	roiManager("Show All with labels");
	run("From ROI Manager");
	saveAs("png", dir+shortImageTitle+"box1-"+count+"_regionOverlay..png");
	
	selectImage(curDupID);	
	run("Close");
}


//macro "Temporal-Color Coder Action Tool C059T3e16T"
//{
/*

************* Temporal-Color Coder *******************************
Color code the temporal changes.

Kota Miura (miura@embl.de) +49 6221 387 404 
Centre for Molecular and Cellular Imaging, EMBL Heidelberg, Germany

!!! Please do not distribute. If asked, please tell the person to contact me. !!!
If you publish a paper using this macro, it would be cratedful if you could acknowledge. 

Edit MA 08/10/2015 to save image with file name
 
---- INSTRUCTION ----

1. Open a stack (8 bit or 16 bit)
2. Run the macro
3. In the dialog choose one of the LUT for time coding.
	select frame range (default is full).
	check if you want to have color scale bar.

History

080212	created ver1 K_TimeRGBcolorcode.ijm
080213	stack slice range option added.
		time color code scale option added.

		future probable addiition: none-linear assigning of gray intensity to color intensity
		--> but this is same as doing contrast enhancement before processing.
*****************************************************************************
*/


var Glut = "Fire";	//default LUT
var Gstartf = 1;
var Gendf = 10;
var GFrameColorScaleCheck=0;

macro "Time Lapse Color Coder Action Tool - C059T3e16T"{
		
	dir=getInfo("image.directory");
	curID=getImageID();
	imageTitle=getTitle();
	shortImageTitle=replace(imageTitle,".tif","_");
	
	
	Gendf = nSlices;
	Glut = ChooseLut();
	run("Duplicate...", "title=listeriacells-1.stk duplicate");
	run("Enhance Contrast", "saturated=0.5");

	hh = getHeight();
	ww = getWidth();
	totalslice = nSlices;
	calcslices = Gendf - Gstartf +1;
	run("8-bit");
	imgID = getImageID();

	newImage("colored", "RGB White", ww, hh, calcslices);
	newimgID = getImageID();

	setBatchMode(true);

	newImage("stamp", "8-bit White", 10, 10, 1);
	run(Glut);
	getLut(rA, gA, bA);
	close();
	nrA = newArray(256);
	ngA = newArray(256);
	nbA = newArray(256);

	for (i=0; i<calcslices; i++) {
		colorscale=floor((256/calcslices)*i);
		//print(colorscale);
		for (j=0; j<256; j++) {
			intensityfactor=0;
			if (j!=0) intensityfactor = j/255;
			nrA[j] = round(rA[colorscale] * intensityfactor);
			ngA[j] = round(gA[colorscale] * intensityfactor);
			nbA[j] = round(bA[colorscale] * intensityfactor);
		}
		newImage("temp", "8-bit White", ww, hh, 1);
		tempID = getImageID;

		selectImage(imgID);
		setSlice(i+Gstartf);
		run("Select All");
		run("Copy");

		selectImage(tempID);
		run("Paste");
		setLut(nrA, ngA, nbA);
		run("RGB Color");
		run("Select All");
		run("Copy");
		close();

		selectImage(newimgID);
		setSlice(i+1);
		run("Select All");
		run("Paste");
	}

	selectImage(imgID);
	close();
	selectImage(newimgID);
	op = "start=1 stop="+totalslice+" projection=[Max Intensity]";
	run("Z Project...", op);
	setBatchMode(false);	
	if (GFrameColorScaleCheck) CreatGrayscale256(Glut, Gstartf, Gendf);
	temporalColorID=getImageID();

	selectImage(newimgID);
	close();
		
	selectImage(temporalColorID);

	saveAs("TIF", dir+shortImageTitle+"temporalColor__t"+Gstartf+"-"+Gendf+".tif");
	saveAs("png", dir+shortImageTitle+"temporalColor_t"+Gstartf+"-"+Gendf+".png");
		
	
}

/*
run("Spectrum");
run("jet"); //MA
run("Fire");
run("Ice");
run("3-3-2 RGB");
run("brgbcmyw");
run("Green Fire Blue");
run("royal");
run("thal");
run("smart");
run("unionjack");
run("jet"); //MA 10-28-14
run("5_ramps");//MA 6-4-13
run("phase");//MA
14 luts
*/

function ChooseLut() {
	lutA=newArray(10);
	lutA[0] = "jet"; //MA
	lutA[1] = "Fire";
	lutA[2] = "Ice";
	lutA[3] = "3-3-2 RGB";
	lutA[4] = "brgbcmyw";
	lutA[5] = "Green Fire Blue";
	lutA[6] = "royal";
	lutA[7] = "Spectrum"; //MA
	//lutA[7] = "thal";
	lutA[8] = "smart";
	//lutA[9] = "unionjack";
	lutA[9] = "Red Hot";

 	Dialog.create("Color Code Settings");
	Dialog.addChoice("LUT", lutA);
	Dialog.addNumber("start frame", Gstartf);
	Dialog.addNumber("end frame", Gendf);
	Dialog.addCheckbox("create Time Color Scale Bar", GFrameColorScaleCheck);
 	Dialog.show();
 	Glut = Dialog.getChoice();
	Gstartf= Dialog.getNumber();
	Gendf= Dialog.getNumber();
	GFrameColorScaleCheck = Dialog.getCheckbox();
	print("selected lut:"+ Glut);
	return Glut;
}

function CreatGrayscale256(lutstr, beginf, endf){
	ww = 256;
	hh=32;
	newImage("color time scale", "8-bit White", ww, hh, 1);
	for (j=0; j<hh; j++) {
		for (i=0; i<ww; i++) {
			setPixel(i, j, i);
		}
	}
	run(lutstr);
	//setLut(nrA, ngA, nbA);
	run("RGB Color");
	op = "width="+ww+" height="+hh+16+" position=Top-Center zero";
	run("Canvas Size...", op);
	setFont("SansSerif", 12, "antiliased");
	run("Colors...", "foreground=white background=black selection=yellow");
	drawString("frame", round(ww/2)-12, hh+16);
	drawString(leftPad(beginf, 3), 0, hh+16);
	drawString(leftPad(endf, 3), ww-24, hh+16);

}

function leftPad(n, width) {
    s =""+n;
    while (lengthOf(s)<width)
        s = "0"+s;
    return s;
}
/*
macro "drawscale"{
	CreatGrayscale256("Fire", 1, 100);
}
*/


macro "Adjust registration Action Tool - C059T3e16R"{
	
	if(isKeyDown("shift"))
		{greenRed=1;
		print("You're making a green/red composite");
		}
	
	else
		{greenRed = 0;
		}
	
	regFile=File.openDialog("Select file to change registration");

	// Change "xShift" and "yShift" values to change registration of selected image, in pixels
	
	xShift = 5;
	yShift = 40;
	
	open(regFile);
	regImage = getTitle();
		
	dir=getInfo("image.directory");
	curID=getImageID();
	imageTitle=getTitle();
	shortImageTitle=replace(imageTitle,".tif","_");
	
	run("Translate...", "x="+xShift+" y="+yShift+" interpolation=None stack");
	
	// This sets the contrast of the two colors such that the red channel is 3x brighter than the green channel.
	
	saveAs("TIF", dir+shortImageTitle+"translateX"+xShift+"Y"+yShift+".tif");
	
}

macro "Measure line from kymograph Action Tool - C059T3e16L"
// This tool measures lines you have drawn on kymographs and saves them with the same name as the kymograph.
// Option for measuring just between the top and bottom most lines. If you've measured just a subest of the tracks in the kymograph.
{
	
	dir=getInfo("image.directory");
	curID=getImageID();
	imageTitle=getTitle();
	shortImageTitle=replace(imageTitle,".tif","_");
	roiNb=roiManager("count");
	
	curSlice = getSliceNumber();
	channelImageTitle = shortImageTitle+"c"+curSlice+"_";
	
	roiManager("deselect");
	roiManager("Set Color", "#4dffff00");
	roiManager("Save", dir+channelImageTitle+"lineROIs.zip");
	print(dir+channelImageTitle+"lineROIs.zip saved");	
	
	run("Set Scale...", "distance=0");
	
	if(isOpen("Results")){;
		//			IJ.deleteRows(0,nResults);
		selectWindow("Results"); 
		run("Close"); 
	}
	
	
	
	roiManager("deselect");
	roiManager("Measure");
	
	saveAs("measurements", dir+channelImageTitle+"lineROIs.xls");
	saveAs("measurements", dir+channelImageTitle+"lineROIs.txt");
	
	print( dir+channelImageTitle+"lineROIs.xls and .txt saved");
	roiManager("deselect");
	roiManager("Delete");
	
	// Get all the Y coordinates of the lines measured in results
	
	ys = newArray(nResults);
	
	for(i=0;i<nResults;i++){ 
		ys[i] = getResult("Y", i);
	}
	
	Array.getStatistics(ys,ymin,ymax);
	
	
	//// Measure and save density of events. In units of events/um^2/minute.
	
	if(measureInitiations>0)
     	{
     		heightpx = getHeight();
     		widthpx  = getWidth();
     		nEvents = nResults;
     		
     		subsetheightpx = ymax-ymin; // This is the distance between the top and bottom line ROI
     		
     		unscaleddensity = nEvents/heightpx/widthpx/kymographWidth; // in events per pixels squared per time interval
     		unscaleddensitysubset = nEvents/subsetheightpx/widthpx/kymographWidth; // in events per pixels squared per time interval. For just region between top and botom line ROI, not entire kymograph.

     		// convert to microns, seconds
     		
     		heightmicron = heightpx * micronsperpixel; // e.g. .108 microns per pixel
     		subsetheightmicron = subsetheightpx * micronsperpixel;
     		widthminutes  = widthpx  * timeinterval /60;    // e.g. 2 seconds per time point (pixel) divided by 60 seconds/ min
     		kymographwidthmicron = kymographWidth * micronsperpixel;
     		
     		density = nEvents/heightmicron/widthminutes/kymographwidthmicron; // in events per pixels squared per time interval
     		
     		densitysubset = nEvents/subsetheightmicron/widthminutes/kymographwidthmicron; // in events per pixels squared per time interval
     		
     		
     		
     		//run("Text Window...", "name="+tableTitle+" width=72 height=8 menu");
     		
     		//print(f, "\\Clear");
     		//print(f, density);
     		
     		
     		if(isOpen("Results")){;	
     			selectWindow("Results"); 
     			run("Close"); 
     		}
		
     		if(useSubset>0){ // Report the initiation density for a subset between Y max and min of line ROIs, if "var useSubset" is 1 at top of this file
     		
     			setResult("Track density pixels", 0, unscaleddensitysubset);
     			setResult("Track density", 0, densitysubset);
     		}
     		else {	
			setResult("Track density pixels", 0, unscaleddensity);
			setResult("Track density", 0, density);
		}	
     		setResult("Number events", 0, nEvents);
     		setResult("Time interval (s)", 0, timeinterval);
     		setResult("Microns per pixel", 0, micronsperpixel);
     		setResult("Channel", 0, curSlice);
     		setResult("Kymograph width", 0, kymographWidth);
     		setResult("Subset or kymograph", 0, useSubset); 
     		
     		if(useSubset>0){
     			saveAs("measurements", dir+channelImageTitle+"eventDensitySubset.xls");
     			saveAs("measurements", dir+channelImageTitle+"eventDensitySubset.txt");	
     		}
     		else {	
			saveAs("measurements", dir+channelImageTitle+"eventDensity.xls");
			saveAs("measurements", dir+channelImageTitle+"eventDensity.txt");
     		}
	}
}
	

