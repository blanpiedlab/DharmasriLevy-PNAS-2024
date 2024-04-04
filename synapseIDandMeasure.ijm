//Script to process folders of tifs for synapse identification and measurement by SynQuant used in
//"Differential nanoscale organization of excitatory synapses onto excitatory vs. inhibitory neurons"
//Poorna A Dharmasri, Aaron D Levy, Thomas A Blanpied, PNAS 2024
//author: Aaron Levy
//copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine


//This script can process a one-level tree of folders full of tif files.
//ie, directory/subdirectory1; directory/subdirectory2 through directory/subdirectory(n)

//Check for needed plugins and check for mesaurement settings
List.setCommands;
if (List.get("Read and Write Excel")=="") {exit("You need to install the Read and Write Excel plugin. Go to Help>Update>Manage Update Sites, add ResultsToExcel, restart FIJI.");}
if (List.get("SynQuantVid ")=="") {exit("You need to install the SynQuant plugin. Go to https://github.com/yu-lab-vt/SynQuant and follow the directions, and restart FIJI");}
run("Set Measurements...");
print("\\Clear"); //clear the log

//Get path info 
path = getDirectory("Input directory"); //prompts user to input a base folder to analyze ("directory" above)
path = replace(path,"/",File.separator); //replace any / with \ in path name

//Loop over folders and threshold them
setBatchMode(false);
list = getFileList(path);
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], "/")) {
		folder = replace(list[i],"/",File.separator);
		thresholdfolder(path,folder); 
	}
}

//Now go back and do the synquant
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], "/")) {
		folder = replace(list[i],"/",File.separator);
		sqfolder(path,folder);
	}
}

function thresholdfolder(path,folder) {
//First loop through all the tifs to make thresholds
//If you want to skip the thresholding step (ie already have thresholds with appropriate binary files) then comment out this entire loop
setBatchMode(false);
list = getFileList(path+folder); 
for (i=0; i<list.length; i++) {
    if (endsWith(list[i], ".tif")) {//loop will find all tifs
    	filename = list[i];
    	print("Thresholding " + folder + filename);
    	savefolder = File.getNameWithoutExtension(path + list[i]) + "_Results" + File.separator;
        processbinaries(path,folder,filename,savefolder); //runs function processfiles on the file
	}
}
}

function sqfolder(path,folder) {
//Now loop through all the tifs to do the synaptic quantification
count = 1;
list = getFileList(path+folder);
for (i=0; i<list.length; i++) {
    if (endsWith(list[i], ".tif")) {//loop will find all tifs again
    	filename = list[i];
    	print("Now processing " + folder + filename);
    	savefolder = File.getNameWithoutExtension(path + folder+ filename) + "_Results" + File.separator;
        count = processfiles(path,folder,filename,savefolder,count); //runs function processfiles on the file

	}
}
}
print("Analysis is done");
beep();
setBatchMode(false);

function processbinaries(root,folder,filename,savefolder) { 

	if (File.exists(root + folder + savefolder + "binarybase.tif")) {
		continue;
	}
	path = root + folder;
	open(path + savefolder + "registered.tif");
	
	title = getTitle();
	waitForUser("If you want to crop, do it now, then hit OK");
	saveAs("Tiff",path + savefolder + "imstack_crop_nobkgsub.tif");
	rename(title);
	
	//Split channels and save background subtracted images
	print("     Background subtracting");
	getDimensions(width, height, channels, slices, frames);
	if (channels != 3) {exit("You need to load a 3 channel image for this to work");}
	
	run("Split Channels");
	gfp = "C1-" + title;
	post = "C2-" + title;
	pre = "C3-" + title;
	ims = newArray(gfp,post,pre);
	for (i = 0; i < ims.length; i++) {
		cutdecimal = 0.01; //decimal value of percentage to subtract, ie 0.01 will subtract the pixel value that 1% of pixels fall under
		bkg_sub_pct(path, savefolder, ims[i], cutdecimal); //use this option to bkg sub based on pixel intensity below cutdecimal %
	}

	
	//Make a dendrite mask with user input threshold
	print("     Making dendrite mask");
	open(path + savefolder + gfp);
	setAutoThreshold("Default dark"); //pops the threshold window if not already open
	waitForUser("Set a threshold on the GFP image");
	run("Convert to Mask");
	process = getNumber("Process this image? 1 for yes 2 for no", 1);
	if (process != 1) {
		run("Close All");
		return;  //end the loop if user chooses not to process this image
	}
	else {	
		saveAs("tif",path + savefolder + "binarybase.tif");
		cleanup();
	}
	
}

function processfiles(root,folder,filename,savefolder,count) {
	path = root + folder;
	//skip this folder if the binary doesn't exist
	if (File.exists(path + savefolder + "binarybase.tif") != 1) {print("     Skipping folder " + savefolder + " due to no binary"); count++; return count;}

	//note that these are setup for registered images, change as needed
	gfp = "C1-" + "registered.tif";
	post = "C2-" + "registered.tif";
	pre = "C3-" + "registered.tif";
	ims = newArray(gfp,post,pre);
	
	open(path + savefolder + gfp);
	setBatchMode(true); //runs slightly faster in batch
	dend_mask(path, savefolder, gfp); 
	cropbinary(path, savefolder); //Crop binary to just the ROI for later use
	setBatchMode(false);
	cleanup();

	

	
	//ID synapses in both channels using SynQuant
	im2 = newArray(post,pre);
	pp = newArray("Post","Pre");
	for (i = 0; i < im2.length; i++) {
		
		close("ROI manager"); //SynQuant really doesn't tolerate having ROI manager open, so make sure it's closed
		//Find puncta with synquant and save
		print("     Running SynQuant on " + pp[i] + " image");
		open(path + savefolder + im2[i]);
		minfill = 0.65;
		whratio = 6;
		run("SynQuantVid ", "z-score=10 min=10 max=200 min_0=" + minfill + " max_0=" + whratio + " post-synapse=" + im2[i] + " pre-synapse=Null way=Null dendrite=Null extended=0 z=1 zscore=10");
		roiManager("Save", path + savefolder + "SynQuant" + pp[i] + ".zip");
		run("Close"); //needed to nix the synquant results window, which is annoying
		cleanup();
		
	}
	
	//Filter postsynaptic ROIs by GFP mask
	print("     Finding postsynaptic ROIs that intersect GFP mask");
	imfile = "clean_binary.tif";
	roifile = "SynQuantPost.zip";
	savefile = "Post_FilteredByGFP.zip";
	mask_intersect(path, savefolder, imfile, roifile, savefile);
	cleanup();
	
	//Threshold filtered post ROIs to 40% max
	print("     Thresholding filtered postsynaptic ROIs");
	setBatchMode(true); //batch mode is critical for these steps for speed
	//first dilate the ROIs by dx pixels (2 seems ideal)
	imfile = post;
	roifile = "Post_FilteredByGFP.zip";
	savefile = "Post_FilteredByGFPDilated.zip";
	dx = 2;
	resize_roi(path, savefolder, imfile, roifile, savefile, dx);
	cleanup();
	//Now threshold within that dilated ROI to 40% of max pixel intensity in ROI
	roifile = "Post_FilteredByGFPDilated.zip";
	savefile = "Post_FilteredByGFPThresholded.zip";
	threshold_roi(path, savefolder, imfile, roifile, savefile);
	setBatchMode(false);
	cleanup();
	
	//Filter presynaptic ROIs by thresholded, gfp-overlapping postsynaptic ROIs
	print("     Finding presynaptic ROIs that intersect postsynaptic mask");
	//Make postynaptic binary
	roifile = "Post_FilteredByGFPThresholded.zip";
	savefile = "Post_FilteredByGFPThresholded_Binary.tif";
	fill_ROIs(path, savefolder, roifile, savefile);
	cleanup();
	//Filter the presynaptic ROIs to those that overlap the post binary
	roifile = "SynQuantPre.zip";
	imfile = "Post_FilteredByGFPThresholded_Binary.tif";
	savefile = "Pre_FilteredByPost.zip";
	mask_intersect(path, savefolder, imfile, roifile, savefile);
	cleanup();
	
	//Repeat the dilation and thresholding on the presynaptic side
	print("     Thresholding filtered presynaptic ROIs");
	setBatchMode(true); //batch mode is critical for these steps for speed
	//first dilate the ROIs by dx pixels (2 seems ideal)
	imfile = pre;
	roifile = "Pre_FilteredByPost.zip";
	savefile = "Pre_FilteredByPostDilated.zip";
	dx = 2;
	resize_roi(path, savefolder, imfile, roifile, savefile, dx);
	cleanup();
	//Now threshold within that dilated ROI to 40% of max pixel intensity in ROI
	roifile = "Pre_FilteredByPostDilated.zip";
	savefile = "Pre_FilteredByPostFINAL.zip";
	threshold_roi(path, savefolder, imfile, roifile, savefile);
	setBatchMode(false);
	cleanup();
	
	//Filter postsynaptic ROIS by the filtered presynaptic ROIs
	print("     Finding postynaptic ROIs that intersect the presynaptic mask");
	//Make presynaptic binary
	roifile = "Pre_FilteredByPostFINAL.zip";
	savefile = "Pre_FilteredByPostThresholded_Binary.tif";
	fill_ROIs(path, savefolder, roifile, savefile);
	cleanup();
	//Filter postsynaptic ROIs by those that overlap the pre binary
	roifile = "Post_FilteredByGFPThresholded.zip";
	imfile = "Pre_FilteredByPostThresholded_Binary.tif";
	savefile = "Post_FilteredByPreFINAL.zip";
	mask_intersect(path, savefolder, imfile, roifile, savefile);
	cleanup();
		
	//Measure intensity in puncta
	print("     Measuring intensity");
	p1 = newArray("Post","Pre");
	p2 = newArray("Pre","Post");
	cond = getConditions(filename);
	for (i = 0; i < im2.length; i++) {
		
		open(path + savefolder + im2[i]);
		roiManager("open", path + savefolder + p1[i] + "_FilteredBy" + p2[i] + "FINAL.zip");
		roiManager("multi-measure");
		for (row = 0; row < nResults; row++) {
			//setResult("Image Count", row, count);
			setResult("Week",row,substring(folder,0,indexOf(folder,File.separator)));
			setResult("Condition",row,cond[0]);
			setResult("PrePost", row, p1[i]);
			if (p1[i] == "Post") {
				setResult("Stain",row,"PSD95");
			}
			else if (p1[i] == "Pre") {
				setResult("Stain",row,cond[1]);
			}
			setResult("Image Name", row, filename);
			
		}
		updateResults();
		saveAs("Results", path + savefolder + p1[i] + "Results_Count" + count + ".csv");

		week = substring(folder,0,indexOf(folder,File.separator));
		week = replace(week," ","");
		saveAs("Results", root + week + "_" + cond[0] + "_" + p1[i] + "_" + cond[1] + "_" + count + ".csv");
		//run("Read and Write Excel","file=[" + input + File.separator + "results.xlsx] dataset_label=[" + pp[i] + "] sheet=[" + pp[i] + filecount + "]");
	
		cleanup();
	
	}

	count++;
	return count;
}





//subfunctions



function bkg_sub_pct(path, savefolder, im, cutdecimal) {

	//Background subtract finds the pixel intensity value in which threshold number of 
	//pixels have values below it (ie, if cutdecimal = 0.01, we're finding the 
	//pixel intensity at which 1% of pixels have a value below that, and subtracting that value
	selectWindow(im);
	getDimensions(width, height, channels, slices, frames);
	totalpx = width * height;
	threshold = cutdecimal * totalpx;
	getHistogram(0, counts, Math.pow(2,bitDepth()));
	cutoff = 0;
	i = 0;
	do {
		cutoff = cutoff + counts[i];
		i++;
	} while (cutoff < threshold);
	run("Subtract...", "value=" + i);
	saveAs("tiff", path + savefolder + im);
	close(im);

}

function dend_mask(path, savefolder, im) {
	//makes a binary dendrite selection and then deletes out any particles 
	//smaller than 10 units, resulting in main dendrite selected.
	
	//Binary process
	selectWindow(im);
	run("Convert to Mask");
	run("Dilate");
	run("Dilate");
	run("Close-");

	//Select the binary
	run("Create Selection");
	roiManager("Add");

	//Remove parts of the binary that are smaller than 10
	run("Analyze Particles...", "size=0-10 add"); //Find list of ROIs with small size
	n = roiManager("count");
	origrois = Array.getSequence(n);
	roiManager("select", origrois);
	roiManager("XOR"); //This will delete overlapping ROIs, so small ROIs will be filtered out of the selection
	roiManager("Add");
	roiManager("select",origrois);
	roiManager("delete");
	saveAs("Tiff", path + savefolder + "binary.tif"); //save the tif in the new folder as project_i
	roiManager("Save", path + savefolder + "ImMask.zip");
	cleanup();

}

function cropbinary(path, savefolder) {
	
	//Cleans up the binary to just what's in the saved mask
	open(path + savefolder + "binary.tif");
	roiManager("open", path + savefolder + "ImMask.zip");
	roiManager("select",0);
	run("Colors...", "foreground=black background=white selection=yellow");
	run("Clear Outside");
	roiManager("deselect");
	saveAs("tiff", path + savefolder + "clean_binary.tif");
	cleanup();
}

function resize_roi(path, savefolder, imfile, roifile, savefile, dx) {

	//Increase size of ROIs by dx in all directions
	//Uncomment 4 lines to instead of dialog ask for user input of dx
	open(path + savefolder + imfile);
	roiManager("open", path + savefolder + roifile);
	//dx = 2;
	//Dialog.create("Adjust size of your ROIs");
	//Dialog.addNumber("Increase by:", dx);
	//Dialog.show();
	//dx = Dialog.getNumber();
	n = roiManager("count");
	if (n==0)
		exit("The ROI Manager is empty");
	for (i=0; i<n; i++) {
		roiManager("select", i);
		run("Enlarge...", "enlarge="+ dx + " pixel");
		roiManager("update");
	}

	roiManager("Save", path + savefolder + savefile);
      
}

function threshold_roi(path, savefolder, imfile, roifile, savefile) {

	//Thresholds within each ROI to 40% of max intensity, like metamorph
	open(path + savefolder + imfile);
	roiManager("open", path + savefolder + roifile);
	img = getTitle();
	n = roiManager("count");
	origrois = Array.getSequence(n);
	setBackgroundColor(0,0,0);
	
	for (i = 0; i < n; i++) {

		showProgress(i, n);
		selectWindow(img);
		run("Select All");
		run("Duplicate...", "title=temp");
		roiManager("select", i);
		run("Clear Outside");
		getStatistics(area, mean, min, max, std, histogram);
		threshval = (max-min) * 0.4 + min;
		setThreshold(threshval, Math.pow(2,bitDepth()));
		run("Convert to Mask");
		run("Fill Holes");
		run("Create Selection");
		roiManager("add");
		
		if (selectionType() == 9) { //This section will crop any composite ROIs down to just the largest component		
		
			roiManager("split");
			roiManager("select", n+i);
			roiManager("delete");
			size = 0;
			for (j = n+i; j < roiManager("count"); j++) {
				roiManager("select",j);
				if (Roi.size > size) {
					index = j;
					size = Roi.size;
				}
			}
			inds = newArray(0);
			for (j = n+i; j < roiManager("count"); j++) {
				if (index == j) {continue;}
				inds = Array.concat(inds,j);
			}
			roiManager("select",inds);
			roiManager("delete");
			
		} 
	
		close("temp");

	}

	roiManager("select", origrois);
	roiManager("delete");
	roiManager("save", path + savefolder + savefile);

}


function mask_intersect(path, savefolder, imfile, roifile, savefile) {
	
	//Measure whether there is any signal from binary mask in each ROI and keeps those with signal
	open(path + savefolder + imfile);
	roiManager("open", path + savefolder + roifile);
	roiManager("multi-measure measure_all");
	n = roiManager("count");
	tosave = newArray(0);
	
	for (i = 0; i <n; i++) {

		val = getResult("Mean",i);
		if (val > 0) {
			tosave = Array.concat(tosave,i);
		}
			
	}

	roiManager("select",tosave);
	roiManager("save selected", path + savefolder + savefile);

}

function fill_ROIs(path,savefolder,roifile,savefile) {
	
	//Makes a binary image from a bunch of synapse ROIs
	open(path + savefolder + "binary.tif");
	getDimensions(width, height, channels, slices, frames);
	close("binary.tif");
	newImage("Untitled", "16-bit white", width, height, 1);
	roiManager("Open", path + savefolder + roifile);
	roiManager("select", Array.getSequence(roiManager("count")));
	roiManager("Fill");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	roiManager("reset");
	saveAs("Tiff", path + savefolder + savefile);

}

function cleanup() {
	
	//Close any image windows and results windows, reset the ROI manager
	list = getList("image.titles");
	if (list.length>0) {
		run("Close All");
	}
	if (isOpen("Results")) {
		close("Results");
	}
	roiManager("reset");
}

function getConditions(filename) {
	
		isCam = indexOf(filename,"488-pFCaGW");
		if (isCam==-1) {isCam = indexOf(filename,"pFCaGW-488"); }
		isPV = indexOf(filename,"488-PV");
		if (isPV==-1) {isPV = indexOf(filename,"PV-488");}
		isMunc = indexOf(filename,"Munc");
		isRIM = indexOf(filename,"RIM");
		if (isCam > 1 && isPV == -1) {
			expv = "pCaMKII";
		}
		else if (isCam ==-1 && isPV >1) {
			expv = "PV";
		}
		if (isMunc > 1 && isRIM == -1) {
			prestain = "Munc";
		}
		else if (isMunc ==-1 && isRIM >1) {
			prestain = "RIM";
		}

		return newArray(expv,prestain);
}

