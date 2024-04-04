//Script to process folders of tifs for chromatic aberration correction used in
//"Differential nanoscale organization of excitatory synapses onto excitatory vs. inhibitory neurons"
//Poorna A Dharmasri, Aaron D Levy, Thomas A Blanpied, PNAS 2024
//author: Aaron Levy
//copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine

//File choosing setup
print("\\Clear"); //clear the log
Dialog.create("Check the pixel conversion");
Dialog.addNumber("pixel distance (pixels per micron)", 9.9508, 4, 6, "px/um");
Dialog.show();
px = Dialog.getNumber();
path = getDirectory("Input directory"); //grab the directory
path = replace(path,"/",File.separator); //replace any / with \ in path name
transmaskpath = File.openDialog("Select the translation mask file");

//Structure to loop through folders within parent directory.
setBatchMode(true);
list = getFileList(path);
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], "/")) {
		folder = replace(list[i],"/",File.separator);
		print("Processing " + folder);
		processfiles(path,folder,px); 
	}
}
setBatchMode(false);

//Function to process each file within the subfolders
function processfiles(path,folder,px) {
	list = getFileList(path+folder);
	for (i=0; i<list.length; i++) {
	    if (endsWith(list[i], ".tif")) {//loop will find all tifs 
	    	filename = list[i];
	    	print("     Processing " + filename);
	    	savefolder = File.getNameWithoutExtension(path + folder + filename) + "_Results" + File.separator; //generate a savefolder for the registered tif
	    	if (!File.exists(path + folder + savefolder)) { // make the folder if it's not there
	    		File.makeDirectory(path + folder + savefolder); 
	    	}
	      
	      	//do the actual chromatic abberation correction
			open(path + folder + filename);
			origim = getTitle();
			origim = replace(origim,".tif","");
			run("Hyperstack to Stack"); //nanoj requires stack not hyperstack
			run("Register Channels - Apply", "open=[" + transmaskpath + "] don't");
			
			//turn stack back into hyperstack and reset LUTs/contrast; adjust for your own dat if needed
			//This is hardcoded for the PV stuff
			selectWindow(origim + " - Registered");
			run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Color");
			run("Conversions...", " "); //nanoj converts to 16bit and you need to go back, this is required for keeping px intensities right
			run("16-bit");
			run("Conversions...", "scale"); //turn it on for later
			run("Set Scale...", "distance=" + px + " known=1 unit=micron"); //nanoj also gets rid of the scale, adjust distance as appropriate
			Property.set("CompositeProjection", "Max");
			Stack.setDisplayMode("composite");
			run("OPF fresh");
			run("Enhance Contrast", "saturated=0.35");
			run("Next Slice [>]");
			run("Magenta");
			run("Enhance Contrast", "saturated=0.35");
			run("Next Slice [>]");
			run("Yellow");
			run("Enhance Contrast", "saturated=0.35");
			
			//save "registered.tif" into the savefolder within each week
			saveAs("tif",path + folder + savefolder + "registered.tif");
			run("Close All");
			
	    }
	}
}

