%% Prepare data tables for R
% Script to process file structure output from pairSyanpses.m to
% make datatables ready for R and further analysis for
% "Differential nanoscale organization of excitatory synapses onto excitatory vs. inhibitory neurons"
% Poorna A Dharmasri, Aaron D Levy, Thomas A Blanpied, PNAS 2024
% authors: Aaron Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine


% identify week folders
f1 = pwd;
files = dir(f1);
dirflags = [files.isdir];
weekfolders = files(dirflags);
weekfolders = {weekfolders(3:end).name};

for subf = 1:length(weekfolders) % begin loop over week folders

    % identify results folders 
    thissubf = fullfile(f1,weekfolders{subf});
    subfiles = dir(thissubf);
    subfflags = [subfiles.isdir];
    resultsfolders = subfiles(subfflags);
    resultsfolders = {resultsfolders(3:end).name};

    fieldcount = 1;
    for resultsf = 1:length(resultsfolders) % begin loop over results subfolders

        % Load the data and trim the structure
        thisfolder = fullfile(thissubf,resultsfolders{resultsf});
        if exist(fullfile(thisfolder,'ROI_key.txt'),'file') == 0
            continue % skip the loop if an ROI key doesn't exist, ie that folder was not processed
        end
        roikey = readmatrix(fullfile(thisfolder,'ROI_key.txt'));

        postfile = dir(fullfile(thisfolder,'PostResults*.csv'));
        postsyn = readtable(fullfile(thisfolder,postfile.name));
        postsyn.Properties.VariableNames{1} = 'ROInum';
        newNames = append(postsyn.Properties.VariableNames([1:9,13]),'_post');
        postsyn = renamevars(postsyn,[1:9,13],newNames);
        postsyn = postsyn(:,[1:8,12,13,15,16]);

        prefile = dir(fullfile(thisfolder,'PreResults*.csv'));
        presyn = readtable(fullfile(thisfolder,prefile.name));
        presyn.Properties.VariableNames{1} = 'ROInum';
        newNames = append(presyn.Properties.VariableNames([1:8,15]),'_pre');
        presyn = renamevars(presyn,[1:8,15],newNames);
        presyn = presyn(:,[1:8,15]);

        % Pair them and add a synapse number column
        T = table();
        for s = 1:length(roikey)

            thispair = roikey(s,:);
            if ~isnan(thispair(2))

                thispost = postsyn(postsyn.ROInum_post==thispair(1),:);
                thispre = presyn(presyn.ROInum_pre==thispair(2),:);
                T = [T; table(s,'VariableNames',{'synapse_num'}) thispost(:,[1:8,11]) thispre thispost(:,[9,10,12])];

            end

        end

        [~,week,~] = fileparts(thissubf);
        week = replace(week,' ','');
        if contains(thisfolder,{'488-pFCaGW','pFCaGW-488'})==1
            cond = 'pCaMKII';
        elseif contains(thisfolder,{'488-PV','PV-488'})==1
            cond = 'PV';
        end

        writetable(T,fullfile(thisfolder,'pairedResults.csv')); % save into the results folder
        writetable(T,fullfile(f1,[week '_' cond '_paired_' num2str(fieldcount) '.csv'])) % save into the parent folder

        fieldcount = fieldcount + 1;

    end % results loop

end % week loop

%% merge
csvlist = dir('*.csv');
T = table();

for thiscsv = 1:length(csvlist)
    
    temp = readtable(fullfile(csvlist(thiscsv).folder,csvlist(thiscsv).name));
    T = [T;temp];

end

writetable(T,fullfile(csvlist(1).folder,'mergedPairedData.csv'));



