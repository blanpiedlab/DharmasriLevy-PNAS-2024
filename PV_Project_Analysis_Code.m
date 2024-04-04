% "Differential nanoscale organization of excitatory synapses onto excitatory vs. inhibitory neurons"
% Poorna A Dharmasri, Aaron D Levy, Thomas A Blanpied, PNAS 2024
% authors: Poorna A Dharmasri
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine

%% SMAP data conversion
clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Master')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   for v = 1:size(files,1)
                       if contains(files(v).name, 'Cali')
                          [~,~] = smaploader(files(v).name, 'convertToPixels',1, 'savehdf',1); 
                       elseif ~contains(files(v).name, 'Cali')
                           [~,~] = smaploader(files(v).name, 'convertToPixels',1, 'savehdf',1);
                           fname = strrep(files(v).name,'.mat','.hdf5');
                           [hs,loc] = Omniloader(fname);
                           col = getColumns(hs);
                           loc(loc(:,col.lpx)>0.5,:) = [];
                           loc(loc(:,col.lpy)>0.5,:) = [];
                           fname = fname(1:end-5);
                           RPTP(loc,hs,fname,'_filter');
                       end
                   end
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Master')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% Drift Correct region files only
clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Region')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   for v = 1:size(files,1)
                       if contains(files(v).name,'filter.hdf5')
                          picassoUndrift([pwd,'\',files(v).name]);
                       end
                   end
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Region')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% make duv hdf5
clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Cali')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   fnnum = 1; 
                   filestocombo = {};
                   for fn = 1:size(files,1)
                       if contains(files(fn).name,'hdf5')
                           filestocombo(fnnum) = {files(fn).name};
                           fnnum = fnnum+1;
                       end
                   end
                   concatduv = [];
                   for v = 1:size(filestocombo,2)
                       
                      [hs,loc] = Omniloader(filestocombo{v});
                      col = getColumns(hs);
                      if v == 1
                          concatduv = [concatduv;loc];
                      elseif v>1
                          loc(:,col.frame) = loc(:,col.frame)+((v-1)*100);
                          concatduv = [concatduv;loc];
                      end
                       
                       
                   end
                   fname = filestocombo{v}(1:end-5);
                   if contains(targetfolder(iv).name,'561')
                        RPTP(concatduv,hs,fname,'null','changename','DUVCali561');
                   elseif contains(targetfolder(iv).name,'488')
                       RPTP(concatduv,hs,fname,'null','changename','DUVCali488');
                   end
                       
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Cali')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% make duv tforms
clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Cali')
                   round = cd(targetfolder(iv).name);
                   
                   if contains(targetfolder(iv).name,'488')
                      ADLmakeTform_poly_bidirectional_test('DUVCali488.hdf5',100,'l2r', 'SplitPoint', 804); 
                   elseif contains(targetfolder(iv).name,'561')
                      ADLmakeTform_poly_bidirectional_test('DUVCali561.hdf5',100,'l2r', 'SplitPoint', 804);
                   end
                   
                   
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Cali')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% apply duv

clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Region')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   for v = 1:size(files,1)
                       if contains(files(v).name,'_undrift.hdf5')
                          fname = files(v).name;
                           [header,all_locs] = ADLcorrectDUV_poly_bidirectional(files(v).name,100,'SplitPoint',804,'dontsave',1);
                           col = getColumns(header);
                           all_locs(:,col.x) = all_locs(:,col.x)-804;
                           fname = fname(1:end-5); 
                          RPTP(all_locs,header,fname,'_DUV','widthupdate',804);
                          
                       end
                   end
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Region')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% split into channels

clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Region')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   for v = 1:size(files,1)
                       if contains(files(v).name,'_DUV.hdf5')
                          fname = files(v).name;
                          if contains(fname,'PV')
                            snl = 'PV_Munc';
                            snr = 'PV_PSD95';
                          elseif contains(fname,'CaMKII')
                              snl = 'CaMKII_Munc';
                              snr = 'CaMKII_PSD95';
                          end
                           [hs,loc] = Omniloader(fname, 'verbose',0);
                            col = getColumns(hs);
                            f1l = loc(loc(:,col.channel)==1,:);
                            f1r = loc(loc(:,col.channel)==2,:);
                            f1 = files(v).name(1:end-5);
                            RPTP(f1l,hs,f1,'null','changename',snl);
                            RPTP(f1r,hs,f1,'null','changename',snr);
                            
                          
                          
                       end   
                   end
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Region')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% xcorr and filter

files = dir();
files = files(3:end,:);
filestorun = {};
filecount = 1;
for v = 1:size(files,1)
    if contains(files(v).name,'Munc.hdf5')
        filestorun{filecount,1} = files(v).name;
        filestorun{filecount,2} = strrep(files(v).name,'Munc','PSD95');
        filecount = filecount+1;
    end
end
for vi = 1:size(filestorun,1)
    fname1 = filestorun{vi,1};
    fname2 = filestorun{vi,2};
    [hs,locs1] = Omniloader(fname1,'verbose',0);
    [~,locs2] = Omniloader(fname2,'verbose',0);
    col = getColumns(hs); %For use in filtering step below
    
    % the last variable "distance" sets the min/max allowable shift in pixels
    shiftV = PAINTshift(locs1(:,[col.x col.y]), locs2(:,[col.x col.y]), 5, 100, 50, [], [0 24]);
    locs2(:,col.x) = locs2(:,col.x)+shiftV(1);
    locs2(:,col.y) = locs2(:,col.y)+shiftV(2);
    
    ph1=histogram(locs1(:,col.phot),'BinWidth',100,'BinCountsMode','auto');
    locs1photmode=ph1.BinEdges(ph1.BinCounts==max(ph1.BinCounts));
    
    ph2=histogram(locs2(:,col.phot),'BinWidth',100,'BinCountsMode','auto');
    locs2photmode=ph2.BinEdges(ph2.BinCounts==max(ph2.BinCounts));
    
    
    pixelsize = 100;
    lpcutoff = 20/pixelsize;
    
    locs1(locs1(:,col.phot)<locs1photmode,:)= [];
    locs1(locs1(:,col.sigmax)>2,:)= [];
    locs1(locs1(:,col.sigmay)>2,:)= [];
    locs1(locs1(:,col.lpx)>lpcutoff,:) = [];
    locs1(locs1(:,col.lpy)>lpcutoff,:) = [];
    figure
    histogram(locs1(:,col.LLrel));
    LLrelcut = input('Cutoff?');
    locs1(locs1(:,col.LLrel)<LLrelcut,:)=[];
    close all
    
    locs2(locs2(:,col.phot)<locs2photmode,:)= [];
    locs2(locs2(:,col.sigmax)>2,:)= [];
    locs2(locs2(:,col.sigmay)>2,:)= [];
    locs2(locs2(:,col.lpx)>lpcutoff,:) = [];
    locs2(locs2(:,col.lpy)>lpcutoff,:) = [];
    figure
    histogram(locs2(:,col.LLrel));
    LLrelcut = input('Cutoff?');
    locs2(locs2(:,col.LLrel)<LLrelcut,:)=[];
    close all
    
    fname1 = fname1(1:end-5);
    fname2 = fname2(1:end-5);
    
    RPTP(locs1, hs, fname1, '_xcorr_filter');
    RPTP(locs2, hs, fname2, '_xcorr_filter');
end
%% link

clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Region')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                  for v = 1:size(files,1)
                     if contains(files(v).name,'xcorr')
                         picassoLink([pwd,'\',files(v).name],.3,5);
                     end
                  end
                                     
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Region')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% get rid of extreme z noise

clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Region')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   filestorun = {};
                   filecount = 1;
                   for v = 1:size(files,1)
                      if contains(files(v).name,'Munc') && contains(files(v).name,'link.hdf5')
                          filestorun{filecount,1} = files(v).name;
                          filestorun{filecount,2} = strrep(files(v).name,'Munc','PSD95');
                          filecount = filecount+1;
                      end
                   end
                   for vi = 1:size(filestorun,1)
                       fname1 = filestorun{vi,1};
                       fname2 = filestorun{vi,2};
                        [hs,locs1] = Omniloader(fname1,'verbose',0);
                        [~,locs2] = Omniloader(fname2,'verbose',0);
                        col = getColumns(hs); %For use in filtering step below
                        
                        figure
                        axis equal
                        hold on
                        scatter(locs1(:,col.x), locs1(:,col.z), '.b');
                        scatter(locs2(:,col.x), locs2(:,col.z), '.r');
                        
                        lowerbound = input('what is lower bound?');
                        upperbound = input('what is upper bound?');
                        
                        close all
                        
                        locs1(locs1(:,col.z)>upperbound | locs1(:,col.z)<lowerbound,:)=[];
                        locs2(locs2(:,col.z)>upperbound | locs2(:,col.z)<lowerbound,:)=[];
                        
                        fname1 = fname1(1:end-5);
                        fname2 = fname2(1:end-5);
                        RPTP(locs1,hs,fname1,'_zcut');
                        RPTP(locs2,hs,fname2,'_zcut');

                            
                   end                  
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Region')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% dbscan

clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Region')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   for v = 1:size(files,1)
                      if contains(files(v).name,'zcut.hdf5')
                         picassoDbscan([pwd,'\',files(v).name],.3,5); 
                      end
                   end        
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Region')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% Cluster Filter
clear
clc
tic
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
           week = cd(roundir(iii).name);
           targetfolder = dir();
           targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               if contains(targetfolder(iv).name,'Region')
                   round = cd(targetfolder(iv).name);
                   files = dir();
                   files = files(3:end,:);
                   for v = 1:size(files,1)
                      if contains(files(v).name,'dbscan.hdf5')
                        dbscanfile = files(v).name;
                        clusterfile = strrep(files(v).name,'dbscan','dbclusters');
                        std_range = [2500 20000];    
                        [in,out,dbs_hs] = clusterFilter(dbscanfile,clusterfile,std_range,'minn',15,'spotcheck',0);
                        RPTP(in,dbs_hs,dbscanfile(1:end-5),'_clusterFilter'); 
                          
                      end
                   end        
                  cd(round); 
               elseif ~contains(targetfolder(iv).name,'Region')
                   continue
               end
           end
            cd(week); 
        end
        cd(project);
    end
    cd(main)
end
disp('Done!')
toc
%% 

clear
clc
tic
fincount = 1; 

datasetstrings{1} = {'CaMKII';'PV'};
% filedir = dir();
% filedir = filedir(3:end,:);
% % for i = 1:size(filedir,1)
% %     main = cd(filedir(i).name);
%     weekdir = dir();
%     weekdir = weekdir(3:end,:);
%     for ii = 1:size(weekdir,1)
%         project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
            week = cd(roundir(iii).name);
            targetfolder = dir();
            targetfolder = targetfolder(3:end,:);
           for iv = 1:size(targetfolder,1)
               round = cd(targetfolder(iv).name);
                dataset = 1; 
                for xi = 1:size(datasetstrings{dataset},1)
                    if iii == 3 && xi == 1
                        continue
                    end
                    if iii == 3 && xi == 3
                        continue
                    end
                    regionfileMunc = [datasetstrings{dataset}{xi},'_Munc_filter_xcorr_link_zcut_dbscan_clusterFilter.hdf5'];
                    unfilteredMunc = [datasetstrings{dataset}{xi},'_Munc_filter_xcorr_link_zcut_dbscan.hdf5'];
                    regionfilePSD = [datasetstrings{dataset}{xi},'_PSD95_filter_xcorr_link_zcut_dbscan_clusterFilter.hdf5'];
                    
                    folderdir = dir();
                    folderdir = folderdir(3:end,:);
                    for v = 1:size(folderdir,1)
                        if contains(folderdir(v).name,datasetstrings{dataset}(xi))
                            main = cd(folderdir(v).name);
                            t = load('tforml2r.mat');
                            fname = dir('*.tif');
                            cellfill = imread(fname.name);
                            cellfill = cellfill(:,1:804);
                            cft = imwarp(cellfill,t.tform);
                            
                            
                            pxvals = unique(sort(cft(:),'Descend'));
                            pxcutoff = pxvals(floor(.3*size(pxvals,1)));
                            mask = cft;
                            mask(mask<pxcutoff) = 0;
                            figure
                            axis equal
                            hold on
                            image(cft,'CDataMapping','scaled');
                            colormap('gray');
                        end
                    end
                    cd(main)
                    for vi = 1:size(folderdir,1)
                        if contains(folderdir(vi).name,'Region')
                            main = cd(folderdir(vi).name);
                            [hs,munc,col] = Omniloader(regionfileMunc,'verbose',0);
                            [~,psd] = Omniloader(regionfilePSD,'verbose',0);
                           
                        end
                    end
                    
                    psdgroups = unique(psd(:,col.groups));
                    muncgroups = unique(munc(:,col.groups));
                    muncbygroups = {};
                    for vii = 1:size(muncgroups,1)
                        muncbygroups{vii,1} = munc(munc(:,col.groups)==muncgroups(vii),:);
                    end
                    psdrec = [];
                    checktestcell = [];
                    for viii = 1:size(psdgroups,1)
                        thispsd = psd(psd(:,col.groups)==psdgroups(viii),[col.x col.y col.z]);
                        psdshape = alphaShape(thispsd(:,1), thispsd(:,2), thispsd(:,3));
                        testcell = cellfun(@(x) inShape(psdshape,x(:,col.x),x(:,col.y),x(:,col.z)),muncbygroups,'UniformOutput',0);
                        for ix = 1:size(testcell,1)
                            checktestcell(ix,:) = sum(unique(testcell{ix}));
                        end
                        psdrec = [psdrec;sum(unique(checktestcell))];
                    end
                    psdrec = logical(psdrec);
                    keptpsdgroups = psdgroups(psdrec);
                    keptpsdgroupspercell = {};
                    
                    for keptgroup = 1:size(keptpsdgroups,1)
                        keptpsdgroupspercell{keptgroup,1} = psd(psd(:,col.groups)==keptpsdgroups(keptgroup),:);
                        scatter3(psd(psd(:,col.groups)==keptpsdgroups(keptgroup),col.x), psd(psd(:,col.groups)==keptpsdgroups(keptgroup),col.y), psd(psd(:,col.groups)==keptpsdgroups(keptgroup),col.z)+10,'.m');
                    end
                    muncrec = [];
                    checktestcell2 = [];
                    for muncgroup = 1:size(muncbygroups,1)
                        thismunc = muncbygroups{muncgroup}(:,[col.x col.y col.z]);
                        muncshape = alphaShape(thismunc(:,1), thismunc(:,2), thismunc(:,3));
                        testcell2 = cellfun(@(x) inShape(muncshape,x(:,col.x),x(:,col.y),x(:,col.z)),keptpsdgroupspercell,'UniformOutput',0);
                        for tcrow = 1:size(testcell2,1)
                            checktestcell2(tcrow,:) = sum(unique(testcell2{tcrow}));
                        end
                        muncrec = [muncrec;sum(unique(checktestcell2))];
                    end
                    muncrec = logical(muncrec);
                    keptmuncgroups = muncgroups(muncrec);
                    keptmuncgroupspercell = {};
                    for kmg = 1:size(keptmuncgroups,1)
                        keptmuncgroupspercell{kmg,1} = munc(munc(:,col.groups)==keptmuncgroups(kmg),:);
                        scatter3(munc(munc(:,col.groups)==keptmuncgroups(kmg),col.x), munc(munc(:,col.groups)==keptmuncgroups(kmg),col.y), munc(munc(:,col.groups)==keptmuncgroups(kmg),col.z)+10, '.g');
                    end
                    
                    
                    save([datasetstrings{dataset}{xi},'_munc.mat'],'keptmuncgroupspercell');
                    save([datasetstrings{dataset}{xi},'_psd.mat'],'keptpsdgroupspercell');
                    saveas(gcf,[datasetstrings{dataset}{xi},'_ExampleImage.fig']);
                    disp(['Finished Count: ',num2str(fincount)]);
                    fincount = fincount+1;
                    close all
                    cd(main)
                end
              cd(round);  
            end
             cd(week)
         end
%         cd(project);
%     end
    %cd(main)
% end
disp('Done!')
toc
%%
clear
clc
dslabels = {'CaMKII';'PV'};
for ii = 1:size(dslabels,1)
    ds = dslabels{ii};
    %Pull in the associated dataset but then convert to a single matrix
    openfig([ds,'_ExampleImage.fig']);
    load([ds,'_munc.mat']);
    load([ds,'_psd.mat']);
    
    rgnmunc = [];
    for i = 1:size(keptmuncgroupspercell,1)
        rgnmunc = [rgnmunc;keptmuncgroupspercell{i}];
    end
    rgnmunc = [rgnmunc ones(size(rgnmunc,1),1)];
    rgnpsd = [];
    for i = 1:size(keptpsdgroupspercell,1)
        rgnpsd = [rgnpsd;keptpsdgroupspercell{i}];
    end
    rgnpsd = [rgnpsd ones(size(rgnpsd,1),1).*2];
    qc = 0;
    while qc == 0 
    
    h = drawfreehand;                                           % allow drawing freehand ROI
    qc = input('If happy, enter 1. if not, just enter or 0');
    if isempty(qc)
        qc=0;
    end
    end
    psdinds = inROI(h, rgnpsd(:,2), rgnpsd(:,3));                  % select locs inside ROI
    muncinds = inROI(h, rgnmunc(:,2), rgnmunc(:,3));
    if isempty(psdinds)
        close all
        continue
    end
    thispsd = rgnpsd(psdinds,:);
    thismunc = rgnmunc(muncinds,:);
    
    save([ds,'_munc_topick.mat'],'thismunc');
    save([ds,'_psd_topick.mat'],'thispsd');
    close all
end
disp('Done! Off to the next one!');
%%

clear
clc
tic
fincount = 1;
datasetstrings{1} = {'CaMKII';'PV'};
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        roundir = dir();
        roundir = roundir(3:end,:);
        for iii = 1:size(roundir,1)
            week = cd(roundir(iii).name);
            targetfolder = dir();
            targetfolder = targetfolder(3:end,:);
            for iv = 1:size(targetfolder,1)
                folderdir = dir();
                folderdir = folderdir(3:end,:);
                
                if contains(folderdir(iv).name,'Region')
                    round = cd(folderdir(iv).name);
                    ds = i;
                    thisds = datasetstrings{ds};
                    for vi = 1:size(thisds,1)
                        if exist([thisds{vi},'_combopickedlocs.mat'],'file') == 2
                            continue
                        end
                        load([thisds{vi},'_munc_topick.mat']);
                        load([thisds{vi},'_psd_topick.mat']);
                        combolocs = [thismunc;thispsd];
                        if isempty(combolocs)
                            continue
                        end
                        
                        syns = dbscan(combolocs(:,[2:3]),.3,4); 
                        
                        locs = [combolocs syns];
                        save([thisds{vi},'_combopickedlocs.mat'],'locs');
                     
                    end
                    
                    cd(round);
                else
                    continue
                end
            end
            disp(['done with project ', num2str(i), ', week ', num2str(ii), ', round ' num2str(iii)]);      
                 cd(week);
            
        end
        cd(project);
        
    end
    cd(main)
end
disp('Done!')
toc
%%
clear
clc
dirname = strsplit(pwd,'\');
groupname = dirname{end};

load([groupname,'_combopickedlocs.mat'])
thismunc = locs(locs(:,15)==1,:);
thispsd = locs(locs(:,15)==2,:);

muncsyns = unique(thismunc(:,16));
psdsyns = unique(thismunc(:,16));
psdstosort = thispsd(:,[2 3 16]);

psdBySyn = sortByField(thispsd,16);

locsperpsd = cellfun(@(x) size(x,1), psdBySyn, 'un',0);
lowlocsyns = cell2mat(locsperpsd)<20;
psdBySyn = psdBySyn(~lowlocsyns);
tic
[psdaxisratio,~,~] = cellfun(@(x) getPSDaxes(x(:,2:3)), psdBySyn, 'un',0,'ErrorHandler',@errorFunc);
toc

revisitPSDs = psdBySyn(isnan(cell2mat(psdaxisratio)));
%%
toredo = cat(1,psdBySyn{~isnan(cell2mat(psdaxisratio))});
psdBySyn2 = sortByField(toredo,16);
tic
[psdaxisratio2,~,~] = cellfun(@(x) getPSDaxes(x(:,2:3)), psdBySyn2, 'un',0,'ErrorHandler',@errorFunc);
toc

tic
[psdarea] = cellfun(@(x) area(alphaShape(x(:,2),x(:,3),'HoleThreshold',100000000)), psdBySyn2,'un',0,'ErrorHandler',@errorFunc);
toc

tic
munclocsinpsd = cellfun(@(x) sum(inShape(alphaShape(x(:,2),x(:,3),'HoleThreshold',100000000),thismunc(thismunc(:,16)==x(1,16),2), thismunc(thismunc(:,16)==x(1,16),3))),psdBySyn2,'un',0,'ErrorHandler',@errorFunc);
toc

ratiokeepinds = cell2mat(psdaxisratio2)<=4; 
ratiounkeptinds = cell2mat(psdaxisratio2)>4;
areakeepinds = cell2mat(psdarea)>1.5 & cell2mat(psdarea)<=20;
keepinds = sum([ratiokeepinds,areakeepinds],2)==2;


CheckSyns = psdBySyn2(keepinds);

%%
close all
KeepSyns = {};
for i = 1:size(CheckSyns,1)
    checkmunc = thismunc(thismunc(:,16)==CheckSyns{i}(1,16),:);
    figure
    hold on
    axis equal
    scatter3(CheckSyns{i}(:,2), CheckSyns{i}(:,3), CheckSyns{i}(:,13),'.r'); %z is 13 
    scatter3(checkmunc(:,2),checkmunc(:,3),checkmunc(:,13),'.b');
    title([num2str(i),' of ',num2str(size(CheckSyns,1))]);
    answer = input('0 if disagree, 1 if agree, 2 to redraw');
    if answer == 0
        close all
    elseif answer == 1
        KeepSyns{i} = [CheckSyns{i};checkmunc];
        close all
    elseif answer == 2
         h = drawfreehand;
         psdtokeep = inROI(h,CheckSyns{i}(:,2), CheckSyns{i}(:,3));
         munctokeep = inROI(h,checkmunc(:,2),checkmunc(:,3));
         KeepSyns{i} = [CheckSyns{i}(psdtokeep,:);checkmunc(munctokeep,:)];
        close all
    end
end
%%
KeepSyns = KeepSyns';

SynsToSave = KeepSyns(~cellfun('isempty',KeepSyns));

synconcat = cat(1,SynsToSave{:});

Syn_indtform = [unique(synconcat(:,16)) [1:(size(unique(synconcat(:,16)),1))]'];
for synind = 1:size(Syn_indtform,1)
    synconcat(synconcat(:,16)==Syn_indtform(synind,1),16) = Syn_indtform(synind,2);
end


save([groupname,'_FinalSynPick_byCell.mat'],'SynsToSave');
save([groupname,'_FinalSynPic_byLoc.mat'],'synconcat');

%% 
clear
clc
close all
tic
CaMKII3Dcheck = {};
PV3Dcheck = {};
CaMKIIsynnum = 1;
PVsynnum = 1;
filedir = dir();
filedir = filedir(3:end,:);
for i = 1:size(filedir,1)
    main = cd(filedir(i).name);
    weekdir = dir();
    weekdir = weekdir(3:end,:);
    for ii = 1:size(weekdir,1)
        project = cd(weekdir(ii).name);
        
        targetfolder = dir();
        targetfolder = targetfolder(3:end,:);
        for iv = 1:size(targetfolder,1)
            if contains(targetfolder(iv).name,'Cell') && contains(targetfolder(iv).name,'CaMKII')
                check = load('CaMKII_FinalSynPick_byCell.mat');
                synstomove = check.SynsToSave;
                for i = 1:length(synstomove)
                    CaMKII3Dcheck(CaMKIIsynnum) = synstomove(i);
                    CaMKIIsynnum = CaMKIIsynnum+1;
                end
            elseif contains(targetfolder(iv).name,'Cell') && contains(targetfolder(iv).name,'PV')
                check = load('PV_FinalSynPick_byCell.mat');
                synstomove = check.SynsToSave;
                for i = 1:length(synstomove)
                    PV3Dcheck(PVsynnum) = synstomove(i);
                    PVsynnum = PVsynnum+1;
                end
            end
        end
        cd(project);
        
    end
    cd(main)
end
%%
dataset = PV3Dcheck;
judgment = {};
istart = 1;
if exist('sabbatical.mat','file') == 2
    load('sabbatical.mat');
    istart = i;
end

for i = istart:length(dataset)
    thissyn = dataset{i};
    thispsd = thissyn(thissyn(:,15)==2,:);
    thismunc = thissyn(thissyn(:,15)==1,:);
    
    thispsd = thispsd(~isoutlier(thispsd(:,13)),:);
    thismunc = thismunc(~isoutlier(thismunc(:,13)),:);
    z_zero_psd = thispsd(:,13)-mean(thispsd(:,13));
    z_zero_munc = thismunc(:,13)-mean(thismunc(:,13));
    
    pctpsdoob = (sum([z_zero_psd>1;z_zero_psd<-1])/length(z_zero_psd))*100;
    pctmuncoob = (sum([z_zero_munc>1;z_zero_munc<-1])/length(z_zero_munc))*100;
    
    figure('WindowState','maximized');
    t = tiledlayout(2,4);
    nexttile
    scatter3(thispsd(:,2), thispsd(:,3),z_zero_psd,'.r');
    zlim([-5 5]);
    view(2)
    nexttile
    scatter3(thispsd(:,2), thispsd(:,3), z_zero_psd,'.r');
    zlim([-5 5]);
    view(3)
    nexttile
    scatter3(thismunc(:,2), thismunc(:,3), z_zero_munc,'.b');
    zlim([-5 5]);
    view(2)
    nexttile
    scatter3(thismunc(:,2), thismunc(:,3), z_zero_munc,'.b');
    zlim([-5 5]);
    view(3)
    nexttile
    scatter3(thispsd(:,2), thispsd(:,3), z_zero_psd,'.r');
    zlim([-5 5]);
    view(0,0);
    nexttile
    scatter3(thispsd(:,2), thispsd(:,3), z_zero_psd,'.r');
    zlim([-5 5]);
    view(90,0);
    title(num2str(pctpsdoob))
    nexttile
    scatter3(thismunc(:,2), thismunc(:,3), z_zero_munc,'.b');
    zlim([-5 5]);
    view(0,0);
    nexttile
    scatter3(thismunc(:,2), thismunc(:,3), z_zero_munc,'.b');
    zlim([-5 5]);
    view(90,0);
    title(num2str(pctmuncoob))
    title(t,[num2str(i) ' of ' num2str(length(dataset))]);
    
    call = input('Enter judgment - 1 for good Munc, 2 for good PSD-95, 3 for good both, 0 for trash');
    if call ~= 9
        judgment{i,1} = call;
        judgment{i,2} = [thismunc;thispsd];
    elseif call == 9
        save('sabbatical.mat');
        break
    end
    
    close all
end
PV3DjudgedSyns = judgment;
save('PV3DjudgedSyns.mat','PV3DjudgedSyns');

%%
FinalCaMKIISyns = CaMKII3DjudgedSyns;

for i = 1:length(FinalCaMKIISyns)
    if isempty(FinalCaMKIISyns{i,1})
        continue
    end
    thissyn = FinalCaMKIISyns{i,1};
    thispsd = thissyn(thissyn(:,15)==2,:);
    thismunc = thissyn(thissyn(:,15)==1,:);
    
    if~isempty(thispsd)
        
        thispsdcx = mean(thispsd(:,2));
        thispsdcy = mean(thispsd(:,3));
        thispsdcz = mean(thispsd(:,10));
        thispsdstdx = std(thispsd(:,2));
        thispsdstdy = std(thispsd(:,3));
        thispsdstdz = std(thispsd(:,10));
        
        
        logicalpsdx = thispsd(:,2)<=(thispsdcx+(2*thispsdstdx)) & thispsd(:,2)>=(thispsdcx-(2*thispsdstdx));
        logicalpsdy = thispsd(:,3)<=(thispsdcy+(2*thispsdstdy)) & thispsd(:,3)>=(thispsdcy-(2*thispsdstdy));
        logicalpsdz = thispsd(:,10)<=(thispsdcz+(2*thispsdstdz)) & thispsd(:,10)>=(thispsdcz-(2*thispsdstdz));
        keeppsd = sum([double(logicalpsdx),double(logicalpsdy),double(logicalpsdz)],2)==3;
        
        thispsd = thispsd(keeppsd,:);
    end
    
    if~isempty(thismunc)
        thismunccx = mean(thismunc(:,2));
        thismunccy = mean(thismunc(:,3));
        thismunccz = mean(thismunc(:,10));
        thismuncstdx = std(thismunc(:,2));
        thismuncstdy = std(thismunc(:,3));
        thismuncstdz = std(thismunc(:,10));
        
        logicalmuncx = thismunc(:,2)<=(thismunccx+(2*thismuncstdx)) & thismunc(:,2)>=(thismunccx-(2*thismuncstdx));
        logicalmuncy = thismunc(:,3)<=(thismunccy+(2*thismuncstdy)) & thismunc(:,3)>=(thismunccy-(2*thismuncstdy));
        logicalmuncz = thismunc(:,10)<=(thismunccz+(2*thismuncstdz)) & thismunc(:,10)>=(thismunccz-(2*thismuncstdz));
        keepmunc = sum([double(logicalmuncx),double(logicalmuncy),double(logicalmuncz)],2)==3;
        
        thismunc = thismunc(keepmunc,:);
    end
    
    finalsyn = [thismunc;thispsd];
    
    FinalCaMKIISyns{i,1} = finalsyn;
end

FinalPVSyns = PV3DjudgedSyns;

for i = 1:length(FinalPVSyns)
    if isempty(FinalPVSyns{i,1})
        continue
    end
    thissyn = FinalPVSyns{i,1};
    thispsd = thissyn(thissyn(:,15)==2,:);
    thismunc = thissyn(thissyn(:,15)==1,:);
    
    if~isempty(thispsd)
        
        thispsdcx = mean(thispsd(:,2));
        thispsdcy = mean(thispsd(:,3));
        thispsdcz = mean(thispsd(:,10));
        thispsdstdx = std(thispsd(:,2));
        thispsdstdy = std(thispsd(:,3));
        thispsdstdz = std(thispsd(:,10));
        
        
        logicalpsdx = thispsd(:,2)<=(thispsdcx+(2*thispsdstdx)) & thispsd(:,2)>=(thispsdcx-(2*thispsdstdx));
        logicalpsdy = thispsd(:,3)<=(thispsdcy+(2*thispsdstdy)) & thispsd(:,3)>=(thispsdcy-(2*thispsdstdy));
        logicalpsdz = thispsd(:,10)<=(thispsdcz+(2*thispsdstdz)) & thispsd(:,10)>=(thispsdcz-(2*thispsdstdz));
        keeppsd = sum([double(logicalpsdx),double(logicalpsdy),double(logicalpsdz)],2)==3;
        
        thispsd = thispsd(keeppsd,:);
    end
    
    if~isempty(thismunc)
        thismunccx = mean(thismunc(:,2));
        thismunccy = mean(thismunc(:,3));
        thismunccz = mean(thismunc(:,10));
        thismuncstdx = std(thismunc(:,2));
        thismuncstdy = std(thismunc(:,3));
        thismuncstdz = std(thismunc(:,10));
        
        logicalmuncx = thismunc(:,2)<=(thismunccx+(2*thismuncstdx)) & thismunc(:,2)>=(thismunccx-(2*thismuncstdx));
        logicalmuncy = thismunc(:,3)<=(thismunccy+(2*thismuncstdy)) & thismunc(:,3)>=(thismunccy-(2*thismuncstdy));
        logicalmuncz = thismunc(:,10)<=(thismunccz+(2*thismuncstdz)) & thismunc(:,10)>=(thismunccz-(2*thismuncstdz));
        keepmunc = sum([double(logicalmuncx),double(logicalmuncy),double(logicalmuncz)],2)==3;
        
        thismunc = thismunc(keepmunc,:);
    end
    
    finalsyn = [thismunc;thispsd];
    
    FinalPVSyns{i,1} = finalsyn;
end
save('FinalCaMKIISyns.mat','FinalCaMKIISyns');
save('FinalPVSyns.mat','FinalPVSyns');
%% Analysis script

tic
clear
pathtorequiredfunctions = 'Z:\AAAAAA\Function Repository\Aarons updated NC analysis code\';
addpath(genpath(pathtorequiredfunctions))

% AUTO- AND CROSS-CORRELATION PARAMETERS

minvol = 0.15;          % min synaptic volume (in voxels)
psize = 100;            % camera nm per pixels
shft = [50,180];        % [min max] distance between two synaptic clusters for crossC
renderpixel = 5;        % render pixel size (nm) for autoC and crossC
rmax = 50;              % range of autoC and crossC to correlate (in pixel)

% NANOCLUSTER DETECTION/ENRICHMENT PARAMETERS
radius = (5:10:330)/psize; % center of bins to use for NC enrichment see below for conversion to edges
d = diff(radius)/2;
radius = [radius(1)-d(1), radius(1:end-1)+d, radius(end)+d(end)];
radius(2:end) = radius(2:end)+ eps(radius(2:end));
% The preceding 3 lines shift the bin centers (10:20:300) to bin edges for
% use with histcounts rather than hist. Results are still centered on bins
% defined in initial radius variable.
n_randomizations = 50;      % number of times to randomize position of NCs for determining if an NC is aligned

% INITIALIZE OUTPUT ARRAYS
synoutput647 = nan(10000,11 + rmax+1 + 3 + rmax+1 + 5);
synoutput561 = nan(10000,11 + rmax+1 + 3 + rmax+1+1);
nc647output = nan(10000,10 + (length(radius)-1) + 2 + (length(radius)-1));
nc561output = nan(10000,10 + (length(radius)-1) + 2 + (length(radius)-1));
p2pdist647 = [];
p2pdist647cis = [];
p2pdist561 = [];
p2pdist561cis = [];
closestpairs = [];
% INITIALIZE COUNTERS FOR OUTPUT MATRICES
countsynapses647 = 1;
countsynapses561 = 1;
count647NC = 1;
count561NC = 1;

files = dir();
files = files(3:end);
for i = 1:size(files,1)
   if files(i).isdir
       files(i)=[];
   end
end
dsname = {};
for fnum = 1:size(files,1)
    load(files(fnum).name);
    dsname{fnum}=files(fnum).name(1:end-4);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% BEGIN ANALYSIS
ACfailurerec647 = [];
ACfailurerec561 = [];
% loop through all folders
for dsnum = 1:size(files,1)
    thisdataset = {};
    eval(['thisdataset =',dsname{dsnum},';']);
    disp(['Loaded dataset (',dsname{dsnum},') for analysis']);
    nsyn = size(thisdataset,1);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % FOR EACH SYNAPSE
    for thissyn = 1:nsyn
        
        S = 'Processing synapse %.0f\n';
        fprintf(S,thissyn);
        
        if isempty(thisdataset{thissyn,1})
            continue
        end
        
        synapsetag = thisdataset{thissyn,3}; %Synapse tag 1 = Munc only, 2 = PSD-95 only, 3 = both, 0 = skip
        if synapsetag == 0
            continue
        end
        
        loc = thisdataset{thissyn,1}(:,[15 2 3 10]); %Organize loc file as channel, x, y, z
        
        if synapsetag == 1
            loc(loc(:,1)==2,:) = [];
            loc(loc(:,1)==1,1) = 647;
            xyz647 = loc(loc(:,1) == 647,2:4);                      % xyz in channel 647
            synlocnum647 = size(xyz647,1);
            syn_vol_647 = volume(alphaShape(xyz647,150/psize,'HoleThreshold',10000000));
            syndens647 = synlocnum647/syn_vol_647;
        elseif synapsetag == 2
            loc(loc(:,1)==1,:) = [];
            loc(loc(:,1)==2,1) = 561;
            xyz561 = loc(loc(:,1) == 561,2:4);                      % xyz in channel 561
            synlocnum561 = size(xyz561,1);
            syn_vol_561 = volume(alphaShape(xyz561,150/psize,'HoleThreshold',10000000));
            syndens561 = synlocnum561/syn_vol_561;
        elseif synapsetag == 3
            loc(loc(:,1)==1,1) = 647;
            loc(loc(:,1)==2,1) = 561;
            xyz647 = loc(loc(:,1) == 647,2:4);                      % xyz in channel 647
            synlocnum647 = size(xyz647,1);
            syn_vol_647 = volume(alphaShape(xyz647,150/psize,'HoleThreshold',10000000));
            xyz561 = loc(loc(:,1) == 561,2:4);                      % xyz in channel 561
            synlocnum561 = size(xyz561,1);
            syn_vol_561 = volume(alphaShape(xyz561,150/psize,'HoleThreshold',10000000));
            syndens647 = synlocnum647/syn_vol_647;
            syndens561 = synlocnum561/syn_vol_561;
        end
        
        %~~~~~~~~~~~GET SYNAPSE VOLUMES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % if less than 10 locs in either channel then reject
        if size(xyz647,1) < 10 || size(xyz561,1) < 10
            continue;
        end
        
        
        % %~~~~~~~~~~~DO THE AUTOCORRELATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        S = 'Doing autoC on synapse %.0f\n';
        fprintf(S,thissyn);
        warning('off','MATLAB:alphaShape:DupPointsBasicWarnId'); % suppress a warning from alphaShape about duplicate points, it's fine.
        if synapsetag ~= 2
            try
                autoC647 = get_autocorr_3d_ADL(xyz647, renderpixel, rmax, psize, 0);
            catch
                ACfailurerec647 = [ACfailurerec647;dsnum thissyn];
                autoC647 = NaN(1,length([0:renderpixel:(renderpixel*rmax)]));
            end
        end
        if synapsetag ~= 1
            try
                autoC561 = get_autocorr_3d_ADL(xyz561, renderpixel, rmax, psize, 0);
            catch
                ACfailurerec561 = [ACfailurerec561;dsnum thissyn];
                autoC561 = NaN(1,length([0:renderpixel:(renderpixel*rmax)]));
            end
        end
        
        
        %~~~~~~~~~~~DO THE CROSSCORRELATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Note I got rid of the option to load old ShiftV because I
        % want it to recalculate the whole crosscorrelation each time
        if synapsetag == 3
            S = 'Doing crossC on synapse %.0f\n';
            fprintf(S,thissyn);
            CC = get_crosscorr_3d_ADL(xyz647 ,xyz561, renderpixel, psize, rmax, [], shft, 0); % get crossC
            C = mean(CC(6:14));      % bin(n) is centered on ((bin#-1)*renderpixel) and includes (renderpixel/2) nm. ie bin 9 is 37.5-42.5 nm. 6:14 is bins 1:9
            shiftV = CC(1:3);        % shift vector to get the maximal overlap [dx dy dz]
            
            
            warning('on','MATLAB:alphaShape:DupPointsBasicWarnId');     % turn the warnings back on
            
            %~~~~~~~~~~~SHIFT 561 LOCS ONTO 647 LOCS BY SHIFTV~~~~~~~~~~~~~~~~~~~~~~~~
            xyz561(:,1) = xyz561(:,1) + shiftV(1);
            xyz561(:,2) = xyz561(:,2) + shiftV(2);
            xyz561(:,3) = xyz561(:,3) + shiftV(3);
        else
            C = NaN;
            CC = NaN(1,length([0:renderpixel:(renderpixel*rmax)])+5);
            shiftV = NaN(1,3);
            
        end
        
        %~~~~~~~~~~~DETECT NANOCLUSTERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        S = 'Finding NCs in synapse %.0f\n';
        fprintf(S,thissyn);
        if synapsetag ~= 2
            NC_647 = mmdbasedDBSCAN3D(xyz647,1.6);
            n_NC647 = max(NC_647(:,1));
            totalncvol647 = 0;
            for ncnum = 1:n_NC647
                thisncvol = volume(alphaShape(xyz647(NC_647==ncnum,:),150/psize));
                totalncvol647= totalncvol647+thisncvol;
            end
        end
        if synapsetag ~= 1
            NC_561 = mmdbasedDBSCAN3D(xyz561,2.7);
            n_NC561 = max(NC_561(:,1));
            totalncvol561 = 0;
            for ncnum = 1:n_NC561
                thisncvol = volume(alphaShape(xyz561(NC_561==ncnum,:),150/psize));
                totalncvol561= totalncvol561+thisncvol;
            end
        end
        
        
        %~~~~~~~~~~~SAVE SYNAPTIC DATA TO OUTPUT VARIABLES FOR EACH CHANNEL~~~~~~~
        if synapsetag ~= 2
            synoutput647(countsynapses647,:) = [dsnum synapsetag thissyn 647 ...
                n_NC647 0 0 autoC647 0 0 C CC(:,6:end) 0 shiftV 0 synlocnum647 syn_vol_647 syndens647 totalncvol647];
            countsynapses647 = countsynapses647+1;
        end
        if synapsetag ~= 1
            synoutput561(countsynapses561,:) = [dsnum synapsetag thissyn 561 ...
                n_NC561 0 0 autoC561 0 0 C CC(:,6:end) 0 synlocnum561 syn_vol_561 syndens561 totalncvol561];
            countsynapses561 = countsynapses561+1;
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~BEGIN NC DISTRIBUTION ANALYSIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~GET RANDOM SYNAPTIC CLUSTERS FOR NORMALIZATION~~~~~~~~~~~~~~~~
        factor = 3;     %density factor to generating the rand cluster
        binnum = length(radius)-1;
        if synapsetag~=2
            rand_cluster_647 = get_cluster_randomized_ADL(xyz647,factor,psize,0);
        end
        if synapsetag~=1
            rand_cluster_561 = get_cluster_randomized_ADL(xyz561,factor,psize,0);
        end
        S = 'Processing nanoclusters in synapse %.0f\n';
        fprintf(S,thissyn);
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~PROCESS 647 NANOCLUSTERS FOR VOLUME AND ENRICHMENT~~~~~~~~~~~~
        if synapsetag ~= 2
            for this647NC = 1:max(NC_647)
                
                NClocnum647 = find(NC_647 == this647NC);
                NClocnum647 = size(NClocnum647,1);
                
                % reject NC if fewer than 4 points
                if NClocnum647 < 5
                    continue;
                end
                
                %~~~~~~~~~~~~~~~GET TESSELATION VOLUME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                NCvol647 = volume(alphaShape(xyz647(NC_647==this647NC,:),150/psize));
                %NCvol647tess = sum(volume_tes_3d(xyz647,find(NC_647==this647NC)));
                
                %~~~~~~~~~~~~~~~CALCULATE DENSITY DISTRIBUTIONS RELATIVE TO NC CENTER~~~~~
                center = findDBSCANNCpeak(xyz647(NC_647==this647NC,:),psize,'Tra',1.5,'is3D',1);
                dis647_647NCctr = pdist2(xyz647,center);
                dis647rand_647NCctr = pdist2(rand_cluster_647,center);
                randhistcount = histcounts(dis647rand_647NCctr, radius);
                attemptflag = 0;
                for thisbin = 1:binnum
                    rng('shuffle')
                    loccheck = 0;
                    attemptnum = 1;
                    while loccheck == 0
                        randbin = randi(binnum);
                        loccheck = randhistcount(randbin)>0;
                        attemptnum = attemptnum + 1;
                        if attemptnum > 50
                            attemptflag = 1;
                            break
                        end
                    end
                    randhistcount(randbin)=randhistcount(randbin)-1;
                end
                if attemptflag == 0
                    randhistcount = randhistcount+1;
                elseif attemptflag == 1
                    randhistcount = zeros(1,33);
                end
                dist647_647NC = histcounts(dis647_647NCctr, radius) ./ randhistcount .* factor;
                if synapsetag == 3
                    dis561_647NCctr = pdist2(xyz561,center);
                    dis561rand_647NCctr = pdist2(rand_cluster_561,center);
                    randhistcountcross = histcounts(dis561rand_647NCctr, radius);
                    attemptflag = 0;
                    for thisbin = 1:binnum
                        rng('shuffle')
                        loccheck = 0;
                        attemptnum = 1;
                        while loccheck == 0
                            randbin = randi(binnum);
                            loccheck = randhistcountcross(randbin)>0;
                            attemptnum = attemptnum + 1;
                            if attemptnum > 50
                                attemptflag = 1;
                                break
                            end
                        end
                        randhistcountcross(randbin)=randhistcountcross(randbin)-1;
                    end
                    if attemptflag == 0
                        randhistcountcross = randhistcountcross+1;
                    elseif attemptflag == 1
                        randhistcountcross = zeros(1,33);
                    end
                    dist561_647NC = histcounts(dis561_647NCctr, radius) ./ randhistcountcross .* factor;
                    
                    %~~~~~~~~~~~~~~~DETERMINE WHETHER GIVEN NANOCLUSTER IS ALIGNED~~~~~~~~~~~~
                    
                    EIrand561dist_647NC = nan(n_randomizations/5,1);
                    for randrunnum = 1:n_randomizations/5
                        
                        this_rand_cluster_561 = get_cluster_randomized_ADL(xyz561,1,psize,0);
                        thisdisp561rand_647NCctr = pdist2(this_rand_cluster_561,center);
                        thisdist561_647NC = histcounts(thisdisp561rand_647NCctr,radius) ./ randhistcountcross .* factor;
                        thisdist561_647NC(isnan(thisdist561_647NC))=0;
                        EIrand561dist_647NC(randrunnum,1) = nanmean(thisdist561_647NC(2:6));
                        
                    end
                    alignthresholds = [(nanmean(EIrand561dist_647NC) + (1.96 * nanstd(EIrand561dist_647NC)))...
                        (nanmean(EIrand561dist_647NC) - (1.96 * nanstd(EIrand561dist_647NC)))];
                    testCE = nanmean(dist561_647NC(2:6));
                    
                    if testCE > alignthresholds(1)
                        isaligned647 = 1;
                    elseif testCE < alignthresholds(2)
                        isaligned647 = -1;
                    elseif testCE <= alignthresholds(1) && testCE >= alignthresholds(2)
                        isaligned647 = 0;
                    end
                elseif synapsetag ~= 3
                    isaligned647 = NaN;
                    dist561_647NC = NaN(1,33);
                end
                
                
                %~~~~~~~~~~~~~~~SAVE THE DATA TO OUTPUT VARIABLE AND TICK THE NC COUNTER~~
                nc647output(count647NC,:) = [dsnum synapsetag thissyn 647 ...
                    this647NC NCvol647 NClocnum647 isaligned647 0 0 dist647_647NC 0 0 dist561_647NC];
                count647NC = count647NC + 1;
                
                
            end %647NC
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~PROCESS 561 NANOCLUSTERS FOR VOLUME AND ENRICHMENT~~~~~~~~~~~~
        if synapsetag ~=1
            for this561NC = 1:max(NC_561)
                
                NClocnum561 = find(NC_561 == this561NC);
                NClocnum561 = size(NClocnum561,1);
                % reject NC if fewer than 4 points
                if NClocnum561 < 5
                    continue;
                end
                
                %~~~~~~~~~~~~~~~GET TESSELATION VOLUME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                NCvol561 = volume(alphaShape(xyz561(NC_561==this561NC,:),150/psize));
               % NCvol561tess = sum(volume_tes_3d(xyz561,find(NC_561==this561NC)));
                %~~~~~~~~~~~~~~~CALCULATE DENSITY DISTRIBUTIONS RELATIVE TO NC CENTER~~~~~
                center = findDBSCANNCpeak(xyz561(NC_561==this561NC,:),psize,'Tra',1.5,'is3D',1);
                dis561_561NCctr = pdist2(xyz561,center);
                dis561rand_561NCctr = pdist2(rand_cluster_561,center);
                randhistcount = histcounts(dis561rand_561NCctr, radius);
                attemptflag = 0;
                for thisbin = 1:binnum
                    rng('shuffle')
                    loccheck = 0;
                    attemptnum = 1;
                    while loccheck == 0
                        randbin = randi(binnum);
                        loccheck = randhistcount(randbin)>0;
                        attemptnum = attemptnum + 1;
                        if attemptnum > 50
                            attemptflag = 1;
                            break
                        end
                    end
                    randhistcount(randbin)=randhistcount(randbin)-1;
                end
                if attemptflag == 0
                    randhistcount = randhistcount+1;
                elseif attemptflag == 1
                    randhistcount = zeros(1,33);
                end
                dist561_561NC = histcounts(dis561_561NCctr, radius) ./ randhistcount .* factor;
                if synapsetag == 3
                    dis647_561NCctr = pdist2(xyz647,center);
                    dis647rand_561NCctr = pdist2(rand_cluster_647,center);
                    randhistcountcross = histcounts(dis647rand_561NCctr, radius);
                    attemptflag = 0;
                    for thisbin = 1:binnum
                        rng('shuffle')
                        loccheck = 0;
                        attemptnum = 1;
                        while loccheck == 0
                            randbin = randi(binnum);
                            loccheck = randhistcountcross(randbin)>0;
                            attemptnum = attemptnum + 1;
                            if attemptnum > 50
                                attemptflag = 1;
                                break
                            end
                        end
                        randhistcountcross(randbin)=randhistcountcross(randbin)-1;
                    end
                    if attemptflag == 0
                        randhistcountcross = randhistcountcross+1;
                    elseif attemptflag == 1
                        randhistcountcross = zeros(1,33);
                    end
                    dist647_561NC = histcounts(dis647_561NCctr, radius) ./ randhistcountcross .* factor;
                    
                    %~~~~~~~~~~~~~~~DETERMINE WHETHER GIVEN NANOCLUSTER IS ALIGNED~~~~~~~~~~~~
                    
                    EIrand647dist_561NC = nan(n_randomizations/5,1);
                    for randrunnum = 1:n_randomizations/5
                        
                        this_rand_cluster_647 = get_cluster_randomized_ADL(xyz647,1,psize,0);
                        thisdisp647rand_561NCctr = pdist2(this_rand_cluster_647,center);
                        thisdist647_561NC = histcounts(thisdisp647rand_561NCctr,radius) ./ randhistcountcross .* factor;
                        thisdist647_561NC(isnan(thisdist647_561NC))=0;
                        EIrand647dist_561NC(randrunnum,1) = nanmean(thisdist647_561NC(2:6));
                        
                    end
                    alignthresholds = [(nanmean(EIrand647dist_561NC) + (1.96 * nanstd(EIrand647dist_561NC)))...
                        (nanmean(EIrand647dist_561NC) - (1.96 * nanstd(EIrand647dist_561NC)))];
                    testCE = nanmean(dist647_561NC(2:6));
                    
                    if testCE > alignthresholds(1)
                        isaligned561 = 1;
                    elseif testCE < alignthresholds(2)
                        isaligned561 = -1;
                    elseif testCE <= alignthresholds(1) && testCE >= alignthresholds(2)
                        isaligned561 = 0;
                    end
                elseif synapsetag ~= 3
                    isaligned561 = NaN;
                    dist647_561NC = NaN(1,33);
                end
                
                %~~~~~~~~~~~~~~~SAVE THE DATA TO OUTPUT VARIABLE AND TICK THE NC COUNTER~~
                nc561output(count561NC,:) = [dsnum synapsetag thissyn 561 ...
                    this561NC NCvol561 NClocnum561 isaligned561 0 0 dist561_561NC 0 0 dist647_561NC];
                count561NC = count561NC + 1;
                
            end %561 NC
        end
        
        if synapsetag ~=2
            if n_NC647 < 1
                continue
            end
            ncstorage647 = {};
            storageind647 = 1;
            for ncnum647 = 1:n_NC647
                thisnc647 = xyz647(NC_647==ncnum647,:);
                if size(thisnc647,1)<5
                    continue
                end
                ncstorage647{storageind647,1} = thisnc647;
                storageind647 = storageind647+1;
            end
            nccenters647 = cell2mat(cellfun(@(x) findDBSCANNCpeak(x,psize,'Tra',1.5,'is3D',1),ncstorage647,'UniformOutput',0));
            
            intraNCdistances_cis_munc = pdist2(nccenters647,nccenters647);
            
            
            if n_NC647 <2
                for i = 1:size(intraNCdistances_cis_munc,1)
                    p2pdist647cis = [p2pdist647cis;dsnum synapsetag thissyn 647 ...
                        i NaN];
                end
            elseif n_NC647 >=2
                for i = 1:size(intraNCdistances_cis_munc,1)
                    p2pdist647cis = [p2pdist647cis;dsnum synapsetag thissyn 647 ...
                        i min(intraNCdistances_cis_munc(i,intraNCdistances_cis_munc(i,:)>0,:))];
                end
            end
            
        end
        
        if synapsetag ~=1
            if n_NC561 < 1
                continue
            end
            ncstorage561 = {};
            storageind561 = 1;
            for ncnum561 = 1:n_NC561
                thisnc561 = xyz561(NC_561==ncnum561,:);
                if size(thisnc561,1)<5
                    continue
                end
                ncstorage561{storageind561,1} = thisnc561;
                storageind561 = storageind561+1;
            end
            nccenters561 = cell2mat(cellfun(@(x) findDBSCANNCpeak(x,psize,'Tra',1.5,'is3D',1),ncstorage561,'UniformOutput',0));
            
            intraNCdistances_cis_psd = pdist2(nccenters561,nccenters561);
            
            
            if n_NC561 <2
                for i = 1:size(intraNCdistances_cis_psd,1)
                    p2pdist561cis = [p2pdist561cis;dsnum synapsetag thissyn 561 ...
                        i NaN];
                end
            elseif n_NC647 >=2
                for i = 1:size(intraNCdistances_cis_psd,1)
                    p2pdist561cis = [p2pdist561cis;dsnum synapsetag thissyn 561 ...
                        i min(intraNCdistances_cis_psd(i,intraNCdistances_cis_psd(i,:)>0,:))];
                end
            end
            
        end
        
        if synapsetag == 3
            if n_NC647 < 1 || n_NC561 < 1
                continue
            end
            ncstorage647 = {};
            ncstorage561 = {};
            storageind647 = 1;
            storageind561 = 1;
            for ncnum647 = 1:n_NC647
                thisnc647 = xyz647(NC_647==ncnum647,:);
                if size(thisnc647,1)<5
                    continue
                end
                ncstorage647{storageind647,1} = thisnc647;
                storageind647 = storageind647+1;
            end
            for ncnum561 = 1:n_NC561
                thisnc561 = xyz561(NC_561==ncnum561,:);
                if size(thisnc561,1)<5
                    continue
                end
                ncstorage561{storageind561,1} = thisnc561;
                storageind561 = storageind561+1;
            end
            
            nccenters647 = cell2mat(cellfun(@(x) findDBSCANNCpeak(x,psize,'Tra',1.5,'is3D',1),ncstorage647,'UniformOutput',0));
            nccenters561 = cell2mat(cellfun(@(x) findDBSCANNCpeak(x,psize,'Tra',1.5,'is3D',1),ncstorage561,'UniformOutput',0));
            
            intraNCdistances = pdist2(nccenters647,nccenters561);
            
            
            for i = 1:size(intraNCdistances,1)
                p2pdist647 = [p2pdist647;dsnum synapsetag thissyn 647 ...
                    i min(intraNCdistances(i,:))];
            end
            
            for i = 1:size(intraNCdistances,2)
                p2pdist561 = [p2pdist561;dsnum synapsetag thissyn 561 ...
                    i min(intraNCdistances(:,i))];
            end
            
            % % %            for ncnum = 1:size(intraNCdistances,1)
            % % %                rowmin = min(intraNCdistances(ncnum,:));
            % % %                [~,intraNDcol] = ind2sub(size(intraNCdistances),find(intraNCdistances==rowmin));
            % % %                if min(intraNCdistances(:,intraNDcol))==rowmin
            % % %                   closestpairs = [closestpairs; dsnum synapsetag thissyn 647 ncnum 0 0 ...
            % % %                       crossenrichmentstorage647(ncnum,:) 999 999 dsnum synapsetag thissyn 561 ...
            % % %                       intraNDcol 0 0 crossenrichmentstorage561(intraNDcol,:)];
            % % %                end
            % % %            end
        end
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        S = 'Analysis of synapse %.0f of %.0f in dataset %.0f done!\n';
        fprintf(S,thissyn,nsyn,dsnum);
        
        
    end %synapse
    
end %folder
disp('all done!');


% TRIM EXCESS ROWS FROM OUTPUT VARIABLES
synoutput647(isnan(synoutput647(:,1)),:) = [];
synoutput561(isnan(synoutput561(:,1)),:) = [];
nc647output(isnan(nc647output(:,1)),:) = [];
nc561output(isnan(nc561output(:,1)),:) = [];

% SAVE THE DATA
save('Synaptic_data_Munc.txt', 'synoutput647', '-ascii', '-tabs');
save('Synaptic_data_PSD.txt', 'synoutput561', '-ascii', '-tabs');
save('NC_data_Munc.txt', 'nc647output', '-ascii', '-tabs');
save('NC_data_PSD.txt', 'nc561output', '-ascii', '-tabs');

save('p2pdistances_Munctrans.txt','p2pdist647','-ascii','-tabs');
save('p2pdistances_PSDtrans.txt','p2pdist561','-ascii','-tabs');
save('p2pdistances_MuncCis.txt','p2pdist647cis','-ascii','-tabs');
save('p2pdistances_PSDCis.txt','p2pdist561cis','-ascii','-tabs');
toc
% SYNAPTIC_DATA RESULTS STRUCTURE
% 1: folder
% 2: file
% 3: synapse number
% 4: channel
% 5: synapse alpha volume of this channel
% 6: number of nanoclusters in this channel
% 7-8: zeros/break
% 9-next zeros: autocorrelation for this channel
% after zeros-end: crosscorrelation for this channel (will be
% the same for both channels for a given synapse)

% NC_DATA RESULTS STRUCTURE
% 1: folder
% 2: file
% 3: synapse number
% 4: channel
% 5: nanocluster number (ie identifier, not the number of nc/syn)
% 6: tesselation volume of this NC
% 7: is this nanocluster aligned across the synapse (1 for yes 0 for no)
% 8-9: zeros/break
% 10-next zeros: autoenrichment
% after zeros-end: cross-enrichment

%% Data Management, outlier removal, smoothing, etc.
clear
clc
%This first step is contingent on manually removing synapses based on PRISM
%ROUT 0.1% outlier removal on synapse volume within group and protein 
psddata = xlsread('Synaptic_data_PSD.xlsx');
muncdata = xlsread('Synaptic_data_Munc.xlsx');
%First parse out the cross syn data
psdcross = psddata(psddata(:,2)==3,:);
munccross = muncdata(muncdata(:,2)==3,:);
[finalcrosssyns,indpsd,indmunc] = intersect(psdcross(:,1:3),munccross(:,1:3),'rows','stable');
crosspsdfinal = psdcross(indpsd,:);
crossmuncfinal = munccross(indmunc,:);

xlswrite('Synaptic_data_Munc_Cross.xlsx',crossmuncfinal);
xlswrite('Synaptic_data_PSD_cross.xlsx',crosspsdfinal);
%Next wittle down the NCs by removing the ones in unanalyzed synapses
ncmuncdata = xlsread('NC_data_Munc.xlsx');
[muncncpass] = intersect(ncmuncdata(:,1:3),muncdata(:,1:3),'rows','stable');
indmuncnc = ismember(ncmuncdata(:,1:3),muncncpass,'rows');
ncmuncsynoutlierpass = ncmuncdata(indmuncnc,:);
xlswrite('NC_data_Munc_synoutpass.xlsx',ncmuncsynoutlierpass);

ncpsddata = xlsread('NC_data_PSD.xlsx');
[psdncpass] = intersect(ncpsddata(:,1:3),psddata(:,1:3),'rows','stable');
indpsdnc = ismember(ncpsddata(:,1:3),psdncpass,'rows');
ncpsdsynoutlierpass = ncpsddata(indpsdnc,:);
xlswrite('NC_data_PSD_synoutpass.xlsx',ncpsdsynoutlierpass);

cisp2pmunc = xlsread('p2pdistances_MuncCis.xlsx');
indcisp2pmunc = ismember(cisp2pmunc(:,1:3),ncmuncsynoutlierpass(:,1:3),'rows');
cisp2pmuncsynoutlierpass = cisp2pmunc(indcisp2pmunc,:);
xlswrite('p2pdistances_MuncCis_synoutpass.xlsx',cisp2pmuncsynoutlierpass);

cisp2ppsd = xlsread('p2pdistances_PSDCis.xlsx');
indcisp2ppsd = ismember(cisp2ppsd(:,1:3),ncpsdsynoutlierpass(:,1:3),'rows');
cisp2ppsdsynoutlierpass = cisp2ppsd(indcisp2ppsd,:);
xlswrite('p2pdistances_PSDCis_synoutpass.xlsx',cisp2ppsdsynoutlierpass);

%next do the same as above but using just the cross synapses for the cross
%p2ps
transp2pmunc = xlsread('p2pdistances_Munctrans.xlsx');
indtransp2pmunc = ismember(transp2pmunc(:,1:3),crossmuncfinal(:,1:3),'rows');
transp2pmuncsynoutlierpass = transp2pmunc(indtransp2pmunc,:);
xlswrite('p2pdistances_Munctrans_synoutpass.xlsx',transp2pmuncsynoutlierpass);

transp2ppsd = xlsread('p2pdistances_PSDtrans.xlsx');
indtransp2ppsd = ismember(transp2ppsd(:,1:3),crosspsdfinal(:,1:3),'rows');
transp2ppsdsynoutlierpass = transp2ppsd(indtransp2ppsd,:);
xlswrite('p2pdistances_PSDtrans_synoutpass.xlsx',transp2ppsdsynoutlierpass);

%Next is manual outlier removal in PRISM with ROUT 0.1% for NC vol AFTER
%deciding whether to use tess vol or alpha shape vol (Note: Used AlphaShape
%vol)

%Now we need to filter the p2p NCs to the ones that passed volume outlier
%removal.
finalmuncnc = xlsread('NC_data_Munc_FINAL.xlsx');
cisp2pmunc = xlsread('p2pdistances_MuncCis_synoutpass.xlsx');
transp2pmunc = xlsread('p2pdistances_Munctrans_synoutpass.xlsx');
indcisp2pmunc = ismember(cisp2pmunc(:,1:3),finalmuncnc(:,1:3),'rows');
indtransp2pmunc = ismember(transp2pmunc(:,1:3),finalmuncnc(:,1:3),'rows');
cisp2pmuncfinal = cisp2pmunc(indcisp2pmunc,:);
transp2pmuncfinal = transp2pmunc(indtransp2pmunc,:);
xlswrite('p2pdistances_MuncCis_FINAL.xlsx',cisp2pmuncfinal);
xlswrite('p2pdistances_Munctrans_FINAL.xlsx',transp2pmuncfinal);

finalpsdnc = xlsread('NC_data_PSD_FINAL.xlsx');
cisp2ppsd = xlsread('p2pdistances_PSDCis_synoutpass.xlsx');
transp2ppsd = xlsread('p2pdistances_PSDtrans_synoutpass.xlsx');
indcisp2ppsd = ismember(cisp2ppsd(:,1:3),finalpsdnc(:,1:3),'rows');
indtransp2ppsd = ismember(transp2ppsd(:,1:3),finalpsdnc(:,1:3),'rows');
cisp2ppsdfinal = cisp2ppsd(indcisp2ppsd,:);
transp2ppsdfinal = transp2ppsd(indtransp2ppsd,:);
xlswrite('p2pdistances_PSDCis_FINAL.xlsx',cisp2ppsdfinal);
xlswrite('p2pdistances_PSDtrans_FINAL.xlsx',transp2ppsdfinal);
%% This section is to help smooth the auto and cross enrichments
%First run the ROUT 0.1% for each column and then replace outliers with the
%highest non-outlier value. 
clear
clc
temp = xlsread('enrichmentoutlierremovaltemp.xlsx');
data = temp(2:end,:);
maxvals = temp(1,:);
newdata = [];
for i = 1:size(data,2)
   thiscol = data(:,i);
   thiscol(thiscol>maxvals(i)) = maxvals(i);
   newdata = [newdata thiscol];
end

%% JI
clear
clc

AllMuncData = xlsread("FINALMuncNCsForEIAnalysis_withcenter.xlsx");

SW30data = AllMuncData(:,[1:3,5,8:10,27:38]);
DataByCondition = sortByField(SW30data,1);

DataByConditionAndSynapse = cellfun(@(x) sortByField(x,3),DataByCondition,'uni',0);

DataByConditionAndSynapse = cellfun(@(x) x(~cellfun('isempty',x)),DataByConditionAndSynapse,'uni',0);
JI_CaMKII = cellfun(@(x) JaccardIndex(x),DataByConditionAndSynapse{1},'uni',0);
JI_PV = cellfun(@(x) JaccardIndex(x),DataByConditionAndSynapse{2},'uni',0);
ncpersyn_camkii = cell2mat(cellfun(@(x) max(x(:,2)),JI_CaMKII,'uni',0));
ncpersyn_pv = cell2mat(cellfun(@(x) max(x(:,2)),JI_PV,'uni',0));

%First generate the total synaptic data for all AZs with at least 2 NCs of
%Munc13-1

JI_CaMKII_TotalSyns = JI_CaMKII(ncpersyn_camkii>1);
JI_PV_TotalSyns = JI_PV(ncpersyn_pv>1);

JI_CaMKII_TotalSyns_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_CaMKII_TotalSyns,'uni',0));
JI_PV_TotalSyns_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_PV_TotalSyns,'uni',0));

JI_CaMKII_TotalSynsWithNoZeroes = JI_CaMKII_TotalSyns_meanJIpersyn(JI_CaMKII_TotalSyns_meanJIpersyn>0);
JI_PV_TotalSynsWithNoZeroes = JI_PV_TotalSyns_meanJIpersyn(JI_PV_TotalSyns_meanJIpersyn>0);

JI_CaMKII_AllNCsPooled = [];
JI_PV_AllNCsPooled = [];

for i = 1:size(JI_CaMKII_TotalSyns,1)
    JI_CaMKII_AllNCsPooled = [JI_CaMKII_AllNCsPooled; JI_CaMKII_TotalSyns{i}];
end
JI_CaMKII_AllNCsPooled = JI_CaMKII_AllNCsPooled(:,[3 6]);


for i = 1:size(JI_PV_TotalSyns,1)
    JI_PV_AllNCsPooled = [JI_PV_AllNCsPooled; JI_PV_TotalSyns{i}];
end
JI_PV_AllNCsPooled = JI_PV_AllNCsPooled(:,[3 6]);


JI_CaMKII_AllNCsPooled_zeroes = JI_CaMKII_AllNCsPooled(JI_CaMKII_AllNCsPooled(:,1)==0,:);
JI_PV_AllNCsPooled_zeroes = JI_PV_AllNCsPooled(JI_PV_AllNCsPooled(:,1)==0,:);

JI_CaMKII_AllNCsPooled_nonzeroes = JI_CaMKII_AllNCsPooled(JI_CaMKII_AllNCsPooled(:,1)~=0,:);
JI_PV_AllNCsPooled_nonzeroes = JI_PV_AllNCsPooled(JI_PV_AllNCsPooled(:,1)~=0,:);

CaMKII_JI_quartiles = quantile(JI_CaMKII_AllNCsPooled_nonzeroes,[.25 .5 .75]);
PV_JI_quartiles = quantile(JI_PV_AllNCsPooled_nonzeroes,[.25 .5 .75]);

JI_CaMKII_AllNCsPooled_Q1 = JI_CaMKII_AllNCsPooled_nonzeroes(JI_CaMKII_AllNCsPooled_nonzeroes(:,1)<=CaMKII_JI_quartiles(1,1),:);
JI_CaMKII_AllNCsPooled_Q2 = JI_CaMKII_AllNCsPooled_nonzeroes(JI_CaMKII_AllNCsPooled_nonzeroes(:,1)>CaMKII_JI_quartiles(1,1) & JI_CaMKII_AllNCsPooled_nonzeroes(:,1)<=CaMKII_JI_quartiles(2,1),:);
JI_CaMKII_AllNCsPooled_Q3 = JI_CaMKII_AllNCsPooled_nonzeroes(JI_CaMKII_AllNCsPooled_nonzeroes(:,1)>CaMKII_JI_quartiles(2,1) & JI_CaMKII_AllNCsPooled_nonzeroes(:,1)<=CaMKII_JI_quartiles(3,1),:);
JI_CaMKII_AllNCsPooled_Q4 = JI_CaMKII_AllNCsPooled_nonzeroes(JI_CaMKII_AllNCsPooled_nonzeroes(:,1)>CaMKII_JI_quartiles(3,1),:);

JI_PV_AllNCsPooled_Q1 = JI_PV_AllNCsPooled_nonzeroes(JI_PV_AllNCsPooled_nonzeroes(:,1)<=PV_JI_quartiles(1,1),:);
JI_PV_AllNCsPooled_Q2 = JI_PV_AllNCsPooled_nonzeroes(JI_PV_AllNCsPooled_nonzeroes(:,1)>PV_JI_quartiles(1,1) & JI_PV_AllNCsPooled_nonzeroes(:,1)<=PV_JI_quartiles(2,1),:);
JI_PV_AllNCsPooled_Q3 = JI_PV_AllNCsPooled_nonzeroes(JI_PV_AllNCsPooled_nonzeroes(:,1)>PV_JI_quartiles(2,1) & JI_PV_AllNCsPooled_nonzeroes(:,1)<=PV_JI_quartiles(3,1),:);
JI_PV_AllNCsPooled_Q4 = JI_PV_AllNCsPooled_nonzeroes(JI_PV_AllNCsPooled_nonzeroes(:,1)>PV_JI_quartiles(3,1),:);

CaMKII_MeanByQDdata = [mean(JI_CaMKII_AllNCsPooled_zeroes(:,2)) std(JI_CaMKII_AllNCsPooled_zeroes(:,2)) size(JI_CaMKII_AllNCsPooled_zeroes(:,2),1);
    mean(JI_CaMKII_AllNCsPooled_Q1(:,2)) std(JI_CaMKII_AllNCsPooled_Q1(:,2)) size(JI_CaMKII_AllNCsPooled_Q1(:,2),1); 
    mean(JI_CaMKII_AllNCsPooled_Q2(:,2)) std(JI_CaMKII_AllNCsPooled_Q2(:,2)) size(JI_CaMKII_AllNCsPooled_Q2(:,2),1); 
    mean(JI_CaMKII_AllNCsPooled_Q3(:,2)) std(JI_CaMKII_AllNCsPooled_Q3(:,2)) size(JI_CaMKII_AllNCsPooled_Q3(:,2),1); 
    mean(JI_CaMKII_AllNCsPooled_Q4(:,2)) std(JI_CaMKII_AllNCsPooled_Q4(:,2)) size(JI_CaMKII_AllNCsPooled_Q4(:,2),1);];
CaMKII_MedianByQData = [median(JI_CaMKII_AllNCsPooled_zeroes(:,2)) median(JI_CaMKII_AllNCsPooled_Q1(:,2)) median(JI_CaMKII_AllNCsPooled_Q2(:,2)) median(JI_CaMKII_AllNCsPooled_Q3(:,2)) median(JI_CaMKII_AllNCsPooled_Q4(:,2))];

PV_MeanByQDdata = [mean(JI_PV_AllNCsPooled_zeroes(:,2)) std(JI_PV_AllNCsPooled_zeroes(:,2)) size(JI_PV_AllNCsPooled_zeroes(:,2),1);
    mean(JI_PV_AllNCsPooled_Q1(:,2)) std(JI_PV_AllNCsPooled_Q1(:,2)) size(JI_PV_AllNCsPooled_Q1(:,2),1); 
    mean(JI_PV_AllNCsPooled_Q2(:,2)) std(JI_PV_AllNCsPooled_Q2(:,2)) size(JI_PV_AllNCsPooled_Q2(:,2),1); 
    mean(JI_PV_AllNCsPooled_Q3(:,2)) std(JI_PV_AllNCsPooled_Q3(:,2)) size(JI_PV_AllNCsPooled_Q3(:,2),1); 
    mean(JI_PV_AllNCsPooled_Q4(:,2)) std(JI_PV_AllNCsPooled_Q4(:,2)) size(JI_PV_AllNCsPooled_Q4(:,2),1);];
PV_MedianByQData = [median(JI_PV_AllNCsPooled_zeroes(:,2)) median(JI_PV_AllNCsPooled_Q1(:,2)) median(JI_PV_AllNCsPooled_Q2(:,2)) median(JI_PV_AllNCsPooled_Q3(:,2)) median(JI_PV_AllNCsPooled_Q4(:,2))];

JI_CaMKII_2NCspersyn = JI_CaMKII(ncpersyn_camkii==2);
JI_CaMKII_3NCspersyn = JI_CaMKII(ncpersyn_camkii==3);
JI_CaMKII_4NCspersyn = JI_CaMKII(ncpersyn_camkii==4);
JI_CaMKII_5NCspersyn = JI_CaMKII(ncpersyn_camkii==5);
JI_CaMKII_6NCspersyn = JI_CaMKII(ncpersyn_camkii==6);
JI_CaMKII_7NCspersyn = JI_CaMKII(ncpersyn_camkii==7);

JI_PV_2NCspersyn = JI_PV(ncpersyn_pv==2);
JI_PV_3NCspersyn = JI_PV(ncpersyn_pv==3);
JI_PV_4NCspersyn = JI_PV(ncpersyn_pv==4);
JI_PV_5NCspersyn = JI_PV(ncpersyn_pv==5);
JI_PV_6NCspersyn = JI_PV(ncpersyn_pv==6);
JI_PV_7NCspersyn = JI_PV(ncpersyn_pv==7);

JI_CaMKII_2NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_CaMKII_2NCspersyn,'uni',0));
JI_CaMKII_3NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_CaMKII_3NCspersyn,'uni',0));
JI_CaMKII_4NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_CaMKII_4NCspersyn,'uni',0));
JI_CaMKII_5NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_CaMKII_5NCspersyn,'uni',0));
JI_CaMKII_6NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_CaMKII_6NCspersyn,'uni',0));
JI_CaMKII_7NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_CaMKII_7NCspersyn,'uni',0));

JI_PV_2NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_PV_2NCspersyn,'uni',0));
JI_PV_3NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_PV_3NCspersyn,'uni',0));
JI_PV_4NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_PV_4NCspersyn,'uni',0));
JI_PV_5NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_PV_5NCspersyn,'uni',0));
JI_PV_6NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_PV_6NCspersyn,'uni',0));
JI_PV_7NCspersyn_meanJIpersyn = cell2mat(cellfun(@(x) mean(x(:,3)),JI_PV_7NCspersyn,'uni',0));

camkiidata = [mean(JI_CaMKII_2NCspersyn_meanJIpersyn) std(JI_CaMKII_2NCspersyn_meanJIpersyn) size(JI_CaMKII_2NCspersyn_meanJIpersyn,1);
mean(JI_CaMKII_3NCspersyn_meanJIpersyn) std(JI_CaMKII_3NCspersyn_meanJIpersyn) size(JI_CaMKII_3NCspersyn_meanJIpersyn,1);
mean(JI_CaMKII_4NCspersyn_meanJIpersyn) std(JI_CaMKII_4NCspersyn_meanJIpersyn) size(JI_CaMKII_4NCspersyn_meanJIpersyn,1);
mean(JI_CaMKII_5NCspersyn_meanJIpersyn) std(JI_CaMKII_5NCspersyn_meanJIpersyn) size(JI_CaMKII_5NCspersyn_meanJIpersyn,1);
mean(JI_CaMKII_6NCspersyn_meanJIpersyn) std(JI_CaMKII_6NCspersyn_meanJIpersyn) size(JI_CaMKII_6NCspersyn_meanJIpersyn,1);
mean(JI_CaMKII_7NCspersyn_meanJIpersyn) std(JI_CaMKII_7NCspersyn_meanJIpersyn) size(JI_CaMKII_7NCspersyn_meanJIpersyn,1);];
pvdata = [mean(JI_PV_2NCspersyn_meanJIpersyn) std(JI_PV_2NCspersyn_meanJIpersyn) size(JI_PV_2NCspersyn_meanJIpersyn,1);
mean(JI_PV_3NCspersyn_meanJIpersyn) std(JI_PV_3NCspersyn_meanJIpersyn) size(JI_PV_3NCspersyn_meanJIpersyn,1);
mean(JI_PV_4NCspersyn_meanJIpersyn) std(JI_PV_4NCspersyn_meanJIpersyn) size(JI_PV_4NCspersyn_meanJIpersyn,1);
mean(JI_PV_5NCspersyn_meanJIpersyn) std(JI_PV_5NCspersyn_meanJIpersyn) size(JI_PV_5NCspersyn_meanJIpersyn,1);
mean(JI_PV_6NCspersyn_meanJIpersyn) std(JI_PV_6NCspersyn_meanJIpersyn) size(JI_PV_6NCspersyn_meanJIpersyn,1);
mean(JI_PV_7NCspersyn_meanJIpersyn) std(JI_PV_7NCspersyn_meanJIpersyn) size(JI_PV_7NCspersyn_meanJIpersyn,1);];
alldata = [camkiidata pvdata]

camkiidata_CV = [std(JI_CaMKII_2NCspersyn_meanJIpersyn)/mean(JI_CaMKII_2NCspersyn_meanJIpersyn)
std(JI_CaMKII_3NCspersyn_meanJIpersyn)/mean(JI_CaMKII_3NCspersyn_meanJIpersyn)
std(JI_CaMKII_4NCspersyn_meanJIpersyn)/mean(JI_CaMKII_4NCspersyn_meanJIpersyn)
std(JI_CaMKII_5NCspersyn_meanJIpersyn)/mean(JI_CaMKII_5NCspersyn_meanJIpersyn)
std(JI_CaMKII_6NCspersyn_meanJIpersyn)/mean(JI_CaMKII_6NCspersyn_meanJIpersyn)
std(JI_CaMKII_7NCspersyn_meanJIpersyn)/mean(JI_CaMKII_7NCspersyn_meanJIpersyn)];
pvdata_CV = [std(JI_PV_2NCspersyn_meanJIpersyn)/mean(JI_PV_2NCspersyn_meanJIpersyn)
std(JI_PV_3NCspersyn_meanJIpersyn)/mean(JI_PV_3NCspersyn_meanJIpersyn)
std(JI_PV_4NCspersyn_meanJIpersyn)/mean(JI_PV_4NCspersyn_meanJIpersyn)
std(JI_PV_5NCspersyn_meanJIpersyn)/mean(JI_PV_5NCspersyn_meanJIpersyn)
std(JI_PV_6NCspersyn_meanJIpersyn)/mean(JI_PV_6NCspersyn_meanJIpersyn)
std(JI_PV_7NCspersyn_meanJIpersyn)/mean(JI_PV_7NCspersyn_meanJIpersyn)];
alldata_CV = [camkiidata_CV pvdata_CV]


%%
function out = JaccardIndex(InputData)
    if size(InputData,1) == 1
        out = [1 1 1 1 1 -1000];
    elseif size(InputData,1) > 1
        NCsOfGivenSynapse = InputData(:,8:19);
        NCcenters = InputData(:,5:7);
        out = [];
        NCsOfGivenSynapse(NCsOfGivenSynapse<0)=0;
        Combinations = nchoosek([1:size(NCsOfGivenSynapse,1)],2);
        for i = 1:size(Combinations,1)
            out(i,1:2) = Combinations(i,:);
            thiscombo = [NCsOfGivenSynapse(out(i,1),:);NCsOfGivenSynapse(out(i,2),:)];
            thiscombosum = sum(thiscombo);
            JI_numerator = sum(thiscombosum==2);
            JI_denominator = sum(thiscombosum>0);
            thisJI = JI_numerator/JI_denominator;
            out(i,3) = thisJI;
            out(i,4:5) = Combinations(i,:);
            out(i,6) = pdist2(NCcenters(out(i,4),:),NCcenters(out(i,5),:)).*160;
        end
    end
    out(isnan(out))=0;
end


%%
function [A,B,C] = errorFunc(S,varargin)
warning(S.identifier, S.message);
A = NaN;
B = NaN;
C = NaN;
end
