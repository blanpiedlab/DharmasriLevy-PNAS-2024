%% pairing SQ detected synapses from different channels using centroid coordinates
% Script to process file structure output from synapseIDandmeasure.ijm to
% remove outliers based on area and pair pre/post
% "Differential nanoscale organization of excitatory synapses onto excitatory vs. inhibitory neurons"
% Poorna A Dharmasri, Aaron D Levy, Thomas A Blanpied, PNAS 2024
% authors: Aaron Levy, Sarah Metzbower
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine

clear

%want to make this flexible in terms of what channel is what so creating
%variable so that you can define your channels at the top rather than rely
%on the imaging order.

% ChA
% ChB-this is going to be your primary synapse marker
% ChC


%first need to organize and re-number data so that for each SQ channel its
%only the data from that channel as well as remove outliers

%loop parental folders
xcol=6;
ycol=7;
int_col=8;
parent = dir('*Week*');
flags = [parent.isdir];
fols = parent(flags);
fols = {fols.name};
fols = fols(~ismember(fols,{'.' '..'}));
parentpath = parent(1).folder;

for parental = 1:length(fols)

    %loop results folders- comment out lines n-m and n-m if you want to run in a single
    %folder
    thisweek = fols{parental};
    weekpath = fullfile(parentpath,thisweek);
    resultsfs = dir(fullfile(weekpath,'*_Results'));
    resultsflags = [resultsfs.isdir];
    resultsfols = resultsfs(resultsflags);
    resultsfols = {resultsfols.name};
    resultsfols = resultsfols(~ismember(resultsfols,{'.''..'}));

    for q=1:length(resultsfols)

        thisresultsfolder = resultsfols{q};
        thisresultspath = fullfile(weekpath,thisresultsfolder);

        Cond = {'PreResults*'; 'PostResults*'};
        protein= {'Pre'; 'Post'};
        cond= {'ChA'; 'ChB'}; %changing the order here will assign which channel is which

        if exist(fullfile(thisresultspath,'ImMask.zip'))

            for i = 1:length(Cond)

                D = dir(fullfile(thisresultspath,Cond{i}));
                xx = fullfile(D.folder,D.name);
                alldata = readmatrix(xx,'NumHeaderLines',1);

                areacutoff = 0.5;
                Ch_xx = alldata(alldata(:,2)<areacutoff,:);

                savename = fullfile(thisresultspath,[cond{i}, '_OR', '.txt']);
                save(savename,'Ch_xx', '-ascii', '-tabs');

            end

            cutoff=0.500; %define how far apart the puncta can be and still be matched

            %load data files
            Condx = {'*ChA_OR*'; '*ChB_OR*'};

            %Condx= {'*ChA_OR*'; '*ChB_OR*'; '*ChC_OR*'};
            for h=1:length(Condx)
                
                G = dir(fullfile(thisresultspath,Condx{h}));
                gg = fullfile(G.folder,G.name);
                load(gg);
            
            end

            %need to specify which channels you want to compare

            A=1; %how +/- um from point to define search area

            ROIkey=NaN(length(ChB_OR(:,1)),1);

            ROI_key=[ChB_OR(:,1) ROIkey];

            for s=1:length(ChB_OR)

                ROIs=ChB_OR(s,1);
                p=[ChB_OR(s,xcol),ChB_OR(s,ycol)];
                pxmax=(p(1,1)+A);
                pxmin=(p(1,1)-A);
                pymax=(p(1,2)+A);
                pymin=(p(1,2)-A);
                z=find(ChA_OR(:,xcol)<=pxmax & ChA_OR(:,xcol)>=pxmin & ChA_OR(:,ycol)>=pymin & ChA_OR(:,ycol)<=pymax);

                if length(z)>=1
                    zz=1;
                    pq=[ChA_OR(z(zz),xcol),ChA_OR(z(zz),ycol)];
                    dd=pdist2(p,pq);
                    mind=min(dd);
                    minx=find((dd)==mind);
                    if length(z)>=2
                        for zz=2:length(z)
                            pq=[ChA_OR(z(zz),xcol),ChA_OR(z(zz),ycol)];
                            d=pdist2(p,pq);
                            dd=[dd; d];
                        end
                        mind=min(dd);
                        minx=find((dd)==mind);
                        if mind<=cutoff

                            pp=[ChA_OR(z(minx),xcol),ChA_OR(z(minx),ycol)];
                            ppxmax=(pp(1,1)+A);
                            ppxmin=(pp(1,1)-A);
                            ppymax=(pp(1,2)+A);
                            ppymin=(pp(1,2)-A);

                            zd=find(ChB_OR(:,xcol)<=pxmax & ChB_OR(:,xcol)>=pxmin & ChB_OR(:,ycol)>=pymin & ChB_OR(:,ycol)<=pymax);

                            rs=[ChB_OR(zd(1),xcol),ChB_OR(zd(1),ycol)];
                            testd=pdist2(pp,rs);
                            ddd=testd;
                            if length(zd)>=2
                                for ab=2:length(zd)

                                    rs=[ChB_OR(zd(ab),xcol),ChB_OR(zd(ab),ycol)];
                                    testd=pdist2(pp,rs);
                                    ddd=[ddd; testd];
                                end
                            else
                            end
                            mind_t=min(ddd);
                            if mind_t==mind
                                ROI_key(s,2)=ChA_OR(z(minx),1);
                                ChA_OR(z(minx),1:8)=NaN;
                            else
                            end

                        else
                        end
                    else
                        if dd<=cutoff
                            pp=[ChA_OR(z(zz),xcol),ChA_OR(z(zz),ycol)];
                            ppxmax=(pp(1,1)+A);
                            ppxmin=(pp(1,1)-A);
                            ppymax=(pp(1,2)+A);
                            ppymin=(pp(1,2)-A);

                            zd=find(ChB_OR(:,xcol)<=pxmax & ChB_OR(:,xcol)>=pxmin & ChB_OR(:,ycol)>=pymin & ChB_OR(:,ycol)<=pymax);

                            rs=[ChB_OR(zd(1),xcol),ChB_OR(zd(1),ycol)];
                            testd=pdist2(pp,rs);
                            ddd=testd;

                            if length(zd)>=2
                                for ab=2:length(zd)

                                    rs=[ChB_OR(zd(ab),xcol),ChB_OR(zd(ab),ycol)];
                                    testd=pdist2(pp,rs);
                                    ddd=[ddd; testd];

                                end
                            else


                            end
                            mind_t=min(ddd);

                            if mind_t==dd
                                ROI_key(s,2)=ChA_OR(z(minx),1);
                                ChA_OR(z(minx),1:8)=NaN;
                            else
                            end

                        end
                    end
                else
                end
            end

            %switch to next channel

            savename = fullfile(thisresultspath,['ROI_key.txt']);
            save(savename,'ROI_key', '-ascii', '-tabs');
            Condx= {'*ChA_OR*'; '*ChB_OR*'};

            %Condx= {'*ChA_OR*'; '*ChB_OR*'; '*ChC_OR*'};
            for h=1:length(Condx)
                G = dir(fullfile(thisresultspath,Condx{h}));
                gg = fullfile(G.folder,G.name);
                load(gg);
            end


            c=figure('Color', 'white'); hold on; axis equal;
            plot(ChB_OR(:,xcol), ChB_OR(:,ycol),  'Color', 'blue',...
                'Markersize', 7,...
                'LineStyle', 'none',...
                'MarkerFaceColor', 'none',...
                'Marker', 'o');
            plot(ChA_OR(:,xcol), ChA_OR(:,ycol),  'Color', 'red',...
                'Markersize', 6,...
                'LineStyle', 'none',...
                'MarkerFaceColor', 'none',...
                'Marker', 'o');
            saveas(gcf,fullfile(thisresultspath,'scatter_allsynapses'));
            close(c)

            load(fullfile(thisresultspath,'ROI_key.txt'));

            ROI_zeros=ROI_key;

            ROI_zeros(isnan(ROI_zeros))=0;

            ROI_cat=ROI_zeros;

            ROI=ROI_zeros(:,1);

            ROI_ChB=ROI_cat(:,1);
            ROI_ChB(ROI_ChB>=1)=1;

            ROI_ChA=ROI_cat(:,2);
            ROI_ChA(ROI_ChA>=1)=2;


            ROI_sum=ROI_ChA + ROI_ChB;

            ROI_ChB_cat=[ROI ROI_sum];



            colors = {'r', 'm', 'g', 'k', 'c', 'y', 'b'};

            c=figure('Color', 'white'); hold on; axis equal;

            for n=1:length(ChB_OR(:,1))
                nn=ROI_sum(n);
                plot(ChB_OR(n,xcol), ChB_OR(n,ycol), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'Color', colors{nn})
                hold on
            end


            saveas(gcf,fullfile(thisresultspath,'scatter_categorized_ChB_test'));
            close(c)

        end % processing
        
    end % results folder

end % parent
