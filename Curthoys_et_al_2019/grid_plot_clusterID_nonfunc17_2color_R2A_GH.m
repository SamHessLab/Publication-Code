clear
save_cell_output=0;         % flag to write cluster parameter output as a .mat file 
save_overall_output=1;
outfile_overall='output.mat'; % overall output file
                    % default location is within the same folder as this script
save_figure=0;  

min_densityHA_prelim=1000;   % HA density threshold (#/um2) for clustering; not used if rel2avg=1
min_densityPIP2_prelim=1000;   % PIP2 density threshold (#/um2) for clustering; not used if rel2avg=1
grid_width=35;  % width of square grid in micrometers
np_norm=0;     % if nonzero, scale the weight of each localization by a linear factor
      % equal to the ratio of this number to the actual number of HA localizations in the dataset
      % this option is not recommended, but can help normalize datasets of different sizes
loc_unc=40;  % localization uncertainty (in nanometers) used for blurring the density plot before cluster identification
npn=0;   % if nonzero, normalizes dataset as though this # of localizations (each lo calization counts proportionally to this # / actual #)
np_trunc=0;  % if nonzero, truncates dataset at this number of localizations (for each species separately)
avg_density_ratio=4;  % if nonzero, sets density threshold to this # times the average density over the cell
avg_density_fixed=0;  % if nonzero, forces the average density (HA/um2) to be this value for every cell
rel2avg=1;      % flag to select cluster identification relative to the average density
dx=10;         % grid pixel size in nanometers
titletext='HAwt1';
inpath='/Volumes/NASA4TB1/2018_07_28 Dendra2HA-MUTANTS + PAmKatePH/';
infile='*.mat';

close all

dx_um=dx/1000;
gw=grid_width/dx_um;
loc_unc_um=loc_unc/1000;

loc_unc_pix=loc_unc/dx;
kernw=ceil(loc_unc_pix)+1;
[xk,yk]=meshgrid(-kernw:1:kernw,-kernw:1:kernw);
rk2=xk.*xk+yk.*yk;
kern=exp(-2*rk2/(loc_unc_pix*loc_unc_pix));  % prepare a kernel for convolution with the grid plot
kernn=kern/sum(kern(:));         % normalize kernel to total area of one

inpathsearch=strcat(inpath,infile);
filelist = dir(inpathsearch);
nfiles=length(filelist);
ncells=nfiles;                  % this analysis assumes each file within the directory corresponds to one cell
area_values_cell_HA_HA=cell(ncells,1);
density_values_cell_HA_HA=cell(ncells,1);
density_values_cell_HA_PIP2=cell(ncells,1);
np_values_cell_HA_HA=cell(ncells,1);
np_values_cell_HA_PIP2=cell(ncells,1);
perim_values_cell_HA_HA=cell(ncells,1);

average_density_cell_HA=cell(ncells,1);

area_values_cell_PIP2_PIP2=cell(ncells,1);
density_values_cell_PIP2_PIP2=cell(ncells,1);
density_values_cell_PIP2_HA=cell(ncells,1);
np_values_cell_PIP2_HA=cell(ncells,1);
np_values_cell_PIP2_PIP2=cell(ncells,1);
perim_values_cell_PIP2_PIP2=cell(ncells,1);

average_density_cell_PIP2=cell(ncells,1);

nloc=0;
%    sz = get(0,'ScreenSize');               % uncomment these two lines to make the figure fill the whole screen
%    figure('Position',[1 sz(4) sz(3) sz(4)]);

npix_high_HA_w_PIP2=zeros(ncells,1);
npix_high_HA=zeros(ncells,1);
npix_high_PIP2=zeros(ncells,1);
   

for fileloop=1:nfiles    % load and analyze all ".mat" files within the directory specified above
     
    infile=strcat(inpath,filelist(fileloop).name);
    outfile=strcat(infile,'_ColorClus.fig');            % output file for figure generated
    clear q xf_all yf_all xf_um yf_um xf_nd yf_nd
    load(infile,'xf_all','yf_all','q','nrat','lp','r0_all');    % loads localizations using standard Hess lab format
    % this format contains the following variables:
    % positions (xf_all,yf_all) are the (x,y) locations of each molecule encoded in pixel values 
    % q is the number of micrometers per pixel
    % positions (xf_um,yf_um) are the (x,y) locations of each molecule encoded in micrometers
    % a0_all is the fitted amplitude of the peak of the Gaussian
    % "color" of molecular species is encoded in the variable nrat
    % nrat is the ratio of (# of photons in the red channel)/(total photons in both channels)
    %    and is equivalent to the variable alpha defined in Gunewardene et al. Biophysical Journal 2011
    % variable lp is the calculated localization precision for each molecule (in micrometers)
    % variable r0_all is the fitted width of the Gaussian for each molecule 
    
    if exist('xf_all','var')==0
        load(infile,'x_nd','y_nd');
        if exist('x_nd','var')==0
            load(infile,'xf_um','yf_um');
            q=-1;
        else
            xf_um=x_nd;
            yf_um=y_nd;
            q=-1;
        end
    end
    if exist('q','var')==0
        load(infile,'um_per_pixel');
        q=um_per_pixel;
    end
    if q>0
        xf_um=xf_all*q;
        yf_um=yf_all*q;
    end
    xf_min=min(xf_um);           % calculate minimum value of x and y position from dataset 
    yf_min=min(yf_um);          % and use later to subtract from dataset to maximize use of grid
    if xf_min<0 
        xf_min=0;
    end
    if yf_min<0
        yf_min=0;
    end
    
    xf_um0=xf_um; yf_um0=yf_um; % save the old localizations until they are identified as HA or PIP2
    subset1=find(nrat<0.6 & lp>0.005 & lp<0.1 & r0_all>0.7 & r0_all<3);         % Species 1 (HA)
    subset2=find(nrat>0.65 & lp>0.005 & lp<0.1 & r0_all>0.7 & r0_all<3);        % Species 2 (PIP2)
    % in addition to selecting the two species based on nrat (alpha) values,
    % these boolean conditions essentially apply a simple tolerance function
    % which rejects localizations having too small (<5 nm) or too large (>100 nm) a position uncertainty
    % or too large (>700 nm) or small (<100 nm) a width for their PSF
    
    xf_um_HA=xf_um0(subset1);
    yf_um_HA=yf_um0(subset1);
    xf_um_PIP2=xf_um0(subset2);
    yf_um_PIP2=yf_um0(subset2);
    
    npHA=length(xf_um_HA);
    npPIP2=length(xf_um_PIP2);
    np0=npHA;  % save original number of localizations
    
    nloc=nloc+npHA;
    if np_norm>0            % implementation of "weighted grid plot" 
        % where (if flag enabled) localizations are counted more heavily (in the proportion np_norm/npHA)
        % when the total number of localizations is lower than np_norm
        % however, in general, the threshold relative to average cell density seems to better reduce cell to cell variability
        normfac=np_norm/npHA;
    else
        normfac=1;
    end
 
    average_densityHA=density_avg_1colorG_mask_func1B(xf_um_HA,yf_um_HA);   % measure the average density of HA in the cell
    average_densityPIP2=density_avg_1colorG_mask_func1B(xf_um_PIP2,yf_um_PIP2);  % measure the average density of PIP2 in the cell
    average_density_cell_HA{fileloop}=average_densityHA;
    average_density_cell_PIP2{fileloop}=average_densityPIP2;
    
    xf_um_HAn=xf_um_HA-xf_min;          % try to center the HA and PIP2 localizations within the grid
    yf_um_HAn=yf_um_HA-yf_min;
    xf_um_PIP2n=xf_um_PIP2-xf_min;          
    yf_um_PIP2n=yf_um_PIP2-yf_min;

    if avg_density_ratio==0
        min_densityHA_final=min_densityHA_prelim;
    else
        min_densityHA_final=average_densityHA*avg_density_ratio;
    end
    if avg_density_fixed>0
        np_trunc1=floor(np0*avg_density_fixed/average_densityHA);  % truncate localizations so that target avg density is met
        if npHA>np_trunc1 
            xf_um_HA=xf_um_HA(1:np_trunc1);
            yf_um_HA=yf_um_HA(1:np_trunc1);
            npHA=np_trunc1;
        else
            error_message=strcat('insufficient localizations to match avg density . ',infile)
        end
        
    else
        if npHA>np_trunc && np_trunc>0             % if np_trunc is given, truncate the dataset to a max of np_trunc points
            xf_um_HA=xf_um_HA(1:np_trunc);
            yf_um_HA=yf_um_HA(1:np_trunc);
            npHA=np_trunc;
        end
        
    end
 
    if avg_density_ratio==0
        min_densityPIP2_final=min_densityPIP2_prelim;
    else
        min_densityPIP2_final=average_densityPIP2*avg_density_ratio;
    end

    
    grid1=zeros(gw,gw);         % density grid for HA
    grid2=zeros(gw,gw);         % density grid for PIP2
    for i=1:npHA
        gj=round(xf_um_HAn(i)/dx_um);
        gi=round(yf_um_HAn(i)/dx_um);
        if gi>0 && gj>0 && gi<=gw && gj<=gw
            grid1(gi,gj)=grid1(gi,gj)+normfac;
        end
    end
    for i=1:npPIP2
        gj=round(xf_um_PIP2n(i)/dx_um);
        gi=round(yf_um_PIP2n(i)/dx_um);
        if gi>0 && gj>0 && gi<=gw && gj<=gw
            grid2(gi,gj)=grid2(gi,gj)+normfac;
        end
    end
    gridconv1=conv2(grid1,kernn,'same')/(dx_um*dx_um);
    gridconv2=conv2(grid2,kernn,'same')/(dx_um*dx_um);

    gridconvBW_HA=(gridconv1>=min_densityHA_final);
    gridconvBW_PIP2=(gridconv2>=min_densityPIP2_final);
    gridconv_size=size(gridconvBW_HA);
    gridconvColor=gridconvBW_HA*0;
    
    % begin process of identifying clusters based on boolean array gridconvBW_HA 
    % gridconvBW_HA has values: (0= below threshold; 1= above threshold)
    
    sHA  = regionprops(gridconvBW_HA, 'Area','PixelIdxList','Perimeter');
    area_all_HA_HA = cat(1, sHA.Area)*dx_um*dx_um;
    perim_all_HA_HA = cat(1, sHA.Perimeter)*dx_um;
    nclusHA=length(area_all_HA_HA);
    density_all_HA_HA=0*area_all_HA_HA;
    density_all_HA_PIP2=0*area_all_HA_HA;
    np_all_HA_HA=0*area_all_HA_HA;
    np_all_HA_PIP2=0*area_all_HA_HA;

    sPIP2  = regionprops(gridconvBW_PIP2, 'Area','PixelIdxList','Perimeter');
    % find all contiguous regions within the gridconvBW_HA array and measure their properties
    area_all_PIP2_PIP2 = cat(1, sPIP2.Area)*dx_um*dx_um;
    perim_all_PIP2_PIP2 = cat(1, sPIP2.Perimeter)*dx_um;
    nclusPIP2=length(area_all_PIP2_PIP2);
    density_all_PIP2_PIP2=0*area_all_PIP2_PIP2;
    density_all_PIP2_HA=0*area_all_PIP2_PIP2;
    np_all_PIP2_PIP2=0*area_all_PIP2_PIP2;
    np_all_PIP2_HA=0*area_all_PIP2_PIP2;

    for i=1:nclusHA         % calculate density of, and number of localizations within each cluster (for HA)
       indi= sHA(i).PixelIdxList;
       HA_density_iclus=gridconv1(indi);
       PIP2_density_iclus=gridconv2(indi);
       mean_HA_density_iclus=mean(HA_density_iclus);
       mean_PIP2_density_iclus=mean(PIP2_density_iclus);
       density_all_HA_HA(i)=1000*mean_HA_density_iclus/average_densityHA;          % in HA clusters, density of HA
       density_all_HA_PIP2(i)=1000*mean_PIP2_density_iclus/average_densityPIP2;      % in HA clusters, density of PIP2
       np_all_HA_HA(i)=sum(grid1(indi));                % in HA clusters, # of HAs in each cluster
       np_all_HA_PIP2(i)=sum(grid2(indi));              % in HA clusters, # of PIP2s in each cluster
       gridconvColor(indi)=mean_HA_density_iclus;        % create visualization of average density for each cluster
   end
   for i=1:nclusPIP2        % calculate density of, and number of localizations within each cluster (for PIP2)
       indi= sPIP2(i).PixelIdxList;
       HA_density_iclus=gridconv1(indi);
       PIP2_density_iclus=gridconv2(indi);
       mean_HA_density_iclus=mean(HA_density_iclus);
       mean_PIP2_density_iclus=mean(PIP2_density_iclus);
       density_all_PIP2_HA(i)=1000*mean_HA_density_iclus/average_densityHA;          % in PIP2 clusters, density of HA
       density_all_PIP2_PIP2(i)=1000*mean_PIP2_density_iclus/average_densityPIP2;      % in PIP2 clusters, density of PIP2
       np_all_PIP2_HA(i)=sum(grid1(indi));                % in PIP2 clusters, # of HAs in each cluster
       np_all_PIP2_PIP2(i)=sum(grid2(indi));              % in PIP2 clusters, # of PIP2s in each cluster
   end
   gridconvColor3=zeros(gridconv_size(1),gridconv_size(2),3);          % this will be a green (HA) - magenta (PIP2) merge
   gridconvColor3(:,:,1)=gridconvBW_PIP2;
   gridconvColor3(:,:,2)=gridconvBW_HA;
   gridconvColor3(:,:,3)=gridconvBW_PIP2;
   npix_high_HA_w_PIP2(fileloop)=sum(sum(gridconvBW_HA.*gridconvBW_PIP2));
   npix_high_HA(fileloop)=sum(sum(gridconvBW_HA));
   npix_high_PIP2(fileloop)=sum(sum(gridconvBW_PIP2));
   
    % uncomment these next 6 lines to show visual output of cluster identification (density encoded by color)
    %     imagesc(gridconvColor);          
    %     cmap1=hsv(256);
    %     cmap1(1,:)=[0 0 0];
    %     colormap(cmap1);
    %     caxis([1500 5000]);
    %     colorbar
    
   imagesc(gridconvColor3);        % show a color merge of HA(green) and PIP2(magenta) clustering and spatial overlap
    
   axis image
   titletext2=strcat(titletext,'  Avg Density=',num2str(average_densityHA),'  npHA=',num2str(npHA),'  npPIP2=',num2str(npPIP2));
   title(titletext2);
   drawnow
   
   area_values_cell_HA_HA{fileloop}=area_all_HA_HA;      % save values of area, density, number of localizations per cluster
   perim_values_cell_HA_HA{fileloop}=perim_all_HA_HA;       % in several cell arrays, one cell-array-element for each cell(file)
   area_values_cell_PIP2_PIP2{fileloop}=area_all_PIP2_PIP2;
   perim_values_cell_PIP2_PIP2{fileloop}=perim_all_PIP2_PIP2;
   density_values_cell_HA_HA{fileloop}=density_all_HA_HA;
   density_values_cell_HA_PIP2{fileloop}=density_all_HA_PIP2;
   np_values_cell_HA_HA{fileloop}=np_all_HA_HA;
   np_values_cell_HA_PIP2{fileloop}=np_all_HA_PIP2;
   density_values_cell_PIP2_PIP2{fileloop}=density_all_PIP2_PIP2;
   density_values_cell_PIP2_HA{fileloop}=density_all_PIP2_HA;
   np_values_cell_PIP2_PIP2{fileloop}=np_all_PIP2_PIP2;
   np_values_cell_PIP2_HA{fileloop}=np_all_PIP2_HA;
 
   if save_figure==1
        saveas(gcf,outfile_figure,'fig');
   end
   outfile_this_cell=strcat('cell',num2str(fileloop),'.mat');
   if save_cell_output==1
      save(outfile_this_cell);
   end
end
frac_high_HA_w_PIP2=npix_high_HA_w_PIP2./npix_high_HA;          % calculate fraction of pixels with high HA that also have high PIP2
frac_high_PIP2_w_HA=npix_high_HA_w_PIP2./npix_high_PIP2;       % calculate fraction of pixels with high PIP2 that also have high HA
save(outfile_overall);

