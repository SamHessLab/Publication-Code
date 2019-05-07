function avg_density=density_avg_1colorG_mask_func1B(xf_um1,yf_um1)

box_size=0.02; % size of grid pixels for density plot, units of microns
xmin=min(xf_um1);  % value for x coordinate corresponding to left side of grid
ymin=min(yf_um1);  % value for y coordinate corresponding to bottom of grid
grid_width_pixels=1500;  % grid width for estimation of average density (not the same as for clustering)
density_grid1=zeros(grid_width_pixels,grid_width_pixels);  % this is the density grid itself

for i=1:length(xf_um1)
   xi=xf_um1(i);
   yi=yf_um1(i);
   row_number=floor((xi-xmin)/box_size)+1;          % round the coordinates to nearest pixel
   col_number=floor((yi-ymin)/box_size)+1;
   if row_number>0 && row_number<grid_width_pixels && col_number>0 && col_number<grid_width_pixels
     density_grid1(row_number,col_number)=density_grid1(row_number,col_number)+1;       % increment the grid pixel
   end
end
dens1_max=2;        % any pixels with more than this number of molecules (localizations) will be treated as equal 
density_grid1_norm=(density_grid1.*double(density_grid1<dens1_max)+dens1_max*double(density_grid1>=dens1_max));
       
density=zeros(grid_width_pixels,grid_width_pixels);
density(:,:)=density_grid1_norm;

density_is_nonzero=double(density_grid1_norm>1);
kw=50;          % create disk-like kernel for convolution
% thus any pixels with nonzero density will generate a convolved image with a circle of radius kw 
% and therefore any regions within a distance kw*box_size of a localization will be positive

[X,Y]=meshgrid(-kw:1:kw,-kw:1:kw);
R=double(X.*X+Y.*Y);
kern=double(R<kw);
density_conv=conv2(density_is_nonzero,kern,'same');     % do the convolution
density_mask1=double(density_conv~=0);
density_mask2=imdilate(density_mask1,kern);         % now do several rounds of dilation 
% to fill in any hollow regions within generally positive zones
for j=1:2       % this value can be increased to do more filling in, if desired
    density_mask2=imdilate(density_mask2,kern);
end
np_total_grid1=sum(density_grid1(density_mask2>0));     % find the number of localizations within all positive zones
area_total_grid1=sum(density_mask2(:));             % find the total area of positive zones
avg_density=np_total_grid1/(area_total_grid1*box_size*box_size);   % take the ratio to estimate the average density

