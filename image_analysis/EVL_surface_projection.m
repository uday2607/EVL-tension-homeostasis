% EVL_surface_projection.m
% Paper: Hino, Kapoor, Gubbala, Hannezo & Heisenberg
%        "Apical domain mechanosensation regulates tissue tension homeostasis"
%        Institute of Science and Technology Austria (ISTA)
%
% Generates apical surface-projected images from a confocal z-stack TIFF.
% See image_analysis/README.md for parameter descriptions.

clear
clc
File_name = 'zstack_t1.tif';
Surface_ch = 2;
Target_ch =1;
offset=0;
DeltaZ=6;
Th_value = 4;

Mode=1;%mode=1, thresholding method; mode0, differential method
G_blur=10;%select filtor size. if 0 skip blurring. normally 10 works
Height_blur=0;

t_size=1;
z_size=60;
c_size=2;


tiff_info = imfinfo(File_name); % return tiff structure, one element per image
tiff_stack = imread(File_name, 1) ; % read first image
surf_tiff_stack = tiff_stack;
if G_blur~=0
    surf_tiff_stack =imgaussfilt(surf_tiff_stack,G_blur);
end


%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    ii
    temp_tiff = imread(File_name, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);%image order; xyczt
    if G_blur~=0
        surf_temp_tiff =imgaussfilt(temp_tiff,G_blur);
        surf_tiff_stack = cat(3 , surf_tiff_stack, surf_temp_tiff);
    else
        surf_tiff_stack = cat(3 , surf_tiff_stack, temp_tiff);
    end
end

x_size=length(tiff_stack(1,:,1,1,1));
y_size=length(tiff_stack(:,1,1,1,1));

tiff_stack=reshape(tiff_stack,y_size,x_size,c_size,z_size,t_size);
surf_tiff_stack = reshape(surf_tiff_stack,y_size,x_size,c_size,z_size,t_size);
Height = zeros(size(tiff_stack(:,:,1,1,t_size)));
Output = zeros(size(tiff_stack(:,:,1,1,t_size)));
Surface_output = zeros(size(tiff_stack(:,:,1,1,t_size)));
%Surface_mean_output = zeros(size(tiff_stack(:,:,1,1,t_size)));
Surface_surface = zeros(size(tiff_stack(:,:,1,1,t_size)));

for t=1:t_size
    t
    for y=1:y_size
        for x=1:x_size
            Height(y,x,1,1,t)=min([z_size, find(surf_tiff_stack(y,x,Surface_ch,:,t)>=Th_value, 1 )]);
        end
    end
end

if Height_blur~=0
    Height_smooth =imgaussfilt(Height(:,:,1,1,1),Height_blur);      
for t=2:t_size
    Height_temp_tiff =imgaussfilt(Height(:,:,1,1,t),Height_blur);
    Height_smooth = cat(5 , Height_smooth, Height_temp_tiff);
end
Height=round(Height_smooth);
Height(Height<0)=0;
Height(Height>255)=255;
end


 for t=1:t_size
    t
    for y=1:y_size
        for x=1:x_size           
            temp_hight=Height(y,x,1,1,t);
            if temp_hight+offset+DeltaZ>z_size
                top_posi=z_size;
            else
                top_posi=temp_hight+offset+DeltaZ;
            end
            if temp_hight+offset-DeltaZ<1
                bottom_posi=1;
            elseif temp_hight+offset+DeltaZ>z_size
                bottom_posi=z_size;
            else
                bottom_posi=temp_hight+offset-DeltaZ;
            end
            Output(y,x,1,1,t)=max(tiff_stack(y,x,Target_ch,bottom_posi:top_posi,t));
            Surface_output(y,x,1,1,t)=max(tiff_stack(y,x,Surface_ch,bottom_posi:top_posi,t));
            %Surface_mean_output(y,x,1,1,t)=mean(tiff_stack(y,x,Surface_ch,bottom_posi:top_posi,t));
            Surface_surface(y,x,1,1,t)=tiff_stack(y,x,Surface_ch,temp_hight,t);
        end
    end
end





imwrite(Height(:,:,:,:,1)/255,'Hight_map.tif');
imwrite(Output(:,:,:,:,1)/255,['Target_surface_max_off',num2str(offset),'_delta',num2str(DeltaZ),'_th',num2str(Th_value),'.tif']);
imwrite(Surface_output(:,:,:,:,1)/255,['Reference_surface_max_off',num2str(offset),'_delta',num2str(DeltaZ),'_th',num2str(Th_value),'.tif']);
%imwrite(Surface_mean_output(:,:,:,:,1)/255,['Reference_surface_mean_off',num2str(offset),'_delta',num2str(DeltaZ),'_th',num2str(Th_value),'.tif']);
imwrite(Surface_surface(:,:,:,:,1)/255,'Reference_surface.tif');


for t=2:t_size
    imwrite(Height(:,:,:,:,t)/255,'Hight_map.tif','WriteMode','append');
    imwrite(Output(:,:,:,:,t)/255,['Target_surface_max_off',num2str(offset),'_delta',num2str(DeltaZ),'_th',num2str(Th_value),'.tif'],'WriteMode','append');
    imwrite(Surface_output(:,:,:,:,t)/255,['Reference_surface_max_off',num2str(offset),'_delta',num2str(DeltaZ),'_th',num2str(Th_value),'.tif'],'WriteMode','append');
    %imwrite(Surface_mean_output(:,:,:,:,t)/255,['Reference_surface_mean_off',num2str(offset),'_delta',num2str(DeltaZ),'_th',num2str(Th_value),'.tif'],'WriteMode','append');
    imwrite(Surface_surface(:,:,:,:,t)/255,'Reference_surface.tif','WriteMode','append');
end