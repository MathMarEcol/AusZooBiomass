clear
close all
re_plot = 0;

dir = '/Volumes/WD_MYBOOK/Satellite/Climatology/MODIS';
% 
% %% SEA SURFACE TEMPERATURE
% 
% [files,shortfiles] = Dir2FileName('*SST_sst*.nc',dir);
% 
% for a = 1:length(files)
%     
%     if a == 1
%         nc.lat = ncread(files{a},'lat');
%         nc.lon = ncread(files{a},'lon');
%     end
%     
%     if a == 1
%         csvwrite('Lat.csv',nc.lat);
%         csvwrite('Lon.csv',nc.lon');
%     end
%     
%     nc.var = ncread(files{a},'sst');
%     csvwrite([shortfiles{a}(1:end-2),'csv'],nc.var');
%     
%     if re_plot == 1
%         %% Plotting
%         h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
%         p = pcolorm(nc.lat,nc.lon,nc.var');
%         setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);
%         
%         title(shortfiles{a},'Interpreter','none')
%         colormap jet
%         colorbar
%         caxis([0 35])
%         set(gcf,'color','w')
%         
%         eval(['export_fig ',['',shortfiles{a}(1:end-3)],' -pdf'])
%         close
%     end
%     nc = rmfield(nc,'var');
%     
% end
% 
% 
% %%CHLOROPHYLL A
% 
% [files,shortfiles] = Dir2FileName('*CHL_chlor_a*.nc',dir);
% 
% for a = 1:length(files)
%     nc.var = ncread(files{a},'chlor_a');    
%     csvwrite([shortfiles{a}(1:end-2),'csv'],nc.var');
%     
%     if re_plot == 1
%         %% Plotting
%         h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
%         p = pcolorm(nc.lat,nc.lon,log10(nc.var)');
%         setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);
%         
%         title(shortfiles{a},'Interpreter','none')
%         colormap jet
%         caxis([-2 1])
%         clbr = colorbar;
%         clbr.Ticks = [-2 -1 0 1];
%         clbr.TickLabels = {'0.01','0.1','1','10'};
%         set(gcf,'color','w')
%         eval(['export_fig ',['',shortfiles{a}(1:end-3)],' -pdf'])
%         close
%     end
%     nc = rmfield(nc,'var');
%     
% end
% 
% clear fi_lat fi_lon nc
% 
% %% NOW DO SST CLIMATOLOGY
% 
% nc.lat = csvread('Lat.csv');
% nc.lon = csvread('Lon.csv');
% 
% !rm Month00*
% [files,shortfiles] = Dir2FileName('*SST_sst*.csv',dir);
% 
% for a = 1:length(files)
%     if a == 1
%         temp = csvread(files{a});
%         SST = ones(size(temp,1),size(temp,2),12);
%         clear temp
%     end
%     
%     SST(:,:,a) = csvread(files{a});
% end
% 
% Clim = nanmean(SST,3); clear SST
% csvwrite(['Month00',shortfiles{a}(8:end-3),'csv'],Clim);
% 
% if re_plot == 1
%     %% Plotting
%     h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
%     p = pcolorm(nc.lat,nc.lon,Clim');
%     setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);
%     
%     title(['Clim',shortfiles{a}(8:end-2)],'Interpreter','none')
%     colormap jet
%     colorbar
%     caxis([0 35])
%     set(gcf,'color','w')
%     
%     eval(['export_fig ',['Month00',shortfiles{a}(8:end-4)],' -pdf'])
%     close
% end
% 
% clear Clim
% 
% 
% 
% %% NOW DO Chl CLIMATOLOGY
% 
% [files,shortfiles] = Dir2FileName('*CHL_chlor_a*.csv',dir);
% for a = 1:length(files)
%     
%     if a == 1
%         temp = csvread(files{a});
%         Chl = ones(size(temp,1),size(temp,2),12);
%         clear temp
%     end
%     
%     Chl(:,:,a) = csvread(files{a});
%     
% end
% 
% Clim = nanmean(Chl,3); clear Chl
% csvwrite(['Month00',shortfiles{a}(8:end-3),'csv'],Clim);
% 
% if re_plot == 1
%     %% Plotting
%     h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
%     p = pcolorm(nc.lat,nc.lon,log10(Clim)');
%     setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);
%     
%     title(['Clim',shortfiles{a}(8:end-2)],'Interpreter','none')
%     colormap jet
%     caxis([-2 1])
%     clbr = colorbar;
%     clbr.Ticks = [-2 -1 0 1];
%     clbr.TickLabels = {'0.01','0.1','1','10'};
%     
%     set(gcf,'color','w')
%     
%     eval(['export_fig ',['Month00',shortfiles{a}(8:end-4)],' -pdf'])
%     close
% end
% clear Chl
% 
% 

%% NOW DO BATHYMETRY
bathy_lat = single(ncread('GEBCO_2014_2D.nc','lat'));
bathy_lon = single(ncread('GEBCO_2014_2D.nc','lon'));

Depth = single(ncread('GEBCO_2014_2D.nc','elevation'));
Depth(Depth > 10) = NaN;
Depth = fliplr(Depth);

%%
bathy_lon = repmat(bathy_lon,1,length(bathy_lat));
bathy_lat = repmat(bathy_lat',length(bathy_lon),1);

lat = csvread('Lat.csv');
lon = csvread('Lon.csv')';
  
nc.lon = single(repmat(lon,1,length(lat)));
nc.lat = single(repmat(lat',length(lon),1));

%% NEW
% Try this before running again
% https://au.mathworks.com/help/matlab/math/resample-image-with-gridded-interpolation.html

% And a post from Ant
% https://gis.stackexchange.com/questions/255150/using-resample-vs-aggregate-extend-in-r-to-have-rasters-of-matching-resolutio

rlon = single(repmat(lon',length(lat),1));
rlat = single(repmat(lat,1,length(lon)));

nc.Depth = nc.lon.*NaN;

% A quick and dirty calculation
rsDepth = (Depth(1:5:end,1:5:end))';
csvwrite('Bathymetry.csv',rsDepth);
res = 10;
csvwrite('Bathymetry_LR.csv',rsDepth(1:res:end,1:res:end));

csvwrite('Lat_Grid.csv',rlat);
csvwrite('Lon_Grid.csv',rlon);


csvwrite('Lat_Grid_LR.csv',rlat(1:res:end,1:res:end));
csvwrite('Lon_Grid_LR.csv',rlon(1:res:end,1:res:end));



for i = 5:9
    eval(['chl = csvread(''Month0',num2str(i),'.L3m_MC_CHL_chlor_a_4km.csv'');'])
    eval(['csvwrite(''Month0',num2str(i),'.L3m_MC_CHL_chlor_a_4km_LR.csv'',chl(1:res:end,1:res:end));'])
    
    eval(['sst = csvread(''Month0',num2str(i),'.L3m_MC_SST_sst_4km.csv'');'])
    eval(['csvwrite(''Month0',num2str(i),'.L3m_MC_SST_sst_4km_LR.csv'',sst(1:res:end,1:res:end));'])
    clear sst chl
end


for i = 10:12
    eval(['chl = csvread(''Month',num2str(i),'.L3m_MC_CHL_chlor_a_4km.csv'');'])
    eval(['csvwrite(''Month',num2str(i),'.L3m_MC_CHL_chlor_a_4km_LR.csv'',chl(1:res:end,1:res:end));'])
    
    eval(['sst = csvread(''Month',num2str(i),'.L3m_MC_SST_sst_4km.csv'');'])
    eval(['csvwrite(''Month',num2str(i),'.L3m_MC_SST_sst_4km_LR.csv'',sst(1:res:end,1:res:end));'])
    clear sst chl
end



return

q = interp2(bathy_lon,bathy_lat,Depth,nc.lon,nc.lat);



for i = 1:size(nc.lon,1)
    fi_lon = dsearchn(bathy_lon(:,1),lon(i));    
    for ii = 1:size(nc.lon,2)
        fi_lat = dsearchn(bathy_lat(1,:)',lat(ii));
        
        nc.Depth(i,ii) = Depth(fi_lon,fi_lat);
                
    end
end


%% OLD
% F = scatteredInterpolant(bathy_lon(:),bathy_lat(:),Depth(:));
% nc.Depth = F(nc.lon2(:),nc.lat2(:));

% scatter(nc.lon2(:),nc.lat2(:),10,nc.Depth,'filled'); colorbar; caxis([-2000 0])

% Bathy = reshape(nc.Depth,size(nc.lon2,1),size(nc.lon2,2));
csvwrite('Bathymetry.csv',Bathy);








%%

% Depth = ones(length(bathy_lon),length(bathy_lat)).*NaN;
% for i = 1:length(bathy_lon)
%     for ii = 1:length(bathy_lat)
%     Depth(i,ii) = GEBCO_GetDepth(bathy_lat(ii),bathy_lon(i));
%     end
% end


% bathy_lon = repmat(bathy_lon,1,length(bathy_lat));
% bathy_lat = repmat(bathy_lat',length(bathy_lon),1);

% Depth = ncread('GEBCO_2014_2D.nc','elevation',[ix_lon(1)-1 ix_lat(1)-1],[length(ix_lon) length(ix_lat)]);
% Depth(Depth >0) = NaN;
%
%
% nc.lon2 = repmat(nc.lon,1,length(nc.lat));
% nc.lat2 = repmat(nc.lat',length(nc.lon),1);
%
% F = scatteredInterpolant(bathy_lon(:),bathy_lat(:),Depth(:));
%
% Bathy = F(nc.lon2(:),nc.lat2(:));
%
% Bathy2 = reshape(Bathy,size(nc.lon2,1),size(nc.lon2,2));


if re_plot == 1
    
    % Plot
    h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
    p = pcolorm(nc.lat,nc.lon,Bathy2);
    setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);
    
    title('GEBCO Bathymetry');
    colormap jet
    colorbar
    caxis([-3000 0])
    set(gcf,'color','w')
    
    export_fig Bathy -pdf
    close
    
end


%% CARS

DOY = dayofyear([2017,1,15,0,0,0;...
    2017,2,15,0,0,0;...
    2017,3,15,0,0,0;...
    2017,4,15,0,0,0;...
    2017,5,15,0,0,0;...
    2017,6,15,0,0,0;...
    2017,7,15,0,0,0;...
    2017,8,15,0,0,0;...
    2017,9,15,0,0,0;...
    2017,10,15,0,0,0;...
    2017,11,15,0,0,0;...
    2017,12,15,0,0,0]);

nc_lat = csvread('Lat.csv');
nc_lon = csvread('Lon.csv');

[mesh_lon,mesh_lat] = meshgrid(nc_lon,nc_lat);

var = {'nitrate','silicate','phosphate'};

for v = 1:length(var)
    for a = 1:length(DOY)
        
        out = CARS_Extract(var{v},Lat,Lon,DOY(a),0);
        
        out.var2 = interp2(out.lon,out.lat,out.var,mesh_lon,mesh_lat);
        
        csvwrite(['CARS_',var{v},'_Month',sprintf('%02d',a),'.csv'],out.var2);
        
        if re_plot == 1
            %% Plotting
            h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
            p = pcolorm(mesh_lat,mesh_lon,out.var2);
            setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);
            
            title(['CARS_',var{v},'_Month',sprintf('%02d',a)],'Interpreter','none')
            colormap jet
            colorbar
            
            if v == 1
                caxis([0 20])
            elseif v == 2
                caxis([0 10])
            elseif v == 3
                caxis([0 2])
            end
            
            uistack(h.f, 'top')
            
            set(gcf,'color','w')
            
            eval(['export_fig ',['CARS_',var{v},'_Month',sprintf('%02d',a)],' -pdf'])
            close
            
        end
        clear out
    end
end
