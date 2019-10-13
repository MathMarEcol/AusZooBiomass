

clear
close all

Lat = [-55 -8];
Lon = [105 178];

dir = '/Users/jason/UNSW/Projects/ZooplanktonDP/Projects/BiomassPaper/SatelliteCoVariate/';


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
    
nc_lat = csvread('SatelliteCoVariate/Lat.csv');
nc_lon = csvread('SatelliteCoVariate/Lon.csv');

[mesh_lon,mesh_lat] = meshgrid(nc_lon,nc_lat);

var = {'nitrate','silicate','phosphate'};

for v = 2:length(var)
    for a = 1:length(DOY)
        
        out = CARS_Extract(var{v},Lat,Lon,DOY(a),0);
        
        out.var2 = interp2(out.lon,out.lat,out.var,mesh_lon,mesh_lat);
                
        csvwrite(['SatelliteCoVariate/CARS_',var{v},'_Month',sprintf('%02d',a),'.csv'],out.var2);
        
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
        
        eval(['export_fig ',['SatelliteCoVariate/CARS_',var{v},'_Month',sprintf('%02d',a)],' -pdf'])
        close
        
        clear out
    end
end


return
clear


%% SEA SURFACE TEMPERATURE

[files,shortfiles] = Dir2FileName('*SST_sst*.nc',dir);

for a = 1:length(files)

    if a == 1
        nc.lat = ncread(files{a},'lat');
        nc.lon = ncread(files{a},'lon');
        
        fi_lat = find(nc.lat >= Lat(1) & nc.lat <= Lat(2));
        fi_lon = find(nc.lon >= Lon(1) & nc.lon <= Lon(2));
        
        nc.lat = nc.lat(fi_lat);
        nc.lon = nc.lon(fi_lon);
    end
    
    if a == 1
        csvwrite('SatelliteCoVariate/Lat.csv',nc.lat);
        csvwrite('SatelliteCoVariate/Lon.csv',nc.lon');
    end

    nc.var = ncread(files{a},'sst');
    nc.var = nc.var(fi_lon,fi_lat);

    csvwrite([files{a}(1:end-2),'csv'],nc.var');

    %% Plotting
    h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
    p = pcolorm(nc.lat,nc.lon,nc.var');
    setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);

    title(shortfiles{a},'Interpreter','none')
    colormap jet
    colorbar
    caxis([0 35])
    set(gcf,'color','w')

    eval(['export_fig ',['SatelliteCoVariate/',shortfiles{a}(1:end-3)],' -pdf'])
    close

    nc = rmfield(nc,'var');

end


%%CHLOROPHYLL A

[files,shortfiles] = Dir2FileName('*CHL_chlor_a*.nc',dir);

for a = 1:length(files)
    nc.var = ncread(files{a},'chlor_a');
    nc.var = nc.var(fi_lon,fi_lat);
    
    csvwrite([files{a}(1:end-2),'csv'],nc.var');
    
    %% Plotting
    h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
    p = pcolorm(nc.lat,nc.lon,log10(nc.var)');
    setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);
    
    title(shortfiles{a},'Interpreter','none')
    colormap jet
    caxis([-2 1])
    clbr = colorbar;
    clbr.Ticks = [-2 -1 0 1];
    clbr.TickLabels = {'0.01','0.1','1','10'};    
    set(gcf,'color','w')
    eval(['export_fig ',['SatelliteCoVariate/',shortfiles{a}(1:end-3)],' -pdf'])
    close
    
    nc = rmfield(nc,'var');
    
end

clear fi_lat fi_lon nc

%% NOW DO SST CLIMATOLOGY

nc.lat = csvread('SatelliteCoVariate/Lat.csv');
nc.lon = csvread('SatelliteCoVariate/Lon.csv');

!rm SatelliteCoVariate/Month00*
[files,shortfiles] = Dir2FileName('*SST_sst*.csv',dir);

for a = 1:length(files)
    if a == 1
        temp = csvread(files{a});
        SST = ones(size(temp,1),size(temp,2),12);
        clear temp
    end
    
    SST(:,:,a) = csvread(files{a});
end

Clim = nanmean(SST,3); clear SST
csvwrite(['SatelliteCoVariate/Month00',shortfiles{a}(8:end-3),'csv'],Clim);

%% Plotting
h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
p = pcolorm(nc.lat,nc.lon,Clim');
setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);

title(['Clim',shortfiles{a}(8:end-2)],'Interpreter','none')
colormap jet
colorbar
caxis([0 35])
set(gcf,'color','w')

eval(['export_fig ',['SatelliteCoVariate/Month00',shortfiles{a}(8:end-4)],' -pdf'])
close

clear Clim



%% NOW DO Chl CLIMATOLOGY

[files,shortfiles] = Dir2FileName('*CHL_chlor_a*.csv',dir);
for a = 1:length(files)
    
    if a == 1
        temp = csvread(files{a});
        Chl = ones(size(temp,1),size(temp,2),12);
        clear temp
    end
    
    Chl(:,:,a) = csvread(files{a});
    
end

Clim = nanmean(Chl,3); clear Chl
csvwrite(['SatelliteCoVariate/Month00',shortfiles{a}(8:end-3),'csv'],Clim);


%% Plotting
h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
p = pcolorm(nc.lat,nc.lon,log10(Clim)');
setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);

title(['Clim',shortfiles{a}(8:end-2)],'Interpreter','none')
colormap jet
caxis([-2 1])
clbr = colorbar;
clbr.Ticks = [-2 -1 0 1];
clbr.TickLabels = {'0.01','0.1','1','10'};

set(gcf,'color','w')

eval(['export_fig ',['SatelliteCoVariate/Month00',shortfiles{a}(8:end-4)],' -pdf'])
close

clear Chl

%% NOW DO BATHYMETRY

bathy_lat = ncread('GEBCO_2014_2D.nc','lat');
bathy_lon = ncread('GEBCO_2014_2D.nc','lon');
Depth = ncread('GEBCO_2014_2D.nc','elevation');

ix_lat = find(bathy_lat>=Lat(1) & bathy_lat <= Lat(2));
ix_lon = find(bathy_lon>=Lon(1) & bathy_lon <= Lon(2));

bathy_lat = bathy_lat(ix_lat);
bathy_lon = bathy_lon(ix_lon);

Depth = Depth(ix_lon,ix_lat);

%% Reduce resoltion to half;
bathy_lat = bathy_lat(1:end);
bathy_lon = bathy_lon(1:end);
Depth = Depth(1:end,1:end);

Depth(Depth > 0) = NaN;

%% 
bathy_lon = repmat(bathy_lon,1,length(bathy_lat));
bathy_lat = repmat(bathy_lat',length(bathy_lon),1);

nc.lon2 = repmat(nc.lon,1,length(nc.lat));
nc.lat2 = repmat(nc.lat',length(nc.lon),1);

F = scatteredInterpolant(bathy_lon(:),bathy_lat(:),Depth(:));
nc.Depth = F(nc.lon2(:),nc.lat2(:));
scatter(nc.lon2(:),nc.lat2(:),10,nc.Depth,'filled'); colorbar; caxis([-2000 0])

Bathy = reshape(nc.Depth,size(nc.lon2,1),size(nc.lon2,2));
csvwrite('SatelliteCoVariate/Bathymetry.csv',Bathy);


return
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







% Plot
h = plot_EastCoast(Lon(1),Lon(2),Lat(1),Lat(2),2,0);
p = pcolorm(nc.lat,nc.lon,Bathy2);
setm(h.h,'mlinelocation',10,'plinelocation',10,'mlabellocation',10,'plabellocation',10);

title('GEBCO Bathymetry');
colormap jet
colorbar
caxis([-3000 0])
set(gcf,'color','w')

export_fig SatelliteCoVariate/Bathy -pdf
close






