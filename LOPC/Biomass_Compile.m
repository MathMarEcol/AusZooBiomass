clear
close all

load("/Users/jason/UNSW/Projects/ZooplanktonDP/Projects/MissingLink/Output/MNF_Summary.mat")

% Remove some voyages which aren't ready or have collaborators I need to
% deal with.
fi = find(strcmp('IN2014e03', tbl.Voyage)==0 & ...
    strcmp('IN2016v01', tbl.Voyage)==0 & ...
    strcmp('IN2016v02', tbl.Voyage)==0 & ...
    strcmp('IN2019v02', tbl.Voyage)==0 & ...
    tbl.Biomass ~= 0);% & ... % Remove the Zeros (there shouldn't be many now)
%     tbl.Biomass < 10^5); % And the high values


tbl = tbl(fi,:);

tbl.datetimeUTC = datetime(tbl.datestr,'TimeZone','UTC');
tbl.Unit2 = repmat("LOPC",height(tbl),1);
tbl.Unit2(strcmp(tbl.Unit2,"OPC1T")==1 | strcmp(tbl.Unit2,"OPC2T")==1) = "OPC";

[zd,zltr,zone] = timezone(tbl.Lon);
tbl.datetimeLocal = tbl.datetimeUTC + zd/24;

bio.RECORD_ID = NaN .* tbl.datenum;
bio.PROJECT_ID = tbl.Voyage;
bio.SAMPLE_ID = tbl.cast_no;
bio.LONGITUDE = tbl.Lon;
bio.LATITUDE = tbl.Lat;
bio.SAMPLE_DATE_LOCAL = tbl.datetimeUTC;
bio.SAMPLE_DATE_UTC = tbl.datetimeLocal;
bio.SAMPLE_YEAR  = year(tbl.datetimeLocal);
bio.SAMPLE_MONTH = month(tbl.datetimeLocal);
bio.SAMPLE_DAY = day(tbl.datetimeLocal);
bio.DAY_NIGHT =  hour(tbl.datetimeLocal) + (minute(tbl.datetimeLocal) + second(tbl.datetimeLocal)/60)/60;    
bio.MIN_DEPTH = tbl.Depth_min;
bio.MAX_DEPTH = tbl.Depth_max;
bio.BIOMASS = tbl.Biomass;
bio.BIOMASS_UNITS = repmat("mg m-3 (Wet Weight)",height(tbl),1);
bio.BIOMASS_METHOD = repmat("Estimated Wet Weight from BioVolume",height(tbl),1);
bio.C_UGL = NaN .* tbl.datenum;
bio.REPLICATE = (0 .* tbl.datenum) + 1;
bio.TOW_TYPE = tbl.Orientation;
bio.NET_MESH_SIZE = repmat("min. ESD = 300 microns",height(tbl),1);
bio.NET_TYPE = tbl.Unit2;
bio.BOTTOM_DEPTH_ETOPO = NaN .* tbl.datenum;

bio = struct2table(bio);

writetable(bio,'../Data/MNF_LOPC_Biomass.csv');

    
