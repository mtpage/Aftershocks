function mainshockIndices = FindMainshockIndices(catalog,minMainshockMag,maxMainshockMag,maxDepth,exclusionDistance,exclusionTimeBefore,exclusionTimeAfter,excludeDistanceFormat,excludeEarlyCatalog,excludeLateCatalog)
%FINDMAINSHOCKINDICES        Find indices of mainshocks using given exclusion criteria
%        MAINSHOCKINDICES = FINDMAINSHOCKINDICES(CATALOG,MINMAINSHOCKMAG,
%        MAXMAINSHOCKMAG,MAXDEPTH,EXCLUSIONDISTANCE,EXCLUSIONTIMEBEFORE,
%        EXCLUSIONTIMEAFTER,EXCLUDEDISTANCEFORMAT,EXCLUDEEARLYCATALOG, 
%        EXCLUDELATECATLOG)
%        Returns the row indices of mainshocks within CATALOG.
%        Mainshocks are earthquakes within the magnitude range
%        [MINMAINSHOCKMAG MAXMAINSHOCKMAG] and above MAXDEPTH km
%        for which there are no larger earthquakes within the EXCLUSIONDISTANCE
%        for EXCLUSIONTIMEBEFORE days prior and EXCLUSIONTIMEAFTER days after.
%
%        CATALOG is the earthquake catalog entered in standard
%        10-column format.  The columns are:
%        year/month/day/hour/minute/second/latitude/longitude/depth/magnitude
%
%        Earthquake depths are not used to calculate interevent
%        distances, which are great circle distances.
%  
%        If EXCLUDEDISTANCEFORMAT=1, EXCLUDEDISTANCE 
%        should be entered in fault lengths.  Each earthquake then
%        will exclude other earthquakes within EXCLUDEDISTANCE 
%        number of fault lengths, assuming average Wells & Coppersmith
%        magnitude-length scaling.
% 
%        Otherwise EXCLUDEDISTANCE should be entered in km (and
%        therefore is a constant distance for earthquakes of all magnitudes).
% 
%        If EXCLUDEEARLYCATALOG=1, mainshocks prior to the catalog
%        start date + EXCLUSIONTIMEBEFORE are excluded, since they
%        could potentially be excluded by earthquakes prior to the
%        start of the catalog.
%  
%        If EXCLUDELATECATALOG=1, mainshocks in the last 
%        EXCLUSIONTIMEAFTER days of the catalog are excluded, 
%        since they could potentially be excluded by earthquakes
%        after the end of the catalog.
%  
%        Authors: Morgan Page and Karen Felzer
%                 U. S. Geological Survey
%        Last modified: May 2015
  
  
  
%If excludeDistanceFormat = 1 then exclusionDistance should be entered
%in fault lengths, otherwise it should be entered in km
if(excludeDistanceFormat==1)
    exclusionDistanceKm = exclusionDistance*WellsCoppersmithLength(catalog(:,10));  
else
    exclusionDistanceKm = exclusionDistance*ones(size(catalog,1),1);
end

% Serial date of all eqs in catalog
T = datenum(catalog(:,1:6));

% Indices of earthquakes large enough to potentially exclude others
potentialExcluderInd=find(catalog(:,10)>=minMainshockMag)';

% Initially excluded indices -- earthquakes with magnitudes outside
% input magnitude window
excludedIndices=find(catalog(:,10)<minMainshockMag | catalog(:,10)>maxMainshockMag | catalog(:,9)>maxDepth);

% If excludeEarlyCatalog=1, then exclude earthquakes in first
% exclusionTimeBefore days of catalog
if excludeEarlyCatalog==1
  excludedIndices=union(excludedIndices,find(T<min(T)+exclusionTimeBefore));
end

% If excludeEarlyCatalog=1, then exclude mainshocks in last exclusionTimeAfter days of catalog as well 
if excludeLateCatalog==1
  excludedIndices=union(excludedIndices,find(T>max(T)-exclusionTimeAfter));
end

% Now we loop through all potentially exluding earthquakes to find more earthquakes to exclude as mainshocks
for i=potentialExcluderInd
  
    % Time difference, in days, between potential excluder and other eqs
    dT = T - T(i);  
    
    % Find great circle distance, in km, between potential excluder and other eqs    
    arclen = distance('gc',catalog(i,[7 8]),catalog(:,[7 8]),'deg');
    dist=deg2km(arclen);
      
    % Magnitude difference between potential excluder and other eqs
    dMag = catalog(:,10) - catalog(i,10);
    
    % Indices of smaller earthquakes within spatial/temporal exclusion criteria
    ind = find(dist<exclusionDistanceKm(i) & dT<exclusionTimeBefore & dT>-exclusionTimeAfter & dMag<0);
    
    % Add these to list of excluded earthquake indices
    excludedIndices=union(excludedIndices,ind);
    
end


% return non-excluded indices; these are the mainshock indices
mainshockIndices=setdiff(1:size(catalog,1),excludedIndices)';
    
    

function L = WellsCoppersmithLength(M)
%WELLSCOPPERSMITHLENGTH       Find fault length using Wells & Coppersmith relation
%        L = WELLSCOPPERSMITHLENGTH(M)
% Calculates the average length, in km, an earthquake of magnitude M.  
%
% Uses Wells & Coppersmith (1994) average equations (over all focal mechanisms).
logL = -2.44 + 0.59*M;
L = 10.^logL;



