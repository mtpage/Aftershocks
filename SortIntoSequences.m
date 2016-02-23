function assignedCatalog = SortIntoSequences(catalog,mainshockIndices,startTime,endTime,numFaultLengths,minDist,maxDepth,minMag,maxMagDiff)
%SORTINTOSEQUENCES        Assign aftershocks to mainshocks.
%        ASSIGNEDCATALOG = SORTINTOSEQUENCES(CATALOG,MAINSHOCKINDICES,
%        STARTTIME,ENDTIME,NUMFAULTLENGTHS,MINDIST,MAXDEPTH,MINMAG,MAXMAGDIFF)
%        Returns the input CATALOG with an additional column giving
%        mainshock row index for all aftershocks of that
%        mainshock.  An earthquake is tagged as an aftershock of a
%        mainshock if it has a magnitude greater than MINMAG and
%        greater than the mainshock magnitude minus MAXMAGDIFF,
%        occurs within STARTTIME days and ENDTIME days of the
%        mainshock time, is at a depth above MAXDEPTH km, and is
%        within the maximum of NUMFAULTLENGTHS fault lengths and
%        MINDIST km of the mainshock.
%  
%        CATALOG is the input earthquake catalog entered in standard
%        10-column format.  It should be sorted in time from oldest
%        to newest events.  The columns are:
%        year/month/day/hour/minute/second/latitude/longitude/depth/magnitude
%  
%        MAINSHOCKINDICES gives the indices of rows of CATALOG
%        corresponding to mainshocks.
% 
%        This function returns ASSIGNEDCATALOG.  The first 10
%        columns are the same as the input CATALOG.  It has an
%        additional 11th column which is -1 if that earthquake is a
%        mainshock, 0 if that earthquake is unassigned, and is the
%        row index of the corresponding mainshock if that
%        earthquake is  an aftershock.
%
%        If an earthquake matches the criteria to be an aftershock
%        of more than one mainshock, it is assigned to the latest
%        mainshock.
%
%        Earthquake depths are not used to calculate interevent
%        distances, which are great circle distances.
%
%        Fault length is calculated from magnitude using average
%        Wells and Coppersmith (1994) scaling.  
%
%        Authors: Karen Felzer and Morgan Page
%                 U. S. Geological Survey
%        Last modified: May 2015

  
%  If an earthquake is within criteria to be an aftershock of more
%  than one mainshock, we will assign it to the latest mainshock,
%  to avoid double-counting.
%  This only works if catalog is sorted, so throw error if it's not.
if ~issorted(datenum(catalog(:,1:6)))
    error = MException('SortInfoSequences:Inputs','Catalog is not sorted from oldest to newest events.');
    throw(error)
end  
  
% Sanity check inputs: catalog and mainshockIndices 
if max(mainshockIndices)>size(catalog,1)
    error = MException('SortInfoSequences:Inputs','Some mainshockIndices greater than catalog length.');
    throw(error)
end  


% Initialize assignedCatalog.  All non-assigned earthquakes are tagged with 0's
assignedCatalog = [catalog zeros(length(catalog),1)];
% Mainshocks are tagged with -1
assignedCatalog(mainshockIndices,11)=-1;

for i=mainshockIndices'
        
    % The magnitude cut-off to use for aftershocks of this
    % mainshock is the larger of either minMag or the mainshock
    % magnitude - maxMagDiff 
    minMagThisMainshock = max([minMag catalog(i,10)-maxMagDiff]);
    
    % Find indices of possible aftershocks: earthquakes with the
    % time window [startTime endTime] relative to the mainshock,
    % above the minimum magnitude for the mainshock sequence.
    % Do not include eqs that are already tagged as mainshocks. 
    Tdiff = datenum(catalog(:,1:6)) - datenum(catalog(i,1:6));
    possibleAftInd = setdiff(find(Tdiff>startTime & Tdiff<=endTime & catalog(:,10)>=minMagThisMainshock & catalog(:,9)<=maxDepth),mainshockIndices);

    % Now check that those potential aftershocks are close enough
    % to the mainshock.  They can be within numFaultLengths or
    % minDist of the mainshock, whichever is larger.
    maxDistFromMainshock = max(minDist, numFaultLengths*WellsCoppersmithLength(catalog(i,10)));
    arclen = distance('gc',catalog(possibleAftInd,[7 8]),catalog(i,[7 8]),'deg');
    dist=deg2km(arclen);
    aftershockInd=possibleAftInd(find(dist<=maxDistFromMainshock));
     
    % Tag all aftershocks with the mainshock index 
    assignedCatalog(aftershockInd,11)=i;
   
end




function L = WellsCoppersmithLength(M)
%WELLSCOPPERSMITHLENGTH        Find fault length using Wells & Coppersmith relation
%        L = WELLSCOPPERSMITHLENGTH(M)
% Calculates the average length, in km, an earthquake of magnitude M.  
%
% Uses Wells & Coppersmith (1994) average equations (over all focal mechanisms).
logL = -2.44 + 0.59*M;
L = 10.^logL;

