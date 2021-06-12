% From Burke Rosen
% Modified by Yusi to loop over individuals according to subjectID


%function ATLAS_HCP_dMRIparcel_to_IC100node_connectivity
% Maps HCP rs-fMRI IC100 nodes to HCP MMP1.0 parcels, 
% then construct node-to-node dMRI connectivity matrix by taking weighted sum of
% parcel-to-parcel connectivities within each node. 


cd('~/Documents/fMRI_Real/HCP/DTI_Burke/HCP_IC100_dMRIconn/')
%% Gather dependencies
if isempty(which('ft_read_cifti'))
  addpath('./cifti-matlab')
end

%% Load HCPMMP1.0 parcel indices in grayordinate space
% grayordinates are surface-based vertices, 32,492 per hemipshere, 64,984 total
grayOrdParcel = ft_read_cifti(...
  './Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
grayOrdParcel = grayOrdParcel.indexmax; 

%% Load HCP rs-fMRI IC nodes in brainordinate space
% brainordinates are grayordinates + non-corical voxel (volumetric locations)
brainOrdNode = ft_read_cifti(...
  './Mnet1.pconn.nii');
brainOrdNode = brainOrdNode.brainordinate.parcellation;

%% isolate cortical nodel labels (those in grayordinate space)
% the first 64,984 brainordinates are grayordinates
grayOrdNode = brainOrdNode(1:64984);

%%% trivial case spoof for QC
% % make nodes 1 & 2 spatially equivalent to parcel 1 & 2 
% % Then node connectivity output should equal parcel connectivity input
% grayOrdNode = zeros(64984,1);
% grayOrdNode(grayOrdParcel == 1) = 1; 
% grayOrdNode(grayOrdParcel == 2) = 2;

%% load dMRI connectivity
% Note: the code I have written here is for the average connectivity, but with
% minor modification would work the same for individuals, looping over subject
dMRI = load('./averageConnectivity_Fpt.mat');

dMRI_ind = load('../individualConnectivity_10^Fpt');

%% For each node-to-node connection create weighted sum of parcelwise dMRI connectivity 
nN = 100; % number of nodes
nodeConn = zeros(nN,nN);
nodePairs = nchoosek(1:nN,2);%[src targ] (oneway)
nP = size(nodePairs,1);
for iP = 1:nP
    srcParc = grayOrdParcel(grayOrdNode == nodePairs(iP,1));
    targParc = grayOrdParcel(grayOrdNode == nodePairs(iP,2));
    
    % sanity checks
    if isempty(srcParc) || isempty(targParc)
      % warning('source or target node is entirely non-cortical, skipping pair ...')
      continue
    end
    if any(isnan(srcParc)) || any(isnan(targParc))
      error('node grayordinate does not have parcel label!')
    end
    
    nodeConn(nodePairs(iP,1),nodePairs(iP,2)) = ...
      mean(dMRI.Fpt(srcParc,targParc),'all','omitnan');
      % nans must be omitted because sometimes the source and target nodes share
      % grayordinates of the same parcel
end
% symmetrize (dMRI connectivity is always symetric)
nodeConn = nodeConn + nodeConn';% no need to divide by 2 as original is upper only
nodeConn(nodeConn == 0) = nan;
% Note: nan's indicate that one or both of the nodes is entirely non-cortical,
%       they are also along the (trivial) diagonal


%% 
% Create a dictionary of cortical or subcortical component


parcelIDs = [dMRI.parcelIDs;'non-cortical'];
nodeParc = repmat(struct(),nN,1);
for iN = 1:nN %loop over node
  nodePidx = grayOrdParcel(grayOrdNode == iN);
  if isempty(nodePidx) || numel(nodePidx) < sum(brainOrdNode == iN)
    % add 'non-cortical' parcel to node pool
    nodePidx = [nodePidx;repmat(361,sum(brainOrdNode == iN),1)]; %#ok<AGROW>
  end
  uniPidx = unique(nodePidx);

  % count grayorinates of each parcel, including non-cortical
  cntPidx = histc(nodePidx,uniPidx); %#ok<HISTC>
  [cntPidx,srti] = sort(cntPidx,'descend');
  uniPidx = uniPidx(srti);

  % assemble output struct
  nodeParc(iN).parcels = parcelIDs(uniPidx);
  nodeParc(iN).parcelIdxs = uniPidx;
  nodeParc(iN).parcelFraction = cntPidx./sum(cntPidx);
  nodeParc(iN).parcelGrayordinateCount = cntPidx;
end
completelyCorticalNodes = find(cellfun(@(x) ~ismember(361,x),{nodeParc.parcelIdxs}));
save('HCP_IC100_nodeParcels.mat','nodeParc','parcelIDs','completelyCorticalNodes')

Corticality = [];
for i = 1:nN
    if ~ismember(361,nodeParc(i).parcelIdxs)
        Corticality(i) = 1;
        continue
    else
        idx_tmp = find(nodeParc(i).parcelIdxs==361);
        Corticality(i) = 1 - nodeParc(i).parcelFraction(idx_tmp);
    end
end


% Loop over individuals according to subjectID used in fMRI_Real analysis
dMRI_ind = load('../individualConnectivity_10^Fpt');
tmp = load('~/Documents/fMRI_Real/HCP/HCP_PTN1200/Analysis/Timeseries_IC100.mat','SID_list');
SID_list = tmp.SID_list;
dMRI_ICA100_list = {};

for s_idx = 1:length(SID_list)
  SID = str2num(SID_list{s_idx})
  loc_idx = find(dMRI_ind.subjectIDs==SID);
  if isempty(loc_idx)
    continue
  end
  DTI_Parcel = dMRI_ind.Fpt(:,:,loc_idx);
  nN = 100; 
  nodeConn = zeros(nN,nN);
  nodePairs = nchoosek(1:nN,2);%[src targ] (oneway)
  nP = size(nodePairs,1);
  for iP = 1:nP
      srcParc = grayOrdParcel(grayOrdNode == nodePairs(iP,1));
      targParc = grayOrdParcel(grayOrdNode == nodePairs(iP,2));
      
      % sanity checks
      if isempty(srcParc) || isempty(targParc)
        % warning('source or target node is entirely non-cortical, skipping pair ...')
        continue
      end
      if any(isnan(srcParc)) || any(isnan(targParc))
        error('node grayordinate does not have parcel label!')
      end
      
      nodeConn(nodePairs(iP,1),nodePairs(iP,2)) = ...
        mean(DTI_Parcel(srcParc,targParc),'all','omitnan');
        % nans must be omitted because sometimes the source and target nodes share
        % grayordinates of the same parcel
  end
  % symmetrize (dMRI connectivity is always symetric)
  nodeConn = nodeConn + nodeConn';% no need to divide by 2 as original is upper only
  nodeConn(nodeConn == 0) = nan;
  % Note: nan's indicate that one or both of the nodes is entirely non-cortical,
  %       they are also along the (trivial) diagonal
  dMRI_ICA100_list{s_idx} = nodeConn;
end

save('HCP_dMRI_IC100_ind.mat','dMRI_ICA100_list','SID_list','Corticality')














