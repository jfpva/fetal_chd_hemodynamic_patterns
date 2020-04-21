% script_prep_data
%
% Prepare fetal blood flow and oxygen saturation data for analysis.
%
% See README.md for additional information.


%% Define Data Files

xlsxInFilePath = 'data/data.xlsx';  % NOTE: this file not included in online repo
csvDataFilePath = 'data/data.csv';
csvFlowFilePath = 'data/flow.csv';


%% Extract Data from Excel File

% Load relevant data

variableNames = { 'Group', 'SubGroup', 'CaseNo', 'CirculationType', 'GestationalAge', 'DateOfMRI', ...
                  'EFW', 'MPAFlow', 'AAoFlow', 'SVCFlow', 'DAFlow', 'DAoFlow', 'PBFFlow', 'UVFlow', ...
                  'AAoSO2', 'UVSO2', 'DAoSO2', 'MPASO2', 'SVCSO2' };
              
opts = detectImportOptions( xlsxInFilePath );

opts.SelectedVariableNames = variableNames;
              
T = readtable( xlsxInFilePath, opts );

% Remove empty rows

tf = all( ismissing(T), 2 );
T(tf,:) = [];

% Write selected data to csv

writetable( T, csvDataFilePath )


%% Extract Measured Flows

% Initialize Table of Measured Values

M = table(T.Group,T.SubGroup,T.CaseNo,T.CirculationType,T.MPAFlow,T.AAoFlow,T.SVCFlow,T.DAFlow,...
          T.DAoFlow,T.PBFFlow,T.UVFlow,'VariableNames',...
          {'Group';'SubGroup';'CaseNo';'CirculationType';'MPA';'AAo';'SVC';'DA';'DAo';'PBF';'UV'});  

% Remove Reported Flows That Appear to be Derived

% DA
isDerivedDA = ( T.DAFlow == T.DAoFlow - T.AAoFlow + T.SVCFlow );
M.DA( isDerivedDA ) = NaN;

% PBF
isDerivedPBF = ( T.PBFFlow == T.MPAFlow - T.DAFlow );  % NOTE: if Q_PBF == Q_MPA + Q_DA, it's assumed that Q_PBF and not Q_DA or Q_MPA
M.PBF( isDerivedPBF ) = NaN;

% Display Changes
descStr = sprintf( 'Reported Derived Flows' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
fprintf( 'DA\n--\n\n%i cases where reported Q_DA = Q_DAo - Q_AAo + Q_SVC\n\n', sum(isDerivedDA) )
if sum(isDerivedDA) > 0 
    fprintf( 'Changing Q_DA to NaN.\n\n' )
    disp(T(isDerivedDA,{'SubGroup','CaseNo','MPAFlow','AAoFlow','SVCFlow','DAFlow','DAoFlow','PBFFlow','UVFlow'}))
end
fprintf( 'PBF\n---\n\n%i cases where reported Q_PBF = Q_MPA - Q_DA\n\n', sum(isDerivedPBF) )
if sum(isDerivedPBF) > 0 
    fprintf( 'Changing Q_PBF to NaN.\n\n' )
    disp(T(isDerivedPBF,{'SubGroup','CaseNo','MPAFlow','AAoFlow','SVCFlow','DAFlow','DAoFlow','PBFFlow','UVFlow'}))
end

% Write Measured Data to CSV

writetable( M, csvFlowFilePath )


%% Extract Oxygen Saturations for Each Vessel and Save Formatted for Prism

% Output

outputDirPath = 'results';
outputPrismDirPath = fullfile( outputDirPath, 'prism' );

% SaO2 by Group
vesselNames = {'MPASO2','AAoSO2','SVCSO2','DAoSO2','UVSO2'};
for iV = 1:numel(vesselNames)
    Normal  = T.(vesselNames{iV})(strcmp(T.Group,'Normal'));
    HLHS    = T.(vesselNames{iV})(strcmp(T.Group,'HLHS'));
    TGA     = T.(vesselNames{iV})(strcmp(T.Group,'TGA'));
    TOF     = T.(vesselNames{iV})(strcmp(T.Group,'TOF'));
    EA      = T.(vesselNames{iV})(strcmp(T.Group,'EA'));
    TA      = T.(vesselNames{iV})(strcmp(T.Group,'TA'));
    maxGroupSize = max( [ numel(Normal), numel(HLHS), numel(TGA), numel(TOF), numel(EA), numel(TA) ] );
    Normal((end+1):maxGroupSize) = NaN;
    HLHS((end+1):maxGroupSize) = NaN;
    TGA((end+1):maxGroupSize) = NaN;
    TOF((end+1):maxGroupSize) = NaN;
    EA((end+1):maxGroupSize) = NaN;
    TA((end+1):maxGroupSize) = NaN;
    S = table(Normal,HLHS,TGA,TOF,EA,TA);
    outputFileDirPath = fullfile(outputPrismDirPath,sprintf('sao2_group_%s.csv',lower(strrep(vesselNames{iV},'SO2',''))));
    writetable( S,outputFileDirPath );
end

