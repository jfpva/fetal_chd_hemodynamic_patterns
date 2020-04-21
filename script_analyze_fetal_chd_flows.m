%% script_analyze_fetal_chd_flows.m
%
% Analyisis of fetal blood flow distribution patterns, measured using MRI, in late gestation human 
% fetuses with a spectrum of congenital heart disease (CHD) subtypes presenting with neonatal 
% cyanosis and age-matched fetuses with normal hearts.
%
% See README.md for additional information.

% NOTE: balance_fetal_flow_distribution function requires Matlab R2019b or later due to use of 
%       arguments code block for argument validation


%% Options

isLogResults = true;  % Log results flag    
isShowFigs = true;    % Create figures flag

dataFilePath   = fullfile( 'data', 'flow.csv' );  % Path to file with measured flows

minMeasToSummarize = 3;  % Minimum number of measurements per subgroup to summarize flows

objFn = 'wRMSD';  % Objective function used to model flow distributions


%% Setup

% Groups (CHD Sub-Types)
          
groups    = {   'Normal', ...
                'HLHS', ...
                'TGA', ...
                'TOF', ...
                'EA', ...
                'TA' };

groupToProcess = 1:numel(groups);

% Sub-Groups (Categories within CHD Sub-Types)

subGroups = {   'Normal', ...
                'HLHS RAS', ...
                'HLHS MA AS', ...
                'HLHS MS AS', ...
                'HLHS MS AA', ...
                'HLHS MA AA', ...
                'HLHS DORV', ...
                'TGA IVS', ...
                'TGA VSD', ...
                'TGA VSD PS', ...
                'TGA COA', ...
                'TOF ONLY', ...
                'TOF PA', ...
                'Ebstein''s no Circular Shunt', ...
                'Ebstein''s Circular Shunt', ...
                'TA VA Concordance', ...
                'TA VA Discordance' };

subGroupToProcess = 1:numel(subGroups);

% Output

outputDirPath = 'results';
outputDataDirPath = fullfile( outputDirPath, 'tables' );
outputLogDirPath = fullfile( outputDirPath, 'log' );
outputFigDirPath = fullfile( outputDirPath, 'figures' );
outputPrismDirPath = fullfile( outputDirPath, 'prism' );

logFilePath = fullfile( outputDirPath, 'log.txt' );


%% Create Directories as Required

if ~exist( outputDirPath, 'dir' )
    mkdir( outputDirPath )
end

if ~exist( outputDataDirPath, 'dir' )
    mkdir( outputDataDirPath )
end

if ~exist( outputFigDirPath, 'dir' )
    mkdir( outputFigDirPath )
end

if ~exist( outputLogDirPath, 'dir' )
    mkdir( outputLogDirPath )
end

if ~exist( outputPrismDirPath, 'dir' )
    mkdir( outputPrismDirPath )
end


%% Reset Log File

if isLogResults && exist( logFilePath, 'file' )
    delete( logFilePath );
end


%% Load Data to Table

[ ~, ~, dataFileExt ] = fileparts( dataFilePath );

switch dataFileExt
    
    case '.csv'   % Load from comma-separated values file
        
        M = readtable( dataFilePath );

    otherwise
    
        error( 'data file extensions (%s) not recognized', dataFileExt )
        
end

%% Define SubGroups with Pulmonary Atresia, Aortic Atresia or Ebstein's Anomaly with Circular Shunt

isPA   = strcmp( M.SubGroup, 'TOF PA' ) | strcmp( M.SubGroup, 'Ebstein''s no Circular Shunt' );
isAA   = strcmp( M.SubGroup, 'HLHS RAS' ) | strcmp( M.SubGroup, 'HLHS MS AA' ) | strcmp( M.SubGroup, 'HLHS MA AA' );
isEACS = strcmp( M.SubGroup, 'Ebstein''s Circular Shunt' );


%% Define Outlier Bounds

% Initialize Bounds Tables
L = M( M.CaseNo==1, {'SubGroup','MPA','AAo','SVC','DA','DAo','PBF','UV'} );
L{:,2:end} = nan(size(L{:,2:end}));
U = L;

% Determine Bounds
for iS = 1:numel(subGroups)
    % NOTE: isoutlier returns bounds as median +/- 3 scaled median absolute deviations, see isoutlier for details
    [~,lowerBound,upperBound] = isoutlier( M( strcmp(M.SubGroup,subGroups{iS}), {'MPA','AAo','SVC','DA','DAo','PBF','UV'} ) );
    for iV = 1:numel(lowerBound.Properties.VariableNames)
        vesselName = lowerBound.Properties.VariableNames{iV};
        nMeas = sum( strcmp(M.SubGroup,subGroups{iS}) & ~isnan( M.(vesselName) ) );
        if nMeas >= minMeasToSummarize
            L.(vesselName)(iS) = lowerBound.(vesselName);
            U.(vesselName)(iS) = upperBound.(vesselName);
        end
        switch vesselName
            case 'MPA'
                if strcmp( subGroups{iS}, 'Ebstein''s Circular Shunt' )
                    U.(vesselName)(iS) = min( [ 0, U.(vesselName)(iS) ] );
                else
                    L.(vesselName)(iS) = max( [ 0, L.(vesselName)(iS) ] );
                end
            case 'AAo'
                if strcmp( subGroups{iS}, 'HLHS RAS' ) || strcmp( subGroups{iS}, 'HLHS MS AA' ) || strcmp( subGroups{iS}, 'HLHS MA AA' )
                    U.(vesselName)(iS) = min( [ 0, U.(vesselName)(iS) ] );
                else
                    L.(vesselName)(iS) = max( [ 0, L.(vesselName)(iS) ] );
                end
             case 'DA'
                if strcmp( subGroups{iS}, 'TOF PA' ) || strcmp( subGroups{iS}, 'Ebstein''s no Circular Shunt' ) || strcmp( subGroups{iS}, 'Ebstein''s Circular Shunt' )
                    U.(vesselName)(iS) = min( [ 0, U.(vesselName)(iS) ] );
                end
            case {'DAo','SVC','PBF','UV'}
                L.(vesselName)(iS) = max( [ 0, L.(vesselName)(iS) ] );
        end
    end
end

% Save to Structure
Bounds = struct( 'lower', L, 'upper', U );

% Create New Bounds Tables Matched to Elements of M
L = M(:,{'SubGroup','CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV'});
indVess = find(strcmp(L.Properties.VariableNames,'MPA')):size(L,2);
L{:,indVess} = nan(size(L{:,indVess}));
U = L;
for iR = 1:size(L,1)
    L{iR,indVess} = Bounds.lower{strcmp(L.SubGroup(iR),Bounds.lower.SubGroup),2:end};
    U{iR,indVess} = Bounds.upper{strcmp(U.SubGroup(iR),Bounds.upper.SubGroup),2:end};
end

if isLogResults
    writetable( L, fullfile( outputDataDirPath, 'data_outlier_bounds_lower.csv' ) )
    writetable( U, fullfile( outputDataDirPath, 'data_outlier_bounds_upper.csv' ) )
    diary( logFilePath )
end

descStr = 'Outlier Bounds';
fprintf( '\n%s\n%s\n\n\n', descStr, repmat( '=', size(descStr) ) );

fprintf( 'Lower\n\n' )
display_table(Bounds.lower)
fprintf( '\n' )

fprintf( 'Upper\n\n' )
display_table(Bounds.upper)
fprintf( '\n' )

if isLogResults
    diary off
    clean_log_file( logFilePath )
end


%% Create Table of Derived Values

% Initialize Table
D = M;  % derived flows

% Add Derived Flow Variables
D.FO  = nan(size(D.MPA));
D.IVC = nan(size(D.MPA));
D.CA  = nan(size(D.MPA));
D.CS  = nan(size(D.MPA));
D.CVO = nan(size(D.MPA));


%% Set MPA for Pulmonary Atresia Cases

D.MPA(isPA==1) = 0;


%% Derive Missing Values 

isVerbose = false;

% From Measured Values
D = derive_flows( D, M, isAA, isPA, isEACS, L, U, isVerbose );

% Using Derived Values
D = derive_flows( D, D, isAA, isPA, isEACS, L, U, isVerbose );


% Output Results Command Window / Log
if isLogResults
    diary( logFilePath );
end
descStr = sprintf( 'Outlier Derived Flows' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
D = derive_flows( D, D, isAA, isPA, isEACS, L, U, true );
if isLogResults
    diary off
end


%% Calculate CVO, CA, CS

% CVO in Cases With PA
indPA = isPA | strcmp(D.Group,{'EA'});
D.CVO(indPA) = D.AAo(indPA) / 0.97;

% CVO in Cases With AA
indAA = isAA;
D.CVO(indAA) = D.MPA(indAA);

% CVO in All Other Cases
indNormal = ~indAA & ~indPA;
D.CVO(indNormal) = ( D.MPA(indNormal) + D.AAo(indNormal) ) / 0.97;

% CA, CS
D.CA  = 0.03*D.CVO;
D.CS  = D.CA;


%% Derive Values for Aortic Atresia Cases

% AAo
D.AAo(isAA) = -D.CA(isAA);

% Derive Additional Values Using New AAo Values
D = derive_flows( D, D, isAA, isPA, isEACS, L, U );


%% Derive FO

for iD = 1:size(D,1)
    switch D.SubGroup{iD}
        case {'Normal','HLHS MS AS'}
            D.FO(iD) = D.AAo(iD) + D.CA(iD) - D.PBF(iD);
        case {'TGA IVS','TGA COA'}
            D.FO(iD) = D.SVC(iD) + D.IVC(iD) - D.AAo(iD);
        case {'HLHS RAS','HLHS MS AA','HLHS MA AA'}
            D.FO(iD) = -D.PBF(iD);
        case {'TA VA Concordance','TA VA Discordance','Ebstein''s no Circular Shunt'}
            D.FO(iD) = D.SVC(iD) + D.IVC(iD) + D.CS(iD);
        case 'Ebstein''s Circular Shunt'
            D.FO(iD) = -D.MPA(iD) + D.SVC(iD) + D.IVC(iD) + D.CS(iD);
        case {'HLHS MA AS','HLHS DORV','TGA VSD','TGA VSD PS','TOF ONLY','TOF PA'}
            D.FO(iD) = NaN;
    end
end


%% Write Derived Data to CSV

if isLogResults
    writetable( D, fullfile( outputDataDirPath, 'data_derived.csv' ) )
end


%% Summarize Missing CVO Values

if isLogResults
    diary( logFilePath );
end

descStr = sprintf( 'Missing CVO Values' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
fprintf( '    Cannot derive CVO flow in %i cases.\n\n', sum(isnan(D.CVO)) )
display_table( D(isnan(D.CVO),{'Group','SubGroup','CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','CVO'}) )
fprintf( '\n\n\n' )

if isLogResults
    diary off
    clean_log_file( logFilePath )
end


%% Summarize Sub-Groups

S = struct( 'SubGroup', '', 'NumCases', [], 'CirculationType', cell(0), 'Tmeas', table(), 'Smeas', table(), 'Tcalc', table(), 'Scalc', table(), 'TmeasCvo', table(), 'SmeasCvo', table(), 'TcalcCvo', table(), 'ScalcCvo', table() );

for iS = subGroupToProcess
    
    S(iS).SubGroup = subGroups{iS};
    
    indSubGroup = strcmp(D.SubGroup,S(iS).SubGroup);

    S(iS).NumCases = sum( indSubGroup );
    
    S(iS).CirculationType = unique(D(indSubGroup,:).CirculationType);
       
    % Table (T) and Summary (S) of Measured Values in Absolute Units
    S(iS).Tmeas = M(indSubGroup,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV'});
    S(iS).Smeas = summarize_flows( S(iS).Tmeas, minMeasToSummarize );
    
    % Table (T) and Summary (S) of Derived Values in Absolute Units
    S(iS).Tcalc = D(indSubGroup,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS','CVO'});
    S(iS).Scalc = summarize_flows( S(iS).Tcalc, minMeasToSummarize );
    
    % Table (T) and Summary (S) of Measured Values Relative to CVO
    S(iS).TmeasCvo = S(iS).Tmeas(:,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV'});
    S(iS).TmeasCvo.MPA = 100 * S(iS).Tmeas.MPA ./ S(iS).Tcalc.CVO;
    S(iS).TmeasCvo.AAo = 100 * S(iS).Tmeas.AAo ./ S(iS).Tcalc.CVO;
    S(iS).TmeasCvo.SVC = 100 * S(iS).Tmeas.SVC ./ S(iS).Tcalc.CVO;
    S(iS).TmeasCvo.DA  = 100 * S(iS).Tmeas.DA  ./ S(iS).Tcalc.CVO;
    S(iS).TmeasCvo.DAo = 100 * S(iS).Tmeas.DAo ./ S(iS).Tcalc.CVO;
    S(iS).TmeasCvo.PBF = 100 * S(iS).Tmeas.PBF ./ S(iS).Tcalc.CVO;
    S(iS).TmeasCvo.UV  = 100 * S(iS).Tmeas.UV  ./ S(iS).Tcalc.CVO;
    S(iS).SmeasCvo = summarize_flows( S(iS).TmeasCvo, minMeasToSummarize );
    
    % Table (T) and Summary (S) of Derived Values Relative to CVO
    S(iS).TcalcCvo = S(iS).Tcalc(:,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS'});
    S(iS).TcalcCvo.MPA = 100 * S(iS).Tcalc.MPA ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.AAo = 100 * S(iS).Tcalc.AAo ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.SVC = 100 * S(iS).Tcalc.SVC ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.DA  = 100 * S(iS).Tcalc.DA  ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.DAo = 100 * S(iS).Tcalc.DAo ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.PBF = 100 * S(iS).Tcalc.PBF ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.UV  = 100 * S(iS).Tcalc.UV  ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.IVC = 100 * S(iS).Tcalc.IVC ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.FO  = 100 * S(iS).Tcalc.FO  ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.CA  = 100 * S(iS).Tcalc.CA  ./ S(iS).Tcalc.CVO;
    S(iS).TcalcCvo.CS  = 100 * S(iS).Tcalc.CS  ./ S(iS).Tcalc.CVO;
    S(iS).ScalcCvo = summarize_flows( S(iS).TcalcCvo, minMeasToSummarize );
    
end


%% Summarize Groups

G = struct( 'Group', '', 'NumCases', [], 'CirculationType', cell(0), 'Tmeas', table(), 'Smeas', table(), 'Tcalc', table(), 'Scalc', table(), 'TmeasCvo', table(), 'SmeasCvo', table(), 'TcalcCvo', table(), 'ScalcCvo', table() );

for iG = groupToProcess
    
    G(iG).Group = groups{iG};
    
    indGroup = strcmp(D.Group,G(iG).Group);

    G(iG).NumCases = sum( indGroup );
    
    G(iG).CirculationType = unique(D(indGroup,:).CirculationType);
       
    % Table (T) and Summary (S) of Measured Values in Absolute Units
    G(iG).Tmeas = M(indGroup,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV'});
    G(iG).Smeas = summarize_flows( G(iG).Tmeas, minMeasToSummarize );
    
    % Table (T) and Summary (S) of Derived Values in Absolute Units
    G(iG).Tcalc = D(indGroup,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS','CVO'});
    G(iG).Scalc = summarize_flows( G(iG).Tcalc, minMeasToSummarize );
    
    % Table (T) and Summary (S) of Measured Values Relative to CVO
    G(iG).TmeasCvo = G(iG).Tmeas(:,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV'});
    G(iG).TmeasCvo.MPA = 100 * G(iG).Tmeas.MPA ./ G(iG).Tcalc.CVO;
    G(iG).TmeasCvo.AAo = 100 * G(iG).Tmeas.AAo ./ G(iG).Tcalc.CVO;
    G(iG).TmeasCvo.SVC = 100 * G(iG).Tmeas.SVC ./ G(iG).Tcalc.CVO;
    G(iG).TmeasCvo.DA  = 100 * G(iG).Tmeas.DA  ./ G(iG).Tcalc.CVO;
    G(iG).TmeasCvo.DAo = 100 * G(iG).Tmeas.DAo ./ G(iG).Tcalc.CVO;
    G(iG).TmeasCvo.PBF = 100 * G(iG).Tmeas.PBF ./ G(iG).Tcalc.CVO;
    G(iG).TmeasCvo.UV  = 100 * G(iG).Tmeas.UV  ./ G(iG).Tcalc.CVO;
    G(iG).SmeasCvo = summarize_flows( G(iG).TmeasCvo, minMeasToSummarize );
    
    % Table (T) and Summary (S) of Derived Values Relative to CVO
    G(iG).TcalcCvo = G(iG).Tcalc(:,{'CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS'});
    G(iG).TcalcCvo.MPA = 100 * G(iG).Tcalc.MPA ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.AAo = 100 * G(iG).Tcalc.AAo ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.SVC = 100 * G(iG).Tcalc.SVC ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.DA  = 100 * G(iG).Tcalc.DA  ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.DAo = 100 * G(iG).Tcalc.DAo ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.PBF = 100 * G(iG).Tcalc.PBF ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.UV  = 100 * G(iG).Tcalc.UV  ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.IVC = 100 * G(iG).Tcalc.IVC ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.FO  = 100 * G(iG).Tcalc.FO  ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.CA  = 100 * G(iG).Tcalc.CA  ./ G(iG).Tcalc.CVO;
    G(iG).TcalcCvo.CS  = 100 * G(iG).Tcalc.CS  ./ G(iG).Tcalc.CVO;
    G(iG).ScalcCvo = summarize_flows( G(iG).TcalcCvo, minMeasToSummarize );
    
end


%% Balanced Flow Distribution Models

for iS = subGroupToProcess
    S(iS).Sdist = generate_balanced_flow_distribution_model_using_median_flows( S(iS) );
end

for iG = groupToProcess
    G(iG).Sdist = generate_balanced_flow_distribution_model_using_median_flows( G(iG) );
end


%% Display Sub-Group Results

for iS = subGroupToProcess
    
    fileDesc = sprintf( 'subgroup%02i_%s', iS, strrep( strrep( S(iS).SubGroup, ' ', '_' ), '''', '' ) );
    
    if isLogResults
        fileName = fullfile( outputLogDirPath, sprintf( '%s.txt', fileDesc ) );
        warning off
        delete( fileName );
        warning on
        diary( fileName );
    end

    show_results( S(iS), sprintf( 'Sub-Group %02i: %s', iS, S(iS).SubGroup ), false, minMeasToSummarize )
    
    if isShowFigs
        hFig = show_histogram( S(iS).TmeasCvo, S(iS).SmeasCvo, S(iS).Sdist, strcat( fileDesc, '_measured' ) );
        if isLogResults
            saveas( hFig, fullfile( outputFigDirPath, hFig.Name ), 'png' )
            close( hFig )
        end
        hFig = show_histogram( S(iS).TcalcCvo, S(iS).ScalcCvo, S(iS).Sdist, strcat( fileDesc, '_derived' ) );
        if isLogResults
            saveas( hFig, fullfile( outputFigDirPath, hFig.Name ), 'png' )
            close( hFig )
        end
    end
    
    if isLogResults
        diary off
        clean_log_file( fileName )
    end
    
end


%% Display Group Results

for iG = groupToProcess
    
    fileDesc = sprintf( 'group%02i_%s', iG, strrep( strrep( G(iG).Group, ' ', '_' ), '''', '' ) );
    
    if isLogResults
        fileName = fullfile( outputLogDirPath, sprintf( '%s.txt', fileDesc ) );
        warning off
        delete( fileName );
        warning on
        diary( fileName );
    end

    show_results( G(iG), sprintf( 'Group %02i: %s', iG, G(iG).Group ), false, minMeasToSummarize )
    
    if isShowFigs
        hFig = show_histogram( G(iG).TmeasCvo, G(iG).SmeasCvo, G(iG).Sdist, strcat( fileDesc, '_measured' ) );
        if isLogResults
            saveas( hFig, fullfile( outputFigDirPath, hFig.Name ), 'png' )
            close( hFig )
        end
        hFig = show_histogram( G(iG).TcalcCvo, G(iG).ScalcCvo, G(iG).Sdist, strcat( fileDesc, '_derived' ) );
        if isLogResults
            saveas( hFig, fullfile( outputFigDirPath, hFig.Name ), 'png' )
            close( hFig )
        end
    end
    
    if isLogResults
        diary off
        clean_log_file( fileName )
    end
    
end


%% Display Sub-Group Summary Results

if isLogResults
    fileName = fullfile( outputDirPath, 'summary_subgroup_flow_distribution.txt' );
    warning off
    delete( fileName );
    warning on
    diary( fileName );
end

for iS = subGroupToProcess
    
    show_results( S(iS), sprintf( 'Sub-Group %02i: %s', iS, S(iS).SubGroup ), true, minMeasToSummarize )
    
    if isnan(S(iS).Sdist.MPA(end))
        if S(iS).NumCases < minMeasToSummarize
            fprintf( '\bNote: too few measurements to summarize\n\n\n' )
        else
            fprintf( '\b' )
            generate_balanced_flow_distribution_model_using_median_flows( S(iS) );
            fprintf( '\n\n' )
        end
    end
    
end

if isLogResults
    diary off
    clean_log_file( fileName )
end


%% Display Group Summary Results

if isLogResults
    fileName = fullfile( outputDirPath, 'summary_group_flow_distribution.txt' );
    warning off
    delete( fileName );
    warning on
    diary( fileName );
end

for iG = groupToProcess
    
    show_results( G(iG), sprintf( 'Group %02i: %s', iG, G(iG).Group ), true, minMeasToSummarize )
    
    if isnan(G(iG).Sdist.MPA(end))
        if G(iG).NumCases < minMeasToSummarize
            fprintf( '\bNote: too few measurements to summarize\n\n\n' )
        else
            fprintf( '\b' )
            generate_balanced_flow_distribution_model_using_median_flows( G(iG) );
            fprintf( '\n\n' )
        end
    end
    
end

if isLogResults
    diary off
    clean_log_file( fileName )
end


%% Display Difference Between Measured and Balanced Flows

if isLogResults
    diary( logFilePath );
end

% Groups

SdistQGroup = G(1).Sdist(1,{'MPA','AAo','SVC','DA','DAo','PBF','UV'});
for iG = groupToProcess

    Sdist = G(iG).Sdist([5,6,7],{'MPA','AAo','SVC','DA','DAo','PBF','UV'});
    Sdist{4,:} = Sdist{end,:} - Sdist{end-1,:};
    Sdist.Properties.RowNames{4} = 'difference (%CVO)';
    
    SdistQGroup{iG,:} = Sdist{4,:};
    SdistQGroup.Properties.RowNames{iG} = sprintf( '%02i %s', iG, G(iG).Group );
    
    % descStr = sprintf( 'Group %02i: %s', iG, G(iG).Group );
    % fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
    % display_table( Sdist )
    
end

iG = numel(G)+1;
SdistQGroup{iG,:} = mean( SdistQGroup{1:numel(G),:}, 1, 'omitnan' );
formatStr = sprintf( '%%-%is', max( cellfun( @(x) numel(x), Sdist.Properties.RowNames ) ) );
SdistQGroup.Properties.RowNames{iG} = sprintf( formatStr, 'mean' );
iG = iG+1;
SdistQGroup{iG,:} = median( SdistQGroup{1:numel(G),:}, 1, 'omitnan' );
SdistQGroup.Properties.RowNames{iG} = 'median';
iG = iG+1;
SdistQGroup{iG,:} = max( abs(SdistQGroup{1:numel(G),:}), [], 1 );
SdistQGroup.Properties.RowNames{iG} = 'max abs.';

descStr = sprintf( 'Group Flow Difference, balanced - measured' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );

display_table( SdistQGroup )

% Sub-Groups

SdistQSubGroup = S(1).Sdist(1,{'MPA','AAo','SVC','DA','DAo','PBF','UV'});
for iS = subGroupToProcess

    Sdist = S(iS).Sdist([5,6,7],{'MPA','AAo','SVC','DA','DAo','PBF','UV'});
    Sdist{4,:} = Sdist{end,:} - Sdist{end-1,:};
    Sdist.Properties.RowNames{4} = 'difference (%CVO)';
    
    SdistQSubGroup{iS,:} = Sdist{4,:};
    SdistQSubGroup.Properties.RowNames{iS} = sprintf( '%02i %s', iS, S(iS).SubGroup );
    
    % descStr = sprintf( 'Sub-Group %02i: %s', iS, S(iS).SubGroup );
    % fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
    % display_table( Sdist )
    
end

iS = numel(S)+1;
SdistQSubGroup{iS,:} = mean( SdistQSubGroup{1:numel(S),:}, 1, 'omitnan' );
formatStr = sprintf( '%%-%is', max( cellfun( @(x) numel(x), Sdist.Properties.RowNames ) ) );
SdistQSubGroup.Properties.RowNames{iS} = sprintf( formatStr, 'mean' );
iS = iS+1;
SdistQSubGroup{iS,:} = median( SdistQSubGroup{1:numel(S),:}, 1, 'omitnan' );
SdistQSubGroup.Properties.RowNames{iS} = 'median';
iS = iS+1;
SdistQSubGroup{iS,:} = max( abs(SdistQSubGroup{1:numel(S),:}), [], 1 );
SdistQSubGroup.Properties.RowNames{iS} = 'max abs.';

descStr = sprintf( 'Sub-Group Flow Difference, balanced - measured' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );

display_table( SdistQSubGroup )

if isLogResults
    diary off
    clean_log_file( logFilePath )
end


%% Summarize Results to Excel

% Gather Info to Summarize
iG = 1;
iR = 1;
A = nan(22,51);
desc = cell(numel(S)+numel(G)-1,1);
for iS = 1:numel(S)
    desc{iR} = sprintf( '    %s', S(iS).SubGroup );
    A(iR,1)     = S(iS).NumCases;
    A(iR,2:8)   = S(iS).Smeas{strcmp(S(iS).Smeas.Properties.RowNames,'median'),:};
    A(iR,9:20)  = S(iS).Scalc{strcmp(S(iS).Scalc.Properties.RowNames,'median'),:};
    A(iR,21:27) = S(iS).SmeasCvo{strcmp(S(iS).SmeasCvo.Properties.RowNames,'median'),:}; % 100*S(iS).Smeas{strcmp(S(iS).Smeas.Properties.RowNames,'median'),:}./S(iS).Scalc{strcmp(S(iS).Scalc.Properties.RowNames,'median'),strcmp(S(iS).Scalc.Properties.VariableNames,'CVO')};
    A(iR,28:39) = [ S(iS).ScalcCvo{strcmp(S(iS).Scalc.Properties.RowNames,'median'),:} 100 ]; % 100*S(iS).Scalc{strcmp(S(iS).Scalc.Properties.RowNames,'median'),:}./S(iS).Scalc{strcmp(S(iS).Scalc.Properties.RowNames,'median'),strcmp(S(iS).Scalc.Properties.VariableNames,'CVO')};
    A(iR,40:51) = S(iS).Sdist{strcmp(S(iS).Sdist.Properties.RowNames,'model (%CVO)'),:};
    iR = iR + 1;
    if iS < numel(S)
        if ~strcmp( S(iS).SubGroup(1:3), S(iS+1).SubGroup(1:3) )
            iG = iG + 1;
            desc{iR} = sprintf( 'Group: %s', G(iG).Group );
            A(iR,1)     = G(iG).NumCases;
            A(iR,2:8)   = G(iG).Smeas{strcmp(G(iG).Smeas.Properties.RowNames,'median'),:};
            A(iR,9:20)  = G(iG).Scalc{strcmp(G(iG).Scalc.Properties.RowNames,'median'),:};
            A(iR,21:27) = G(iG).SmeasCvo{strcmp(G(iG).Smeas.Properties.RowNames,'median'),:}; % 100*G(iG).Smeas{strcmp(G(iG).Smeas.Properties.RowNames,'median'),:}./G(iG).Scalc{strcmp(G(iG).Scalc.Properties.RowNames,'median'),strcmp(G(iG).Scalc.Properties.VariableNames,'CVO')};
            A(iR,28:39) = [ G(iG).ScalcCvo{strcmp(G(iG).Scalc.Properties.RowNames,'median'),:} 100 ]; %100*G(iG).Scalc{strcmp(G(iG).Scalc.Properties.RowNames,'median'),:}./G(iG).Scalc{strcmp(G(iG).Scalc.Properties.RowNames,'median'),strcmp(G(iG).Scalc.Properties.VariableNames,'CVO')};
            A(iR,40:51) = G(iG).Sdist{strcmp(G(iG).Sdist.Properties.RowNames,'model (%CVO)'),:};
            iR = iR + 1;
        end
    end
end

% Output Results Command Window / Log
if isLogResults
    diary( logFilePath );
end

descStr = sprintf( 'Summary Flows, absolute units' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
fprintf( '%-32s  %-6s  %-54s  %-54s\n', '', '', 'Median Measured Flow (ml/min/kg)', 'Median Derived Flow (ml/min/kg)' )
fprintf( '%-32s  %-6s   %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s\n', '','Cases','MPA','AAo','SVC','DA','DAo','PBF','UV','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS','CVO')
fprintf( '%-32s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s\n', '','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------' )
for iR = 1:size(A,1)
    fprintf( '%-32s ', desc{iR} )
    for iC = 1:20
        if ~isnan(A(iR,iC))
            fprintf( '%6i  ', round(A(iR,iC)) )
        else
            fprintf( '%6s  ', '' )
        end
    end
    fprintf( '\b\b\n' )
end
fprintf( '\n' )

descStr = sprintf( 'Summary Flows, relative units' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
fprintf( '%-32s  %-6s  %-54s  %-54s\n', '', '', 'Median Measured Flow (%CVO)', 'Median Derived Flow (%CVO)' )
fprintf( '%-32s  %-6s   %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s\n', '','Cases','MPA','AAo','SVC','DA','DAo','PBF','UV','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS','CVO')
fprintf( '%-32s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s\n', '','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------','------' )
for iR = 1:size(A,1)
    fprintf( '%-32s ', desc{iR} )
    for iC = [1,21:39]
        if ~isnan(A(iR,iC))
            fprintf( '%6i  ', round(A(iR,iC)) )
        else
            fprintf( '%6s  ', '' )
        end
    end
    fprintf( '\b\b\n' )
end
fprintf( '\n' )

descStr = sprintf( 'Modelled Flow Distribution' );
fprintf( '\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
fprintf( '%-32s  %-6s  %-54s\n', '', '', 'Modelled Flow (%CVO)' )
fprintf( '%-32s  %-6s   %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s\n', '','Cases','MPA','AAo','SVC','DA','DAo','PBF','UV','ICS','IVC','CA','CS','CVO')
fprintf( '%-32s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s\n', '','------','------','------','------','------','------','------','------','------','------','------','------','------' )
for iR = 1:size(A,1)
    fprintf( '%-32s ', desc{iR} )
    for iC = [1,40:51]
        if ~isnan(A(iR,iC))
            fprintf( '%6i  ', round(A(iR,iC)) )
        else
            fprintf( '%6s  ', '' )
        end
    end
    fprintf( '\b\b\n' )
end
fprintf( '\n' )

if isLogResults
    diary off
end
    
% Save To Excel
resultsTemplateFilePath = fullfile( 'templates', 'results_template.xlsx' ); 
if exist( resultsTemplateFilePath, 'file' )
    resultsFilePath = fullfile( outputDirPath, 'fetal_chd_flows.xlsx' );
    copyfile( resultsTemplateFilePath, resultsFilePath )
    % All Cases
    writetable( M(:,{'MPA','AAo','SVC','DA','DAo','PBF','UV'}), resultsFilePath, 'Sheet', 'Flows', 'Range', 'E2' );
    writetable( D(:,{'MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS','CVO'}), resultsFilePath, 'Sheet', 'Flows', 'Range', 'L2' );
    % Summary
    writematrix( A, resultsFilePath, 'Sheet', 'Summary', 'Range', 'C3' );
end


%% All Measured and Derived Flows

fileName = fullfile( outputLogDirPath, 'all_measured_flows.txt' );
warning off
delete( fileName );
warning on
diary( fileName );
display_table( M(:,{'Group','SubGroup','CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV'}) )
diary off
clean_log_file( fileName )

fileName = fullfile( outputLogDirPath, 'all_derived_flows.txt' );
warning off
delete( fileName );
warning on
diary( fileName );
display_table( D(:,{'Group','SubGroup','CaseNo','MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS','CVO'}) )
diary off
clean_log_file( fileName )


%% Flow Distribution Diagram Information

iG = 1;
iR = 1;
Qmeas = nan(22,8);
Qdist = nan(22,8); 
Qcvo = nan(22,1);
desc = cell(numel(S)+numel(G)-1,1);
for iS = 1:numel(S)
    desc{iR} = sprintf( '    %s', S(iS).SubGroup );
    Qmeas(iR,:) = [ S(iS).SmeasCvo{strcmp(S(iS).Smeas.Properties.RowNames,'median'),:}, NaN ];
    Qdist(iR,:) = S(iS).Sdist{strcmp(S(iS).Sdist.Properties.RowNames,'model (%CVO)'),1:8};
    Qcvo(iR) = S(iS).Scalc.CVO(strcmp(S(iS).Scalc.Properties.RowNames,'median'));
    iR = iR + 1;
    if iS < numel(S)
        if ~strcmp( S(iS).SubGroup(1:3), S(iS+1).SubGroup(1:3) )
            iG = iG + 1;
            desc{iR} = sprintf( 'Group: %s', G(iG).Group );
            Qmeas(iR,:) = [ G(iG).SmeasCvo{strcmp(G(iG).Smeas.Properties.RowNames,'median'),:}, NaN ];
            Qdist(iR,:) = G(iG).Sdist{strcmp(G(iG).Sdist.Properties.RowNames,'model (%CVO)'),1:8};
            Qcvo(iR) = G(iG).Scalc.CVO(strcmp(G(iG).Scalc.Properties.RowNames,'median'));
            iR = iR + 1;
        end
    end
end

if isLogResults
    diary( logFilePath );
end

descStr = sprintf( 'Modelled Flow Distribution Diagram Info' );
fprintf( '\n\n%s\n%s\n\n', descStr, repmat( '=', size(descStr) ) );
fprintf( '%s\n', 'Modelled Flow (%CVO)' )

fprintf( '%-32s %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n', '', 'MPA ', 'AAo ', 'SVC ', 'DA  ', 'DAo ', 'PBF ', 'UV  ', 'FO/ICS' )
fprintf( '%-32s %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n', '', '------', '------', '------', '------', '------', '------', '------', '------' )
for iR = 1:size(Qmeas,1)
    fprintf( '%-32s ', desc{iR} )
    for iC = 1:size(Qmeas,2)
        if round(Qmeas(iR,iC)) == Qdist(iR,iC)
            flowStr = sprintf( ' %3i ', Qdist(iR,iC) );
        elseif isnan(Qmeas(iR,iC)) && ~isnan(Qdist(iR,iC))
            flowStr = sprintf( '[%3i]', Qdist(iR,iC) );
        elseif isnan(Qdist(iR,iC))
            flowStr = sprintf( ' ? ' );
        else
            flowStr = sprintf( '{%3i}', Qdist(iR,iC) );
        end
        fprintf( '%6s  ', flowStr )
    end
    fprintf( '\b\n' )
end
fprintf( '\n\n' )

fprintf( '             ?    flow distribution model undetermined\n' )
fprintf( '            { }   flow distribution model deviates from measured median\n' )
fprintf( '            [ ]   measured median flow unavailable, flow derived by model\n' )
fprintf( '\n' )

fprintf( 'Modelled Flow Distribution Relative to Normal Scaled by Median CVO\n' )
fprintf( '%-32s %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n', '', 'MPA ', 'AAo ', 'SVC ', 'DA  ', 'DAo ', 'PBF ', 'UV  ', 'FO/ICS' )
fprintf( '%-32s %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n', '', '------', '------', '------', '------', '------', '------', '------', '------' )
for iR = 1:size(Qmeas,1)
    fprintf( '%-32s ', desc{iR} )
    r = Qdist(iR,:) ./ Qdist(1,:) * ( Qcvo(iR) / Qcvo(1) );
    if ~all( isnan( r ) )
        fprintf( '%6.2f  ', r )
    end
    fprintf( '\b\n' )
end
fprintf( '\n\n' )

if isLogResults
    diary off
end


%% Save Results Formatted for Manuscript Table

% Define output file path
tableFilePath = fullfile( outputDataDirPath, 'table_manuscript_flow_median_or_range.csv' );

% Create table
varNames = {'type','n','MPA','AAo','SVC','DAo','UV','DA','PBF','FO','CVO'};
varTypes = {'string','singlenan','string','string','string','string','string','string','string','string','string'};
nRow = 1+length(G)+length(S)-1;
nCol = length(varNames);
T = table( 'Size', [nRow,nCol], 'VariableNames', varNames, 'VariableTypes', varTypes );

% Initialize indices
iR = 0; % table row
iG = 0; % group
iS = 1; % subgroup

% Define values used
iR = iR + 1;
T.type(iR) = 'values used';
T.MPA(iR)  = 'meas';
T.AAo(iR)  = 'meas';
T.SVC(iR)  = 'meas';
T.DAo(iR)  = 'meas';
T.UV(iR)   = 'meas';
T.DA(iR)   = 'meas+calc';
T.PBF(iR)  = 'meas+calc';
T.FO(iR)   = 'calc';
T.CVO(iR)  = 'calc';

% Normal group
iR = iR + 1;
iG = iG + 1;
T.type(iR)  = G(iG).Group;
T.n(iR)     = G(iG).NumCases;
T(iR,3:end) = extract_flow_median( G(iG).Smeas, G(iG).Scalc );

% HLHS group
iR = iR + 1;
iG = iG + 1;
T.type(iR)  = G(iG).Group;
T.n(iR)     = G(iG).NumCases;
T(iR,3:end) = extract_flow_median( G(iG).Smeas, G(iG).Scalc );

% HLHS subgroups
% TODO: determine flow summary median or range decision on number of values for vessel instead of 
% number of cases in subgroup; doesn't make a difference for the data here with minMeasToSummarize=3
for iI = 1:6
    iS = iS + 1;
    iR = iR + 1;
    T.type(iR)  = S(iS).SubGroup;
    T.n(iR)     = S(iS).NumCases;
    if T.n(iR) >= minMeasToSummarize  
        T(iR,3:end) = extract_flow_median( S(iS).Smeas, S(iS).Scalc );
    else
        T(iR,3:end) = extract_flow_range( S(iS).Tmeas, S(iS).Tcalc );
    end
end

% TGA group
iR = iR + 1;
iG = iG + 1;
T.type(iR)  = G(iG).Group;
T.n(iR)     = G(iG).NumCases;
T(iR,3:end) = extract_flow_median( G(iG).Smeas, G(iG).Scalc );

% TGA subgroups
for iI = 1:4
    iS = iS + 1;
    iR = iR + 1;
    T.type(iR)  = S(iS).SubGroup;
    T.n(iR)     = S(iS).NumCases;
    if T.n(iR) >= minMeasToSummarize
        T(iR,3:end) = extract_flow_median( S(iS).Smeas, S(iS).Scalc );
    else
        T(iR,3:end) = extract_flow_range( S(iS).Tmeas, S(iS).Tcalc );
    end
end

% TOF group
iR = iR + 1;
iG = iG + 1;
T.type(iR)  = G(iG).Group;
T.n(iR)     = G(iG).NumCases;
T(iR,3:end) = extract_flow_median( G(iG).Smeas, G(iG).Scalc );

% TOF subgroups
for iI = 1:2
    iS = iS + 1;
    iR = iR + 1;
    T.type(iR)  = S(iS).SubGroup;
    T.n(iR)     = S(iS).NumCases;
    if T.n(iR) >= minMeasToSummarize
        T(iR,3:end) = extract_flow_median( S(iS).Smeas, S(iS).Scalc );
    else
        T(iR,3:end) = extract_flow_range( S(iS).Tmeas, S(iS).Tcalc );
    end
end

% EA group
iR = iR + 1;
iG = iG + 1;
T.type(iR)  = G(iG).Group;
T.n(iR)     = G(iG).NumCases;
T(iR,3:end) = extract_flow_median( G(iG).Smeas, G(iG).Scalc );

% EA subgroups
for iI = 1:2
    iS = iS + 1;
    iR = iR + 1;
    T.type(iR)  = S(iS).SubGroup;
    T.n(iR)     = S(iS).NumCases;
    if T.n(iR) >= minMeasToSummarize
        T(iR,3:end) = extract_flow_median( S(iS).Smeas, S(iS).Scalc );
    else
        T(iR,3:end) = extract_flow_range( S(iS).Tmeas, S(iS).Tcalc );
    end
end

% TA group
iR = iR + 1;
iG = iG + 1;
T.type(iR)  = G(iG).Group;
T.n(iR)     = G(iG).NumCases;
T(iR,3:end) = extract_flow_median( G(iG).Smeas, G(iG).Scalc );

% TA subgroups
for iI = 1:2
    iS = iS + 1;
    iR = iR + 1;
    T.type(iR)  = S(iS).SubGroup;
    T.n(iR)     = S(iS).NumCases;
    if T.n(iR) >= minMeasToSummarize
        T(iR,3:end) = extract_flow_median( S(iS).Smeas, S(iS).Scalc );
    else
        T(iR,3:end) = extract_flow_range( S(iS).Tmeas, S(iS).Tcalc );
    end
end

writetable( T, tableFilePath )


%% Save Data Formatted for Prism Import

Q = struct;

% Measured Flows by Sub-Group
vesselNames = {'MPA','AAo','SVC','DA','DAo','PBF','UV'};
subGroupNames = {S.SubGroup};
variableTypes = cell(1,numel(subGroupNames));
for iS = 1:numel(variableTypes)
    variableTypes{iS} = 'double';
end
for iV = 1:numel(vesselNames)
    Q.subgroup.meas.(vesselNames{iV}) = table('Size',[max([S.NumCases]),numel(S)],'VariableNames',subGroupNames,'VariableTypes',variableTypes);
    for iS = 1:numel(subGroupNames)
        q = S(iS).Tmeas.(vesselNames{iV});
        if sum( ~isnan( q ) ) < minMeasToSummarize
            q = nan( size( q ) );
        end
        Q.subgroup.meas.(vesselNames{iV}).(subGroupNames{iS})(1:numel(q)) = q;
        Q.subgroup.meas.(vesselNames{iV}).(subGroupNames{iS})((numel(q)+1):end) = NaN;
    end
    writetable( Q.subgroup.meas.(vesselNames{iV}), fullfile( outputPrismDirPath, sprintf('flow_subgroups_measured_%s.csv',lower(vesselNames{iV})) ) )
end

% Measured Flows by Group
vesselNames = {'MPA','AAo','SVC','DA','DAo','PBF','UV'};
groupNames = {G.Group};
variableTypes = cell(1,numel(groupNames));
for iG = 1:numel(variableTypes)
    variableTypes{iG} = 'double';
end
for iV = 1:numel(vesselNames)
    Q.group.meas.(vesselNames{iV}) = table('Size',[max([G.NumCases]),numel(G)],'VariableNames',groupNames,'VariableTypes',variableTypes);
    for iG = 1:numel(groupNames)
        q = G(iG).Tmeas.(vesselNames{iV});
        if sum( ~isnan( q ) ) < minMeasToSummarize
            q = nan( size( q ) );
        end
        Q.group.meas.(vesselNames{iV}).(groupNames{iG})(1:numel(q)) = q;
        Q.group.meas.(vesselNames{iV}).(groupNames{iG})((numel(q)+1):end) = NaN;
    end
    writetable( Q.group.meas.(vesselNames{iV}), fullfile( outputPrismDirPath, sprintf('flow_groups_measured_%s.csv',lower(vesselNames{iV})) ) )
end

% Derived Flows by Sub-Group
vesselNames = {'MPA','AAo','SVC','DA','DAo','PBF','UV','FO','CVO'};
subGroupNames = {S.SubGroup};
variableTypes = cell(1,numel(subGroupNames));
for iS = 1:numel(variableTypes)
    variableTypes{iS} = 'double';
end
for iV = 1:numel(vesselNames)
    Q.subgroup.calc.(vesselNames{iV}) = table('Size',[max([S.NumCases]),numel(S)],'VariableNames',subGroupNames,'VariableTypes',variableTypes);
    for iS = 1:numel(subGroupNames)
        q = S(iS).Tcalc.(vesselNames{iV});
        if sum( ~isnan( q ) ) < minMeasToSummarize
            q = nan( size( q ) );
        end
        Q.subgroup.calc.(vesselNames{iV}).(subGroupNames{iS})(1:numel(q)) = q;
        Q.subgroup.calc.(vesselNames{iV}).(subGroupNames{iS})((numel(q)+1):end) = NaN;
    end
    writetable( Q.subgroup.calc.(vesselNames{iV}), fullfile( outputPrismDirPath, sprintf('flow_subgroups_derived_%s.csv',lower(vesselNames{iV})) ) )
end

% Derived Flows by Group
vesselNames = {'MPA','AAo','SVC','DA','DAo','PBF','UV','FO','CVO'};
groupNames = {G.Group};
variableTypes = cell(1,numel(groupNames));
for iG = 1:numel(variableTypes)
    variableTypes{iG} = 'double';
end
for iV = 1:numel(vesselNames)
    Q.group.calc.(vesselNames{iV}) = table('Size',[max([G.NumCases]),numel(G)],'VariableNames',groupNames,'VariableTypes',variableTypes);
    for iG = 1:numel(groupNames)
        q = G(iG).Tcalc.(vesselNames{iV});
        if sum( ~isnan( q ) ) < minMeasToSummarize
            q = nan( size( q ) );
        end
        Q.group.calc.(vesselNames{iV}).(groupNames{iG})(1:numel(q)) = q;
        Q.group.calc.(vesselNames{iV}).(groupNames{iG})((numel(q)+1):end) = NaN;
    end
    writetable( Q.group.calc.(vesselNames{iV}), fullfile( outputPrismDirPath, sprintf('flow_groups_derived_%s.csv',lower(vesselNames{iV})) ) )
end


%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract_flow_median
function F = extract_flow_median( Smeas, Scalc )

F = cell2table( { ...
        num2str( round( Smeas(strcmp(Smeas.Properties.RowNames,'median'),strcmp(Smeas.Properties.VariableNames,'MPA')).MPA ) ),...
        num2str( round( Smeas(strcmp(Smeas.Properties.RowNames,'median'),strcmp(Smeas.Properties.VariableNames,'AAo')).AAo ) ),...
        num2str( round( Smeas(strcmp(Smeas.Properties.RowNames,'median'),strcmp(Smeas.Properties.VariableNames,'SVC')).SVC ) ),...
        num2str( round( Smeas(strcmp(Smeas.Properties.RowNames,'median'),strcmp(Smeas.Properties.VariableNames,'DAo')).DAo ) ),...
        num2str( round( Smeas(strcmp(Smeas.Properties.RowNames,'median'),strcmp(Smeas.Properties.VariableNames,'UV')).UV ) ),...
        num2str( round( Scalc(strcmp(Scalc.Properties.RowNames,'median'),strcmp(Scalc.Properties.VariableNames,'DA')).DA ) ),...
        num2str( round( Scalc(strcmp(Scalc.Properties.RowNames,'median'),strcmp(Scalc.Properties.VariableNames,'PBF')).PBF ) ),...
        num2str( round( Scalc(strcmp(Scalc.Properties.RowNames,'median'),strcmp(Scalc.Properties.VariableNames,'FO')).FO ) ),...
        num2str( round( Scalc(strcmp(Scalc.Properties.RowNames,'median'),strcmp(Scalc.Properties.VariableNames,'CVO')).CVO ) )...
    } );

varName = F.Properties.VariableNames;
for iV = 1:numel(varName)
    if isequal(F.(varName{iV}),"NaN")
        F.(varName{iV}) = "N/A";
    end
    % if ismissing(F.(varName{iV}))
    % 	F.(varName{iV}) = "";
    % end    
end

end  % extract_flow_median(...)


%% extract_flow_range
function F = extract_flow_range( Tmeas, Tcalc )

varName = {'MPA','AAo','SVC','DAo','UV','DA','PBF','FO','CVO'};
varType = {'string','string','string','string','string','string','string','string','string'};
F = table('Size',[1,numel(varName)],'VariableNames',varName,'VariableTypes',varType);

for iV = 1:numel(varName)
    
    switch varName{iV}
        case {'MPA','AAo','SVC','DAo','UV'}
            T = Tmeas;
        case {'DA','PBF','FO','CVO'}
            T = Tcalc;
        otherwise
            error( 'variable %s not recognized', varName{iV} )
    end
    
    nValue = sum(~isnan(T.(varName{iV})));
    if nValue == 0      % show emdash if no values
        F.(varName{iV}) = "N/A";
    elseif nValue == 1  % show value if only one value
        F.(varName{iV}) = string( sprintf( '[%i]', round(T.(varName{iV})(~isnan(T.(varName{iV})))) ) );
    elseif nValue > 1   % show range if more than one value
        F.(varName{iV}) = string( sprintf( '[%i,%i]', round(min(T.(varName{iV}))), round(max(T.(varName{iV}))) ) );
    else
        error( 'nValue == %g not recognized', nValue )
    end

end

end  % extract_flow_range(...)


%% derive_flows
function T = derive_flows( T, F, isAA, isPA, isEACS, L, U, isVerbose )

% T - table to store derived flows
% F - table of flows to use in derivation
% L - table of lower outlier thresholds
% U - table of upper outlier thresholds

% NOTE: assumes antegrade flow in DAo, SVC and PBF, as well as AAo in cases without aortic atresia
% (AA) and MPA in cases without pulmonary atresia (PA) or Ebstein's Anomaly with circular shunt (EACS)

if ~exist( 'isVerbose', 'var' )
    isVerbose = false;
end

% IVC
indMissingIVC = isnan(T.IVC);
T.IVC(indMissingIVC) = F.DAo(indMissingIVC);

% Derive Using Arch Junction: Q_AAo - Q_SVC + Q_DA = Q_DAo

deriveStr = 'Derived Using Arch Junction: Q_AAo - Q_SVC + Q_DA = Q_DAo';

% AAo
indMissingAAo = isnan(T.AAo) & ~isAA;
T.AAo(indMissingAAo) = F.DAo(indMissingAAo) - F.DA(indMissingAAo) + F.SVC(indMissingAAo);
if isVerbose && any(indMissingAAo&((T.AAo<L.AAo)|(T.AAo>U.AAo))), fprintf('\n%s\n\n',deriveStr), disp(T(indMissingAAo&((T.AAo<L.AAo)|(T.AAo>U.AAo)),{'SubGroup','CaseNo','AAo'})), end
T.AAo((indMissingAAo)&((T.AAo<=0)|(T.AAo<L.AAo)|(T.AAo>U.AAo))) = NaN;

% TODO: add outlier rejection 

% DAo
indMissingDAo = isnan(T.DAo);
T.DAo(indMissingDAo) = F.AAo(indMissingDAo) - F.SVC(indMissingDAo) + F.DA(indMissingDAo);
if isVerbose && any(indMissingDAo&((T.DAo<L.DAo)|(T.DAo>U.DAo))), fprintf('\n%s\n\n',deriveStr), disp(T(indMissingDAo&((T.DAo<L.DAo)|(T.DAo>U.DAo)),{'SubGroup','CaseNo','DAo'})), end
T.DAo((indMissingDAo)&((T.DAo<=0)|(T.DAo<L.DAo)|(T.DAo>U.DAo))) = NaN;

% SVC
indMissingSVC = isnan(T.SVC);
T.SVC(indMissingSVC) = F.AAo(indMissingSVC) + F.DA(indMissingSVC) - F.DAo(indMissingSVC);
if isVerbose && any(indMissingSVC&((T.SVC<L.SVC)|(T.SVC>U.SVC))), fprintf('\n%s\n\n',deriveStr), disp(T(indMissingSVC&((T.SVC<L.SVC)|(T.SVC>U.SVC)),{'SubGroup','CaseNo','SVC'})), end
T.SVC((indMissingSVC)&((T.SVC<=0)|(T.SVC<L.SVC)|(T.SVC>U.SVC))) = NaN;

% DA
indMissingDA = isnan(T.DA);
T.DA(indMissingDA) = F.DAo(indMissingDA) - F.AAo(indMissingDA) + F.SVC(indMissingDA);
if isVerbose && any(indMissingDA&((T.DA<L.DA)|(T.DA>U.DA))), fprintf('\n%s\n\n',deriveStr), disp(T(indMissingDA&((T.DA<L.DA)|(T.DA>U.DA)),{'SubGroup','CaseNo','DA'})), end
T.DA((indMissingDA)&((T.DA<L.DA)|(T.DA>U.DA))) = NaN;


% Derive Using Pulmonary Junction: Q_MPA = Q_DA + Q_PBF

deriveStr = 'Derived Using Pulmonary Junction: Q_MPA = Q_DA + Q_PBF';

% MPA
indMissingMPA = isnan(T.MPA) & ~isPA;
T.MPA(indMissingMPA) = F.DA(indMissingMPA) + F.PBF(indMissingMPA);
if isVerbose && any(indMissingMPA&((T.MPA<L.MPA)|(T.MPA>U.MPA))), fprintf('\n%s\n\n',deriveStr), disp(T(indMissingMPA&((T.MPA<L.MPA)|(T.MPA>U.MPA)),{'SubGroup','CaseNo','MPA'})), end
T.MPA((indMissingMPA)&(~isEACS)&((T.MPA<=0)|(T.MPA<L.MPA)|(T.MPA>U.MPA))) = NaN;

% DA
indMissingDA = isnan(T.DA);
T.DA(indMissingDA) = F.MPA(indMissingDA) - F.PBF(indMissingDA);
if isVerbose && any(indMissingDA&((T.DA<L.DA)|(T.DA>U.DA))), fprintf('\n%s\n\n',deriveStr), disp(T(indMissingDA&((T.DA<L.DA)|(T.DA>U.DA)),{'SubGroup','CaseNo','DA'})), end
T.DA((indMissingDA)&((T.DA<L.DA)|(T.DA>U.DA))) = NaN;

% PBF
indMissingPBF = isnan(T.PBF);
T.PBF(indMissingPBF) = F.MPA(indMissingPBF) - F.DA(indMissingPBF);
if isVerbose && any(indMissingPBF&((T.PBF<L.PBF)|(T.PBF>U.PBF))), fprintf('\n%s\n\n',deriveStr), disp(T(indMissingPBF&((T.PBF<L.PBF)|(T.PBF>U.PBF)),{'SubGroup','CaseNo','PBF'})), end
T.PBF((indMissingPBF)&((T.PBF<=0)|(T.PBF<L.PBF)|(T.PBF>U.PBF))) = NaN;

end  % derive_flows(...)


%% summarize_flows
function S = summarize_flows( T, minMeasToSummarize )

    iMpa = find(strcmp(T.Properties.VariableNames,'MPA'));
    variableNames = T.Properties.VariableNames(iMpa:end);
    rowNames      = {'N','mean','standard deviation','median','median absolute deviation','16th - 50th pct.','84th - 50th pct.'};
    
    nVessel       = numel(variableNames);
    nRow          = numel(rowNames);
    
    variableTypes = cell(1,nVessel);
    for iV = 1:nVessel
        variableTypes{iV} = 'double';
    end

    S = table( 'Size',[nRow,nVessel],'VariableTypes',variableTypes,'VariableNames',variableNames,'RowNames',rowNames );

    for iV = 1:nVessel
        vesselName = S.Properties.VariableNames{iV};
        S.(vesselName)(1) = sum( ~isnan( T.(vesselName) ) );
        if sum(~isnan(T.(vesselName))) >= minMeasToSummarize 
            for iR = 2:nRow
                switch S.Row{iR}
                    case 'mean'
                        S.(vesselName)(iR) = mean( T.(vesselName), 'omitnan' );
                    case 'standard deviation'
                        S.(vesselName)(iR) = std( T.(vesselName), 'omitnan' );
                    case 'median'
                        S.(vesselName)(iR) = median( T.(vesselName), 'omitnan' );
                    case 'median absolute deviation'
                        S.(vesselName)(iR) = mad( T.(vesselName), 1 );
                    case '16th - 50th pct.'
                        S.(vesselName)(iR) = median( T.(vesselName), 'omitnan' ) - prctile( T.(vesselName), 16 );
                    case '84th - 50th pct.'
                        S.(vesselName)(iR) = prctile( T.(vesselName), 84 ) - median( T.(vesselName), 'omitnan' );
                end
            end
        else
            for iR = 2:size(S,1)
                S.(vesselName)(iR) = NaN;
            end 
        end
    end

end  % summarize_flows(...)


%% generate_balanced_flow_distribution_model_using_median_flows
function Sdist = generate_balanced_flow_distribution_model_using_median_flows( R )

    iR = [ find( strcmp( R.Smeas.Row, 'N' ) ), find( strcmp( R.Smeas.Row, 'median' ) ) ];
    Sdist = R.Smeas(iR,{'MPA','AAo','SVC','DA','DAo','PBF','UV'});
    Sdist.Properties.RowNames{1} = 'measured, N (ml/min/kg)';
    Sdist.Properties.RowNames{2} = 'measured median (ml/min/kg)';    
    Sdist.FO  = [0;NaN];
    Sdist.IVC = [0;NaN];
    Sdist.CA  = [0;NaN];
    Sdist.CS  = [0;NaN];
    Sdist.CVO = [0;NaN];

    
    iR = [ find( strcmp( R.Scalc.Row, 'N' ) ), find( strcmp( R.Scalc.Row, 'median' ) ) ];
    Sdist = [ Sdist; R.Scalc(iR,{'MPA','AAo','SVC','DA','DAo','PBF','UV','FO','IVC','CA','CS','CVO'}) ];
    Sdist.Properties.RowNames{3} = 'derived, N (ml/min/kg)';
    Sdist.Properties.RowNames{4} = 'derived median (ml/min/kg)';
    
    iR = [ find( strcmp( R.SmeasCvo.Row, 'N' ) ), find( strcmp( R.SmeasCvo.Row, 'median' ) ) ];
    Stemp = R.SmeasCvo(iR,{'MPA','AAo','SVC','DA','DAo','PBF','UV'});
    Stemp.Properties.RowNames{1} = 'measured, N (%CVO)';
    Stemp.Properties.RowNames{2} = 'measured median (%CVO)';
    Stemp.FO  = [0;NaN];
    Stemp.IVC = [0;NaN];
    Stemp.CA  = [0;NaN];
    Stemp.CS  = [0;NaN];
    Stemp.CVO = [0;NaN];
    Sdist = [ Sdist; Stemp ];
    
    Q = balance_flow_distribution( Sdist, R );
    CVO = max(Q(1),0) + max(Q(2),0) + Q(10)*(Q(2)>0);
    Sdist = [ Sdist; num2cell([Q,CVO]) ];
    Sdist.Properties.RowNames{end} = 'model (%CVO)'; 
    
    iV = strcmp( Sdist.Properties.VariableNames, 'FO' );
    Sdist.Properties.VariableNames{iV} = 'ICS';
    
end  % generate_balanced_flow_distribution_model_using_median_flows(...)


%% balance_flow_distribution
function Q = balance_flow_distribution( Sdist, R )
    
    maxIter = 10;
    
    q = round(Sdist{end,1:11}); 
    if isfield( R, 'SubGroup' )
        switch R.SubGroup
            case {'HLHS RAS','HLHS MS AA','HLHS MA AA'}
                q(1,2) = -3; 
            case {'TOF PA','Ebstein''s no Circular Shunt'}
                q(1,1) = 0;
        end
    end

    if numel(R.CirculationType) > 1
        Q = nan( size(q) );
        warning( 'cannot balance flows using central tendancy (e.g., mean) of flows from multiple, different circulation types' )
        return
    end
    
    if all(isnan(q))
        Q = nan( size(q) );
        warning( 'cannot balance flow distribution if all flows are NaN' )
        return
    end
    
    iQ = 2;
    [q(iQ,:),~,~,exitFlag] = balance_fetal_flow_distribution( round(q(iQ-1,:)), R.CirculationType, 'objFn', 'wRMSD', 'isPctCvo', true, 'isVerbose', false );
    if exitFlag<0
        Q = q(iQ,:);
        return
    end

    while ~isequal( q(iQ,:), q(iQ-1,:) ) && exitFlag>=0
        iQ = iQ + 1;
        [q(iQ,:),~,~,exitFlag] = balance_fetal_flow_distribution( round(q(iQ-1,:)), R.CirculationType, 'objFn', 'wRMSD', 'isPctCvo', true, 'isVerbose', false );
        if iQ >= maxIter
            Q = q(iQ,:);
            warning( 'maximum iterations, %i, reached', maxIter )
            return
        end
    end

    iQ = iQ + 1;
    [q(iQ,:),~,~,exitFlag] = balance_fetal_flow_distribution( round(q(iQ-1,:)), R.CirculationType, 'objFn', 'wMAD', 'isPctCvo', true, 'isVerbose', false );

    while ~isequal( q(iQ,:), q(iQ-1,:) ) && exitFlag>=0
        iQ = iQ + 1;
        [q(iQ,:),~,~,exitFlag] = balance_fetal_flow_distribution( round(q(iQ-1,:)), R.CirculationType, 'objFn', 'wMAD', 'isPctCvo', true, 'isVerbose', false );
        if iQ >= maxIter
            Q = q(iQ,:);
            warning( 'maximum iterations, %i, reached', maxIter )
            return
        end
    end
    
    Q = q(iQ,:);

end  % balance_flow_distribution(...)


%% display_table
function display_table( T, firstColWidth )

R = T(:,find(strcmp(T.Properties.VariableNames,'MPA')):end);

R = varfun( @(x) round(x), R );

T(:,find(strcmp(T.Properties.VariableNames,'MPA')):end) = R;

minColWidth = 4;
variableNameStr = sprintf( '%%-%is', minColWidth );

for iV = 1:numel(T.Properties.VariableNames)
    T.Properties.VariableNames{iV} = sprintf( variableNameStr, T.Properties.VariableNames{iV} );
end

if exist( 'firstColWidth', 'var' )
    variableNameStr = sprintf( '%%-%is', firstColWidth );
    T.Properties.VariableNames{1} = sprintf( variableNameStr, T.Properties.VariableNames{1} );
end

disp(T)

end  % display_table(...)


%% show_results
function show_results( R, descStr, isShowDistributionOnly, minMeasToSummarize )

    fprintf( '%s\n%s\n\n\n', descStr, repmat( '=', size(descStr) ) )
    
    fprintf( 'Flow Distribution Model\n\n' )
    display_table( R.Sdist )
    fprintf( '\n' )
    
    if ~exist( 'isShowDistributionOnly', 'var' )
        isShowDistributionOnly = false;
    end
    
    if ~exist( 'minMeasToSummarize', 'var' )
        minMeasToSummarize = 4;
    end
    
    if ~isShowDistributionOnly
        
        fprintf( '%s\n%s\n\n\n', 'Summary', repmat( '-', size('Summary') ) )
        
        if R.NumCases >= minMeasToSummarize
        
            fprintf( 'Measured (ml/min/kg)\n\n' )
            display_table( R.Smeas )
            fprintf( '\n' )
            
            fprintf( 'Derived (ml/min/kg)\n\n' )
            display_table( R.Scalc )
            fprintf( '\n' )

            fprintf( 'Measured (%%CVO)\n\n' )
            display_table( R.SmeasCvo )
            fprintf( '\n' )
            
            fprintf( 'Derived (%%CVO)\n\n' )
            display_table( R.ScalcCvo )
            fprintf( '\n' )
        
        else
            
            fprintf( 'Number of cases (%i) < minimum for summary (%i)\n\n\n', R.NumCases, minMeasToSummarize )
        
        end
        
        fprintf( '%s\n%s\n\n\n', 'Flows', repmat( '-', size('Flows') ) )
        
        fprintf( 'Measured (ml/min/kg)\n\n' )
        display_table( R.Tmeas )
        fprintf( '\n' )
        
        fprintf( 'Derived (ml/min/kg)\n\n' )
        display_table( R.Tcalc )
        fprintf( '\n' )
        
        fprintf( 'Measured (%%CVO)\n\n' )
        display_table( R.TmeasCvo )
        fprintf( '\n' )
        
        fprintf( 'Derived (%%CVO)\n\n' )
        display_table( R.TcalcCvo )
        fprintf( '\n' )
        
    end

end  % show_results(...)


%% show_histogram
function hFig = show_histogram( T, S, B, descStr )

hFig = figure( 'Name', descStr );

hFig.Position(3:4) = 2*hFig.Position(3:4);

maxFrq  = 2.5*sqrt(size(T,1));

iMpa = find(strcmp(T.Properties.VariableNames,'MPA'));

subplot(3,3,1)
hAx = gca;
hAx.XLim = [0,1];
hAx.YLim = [0,1];
hAx.XTick = [];
hAx.YTick = [];
title( strrep(descStr,'_',' ') )
text(0.4,0.8, sprintf( '%i cases', size(T,1) ) )
axis off

switch descStr((end-7):end)
    case 'measured' 
        faceColor = '#0072BD';
        legendStr = 'measured';
    case '_derived'
        faceColor = '#EDB120';
        legendStr = 'derived';
end

iP = 2;
for iV = iMpa+(0:7)
    if iV <= size(T,2) && iP <= 9
        subplot(3,3,iP)
        histogram( T.(T.Properties.VariableNames{iV}), 'FaceColor', faceColor );
        title( strrep(T.Properties.VariableNames(iV),'FO','FO/ICS') )
        xlabel('Flow (%CVO)')
        if ismember(iP,[2,4,7])
            ylabel('Frequency')
        end
        set(gca,'XLim',105*[-1,+1],'YLim',[0,maxFrq])
        statStrCell = cell(0);
        statStrCell{1} = sprintf( '%-6s = %3i\n', 'N', round(S.(T.Properties.VariableNames{iV})(1)) );
        statStrCell{2} = sprintf( '%-6s %s %s\n   %3i %s %2i', 'mean', char(177), 'std.', round(S.(T.Properties.VariableNames{iV})(2)), char(177), round(S.(T.Properties.VariableNames{iV})(3)) );
        statStrCell{3} = sprintf( '%-6s %s %s\n   %3i %s %2i', 'median', char(177), 'm.a.d.', round(S.(T.Properties.VariableNames{iV})(4)), char(177), round(S.(T.Properties.VariableNames{iV})(5)) );
        statStrCell{4} = sprintf( '%-6s\n   %3i', 'model', B.(strrep(T.Properties.VariableNames{iV},'FO','ICS'))(end) );
        statsStr = sprintf( '%s\n%s\n%s\n%s', statStrCell{1}, statStrCell{2}, statStrCell{3}, statStrCell{4} );
        text( -100, 0.6*maxFrq, statsStr, 'FontSize', 10, 'FontName', 'FixedWidth' )
        legend(legendStr,'Location','SW')
        iP = iP + 1;
    end
end
            
end  % show_histogram(...)


%% clean_log_file
function clean_log_file( fileName )

fid = fopen(fileName,'r');
f = fread(fid,'*char')';
fclose(fid);
f = strrep( strrep( f, '<strong>', '' ), '</strong>', '');
f = regexprep(f,'In <a href="matlab:matlab.internal.language.introspective.errorDocCallback\(.+line \d+</a>\)','','dotexceptnewline');
fid = fopen(fileName,'w');
fprintf(fid,'%s',f);
fclose(fid);

end  % clean_log_file(...)


%% balance_fetal_flow_distribution
function [ optFlow, unoptFlow, wgt, exitFlag ] = balance_fetal_flow_distribution( measFlow, circType, NVArg )
%BALANCE_FETAL_FLOW_DISTRIBUTION  optimize fetal flow distribution so that flow is conserved 
%
%   optFlow  = BALANCE_FETAL_FLOW_DISTRIBUTION( measFlow ) returns a set flow values that fit the 
%   constraints of normal fetal ciruclation as part of a nonlinear optimization. Any flows that were 
%   not measured should be specified as NaN in measFlow and will be derived instead.
%
%   measFlow = [ mpa aao svc da dao pbf uv ics ivc ca cs ] measured flows in ml/min/kg fetal mass
%   optFlow  = [ mpa aao svc da dao pbf uv ics ivc ca cs ] after optimization
%
%   [ ..., unoptFlow, wgt ] = BALANCE_FETAL_FLOW_DISTRIBUTION( ... ) additionally returns flow initial
%   values and weights used in optimization.
%
%   BALANCE_FETAL_FLOW_DISTRIBUTION( ..., circType ) specifies type of circulation; default 'normal'.
%
%   BALANCE_FETAL_FLOW_DISTRIBUTION( ..., 'Name', 'Value' ) specifies addtional name-value arguments.
%
%   vessels     MPA     main pulmonary artery
%               AAo     ascending aorta
%               SVC     superior vena cava
%               DA      ductus arteriosus
%               DAo     descending aorta
%               PBF     pulmonary blood flow
%               UV      umbilical vein
%               ICS     intra-cardiac shunt, from right heart to left
%               IVC     inferior vena cava
%               CA      coronary artery
%               CS      coronary vein
%
%   The optimization seeks to minimize an objective function, f(x), in this case, the weighted mean 
%   absolute difference between initial flow values (x0=unoptFlow) and balanced flows (x=optFlow).
%   i.e., 
%           f(x) = sum( w .*  abs(x-x0) ) / sum( w );
%
%   Optional Arguments
%   circType =   character string specifying type of circulation, valid options are 'normal', 
%                'transposition', 'ebsteins_circular_shunt', 'double_outlet_right_ventricle' and
%                'hlhs_aortic_atresia'
% 
%   Name-Value Arguments
%   'wgt'       vector of weights corresponding to realiability of values in measFlow; if not
%               specified values typical of circlation will be used, subject to flows measured
%   'caFlowPerCvo' the fraction of combined ventricular output flowing through the coronary artery
%   'csFlowPerCai' the fraction of combined atrial input flowing through the coronary sinus
%   'isVerbose' display results
%
%   Example
%       %          [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
%       measFlow = [                                           NaN   NaN   NaN   NaN ];
%       optFlow  = balance_fetal_flow_distn( measFlow, 'normal' );  
%       % returns  [                                           NaN   NaN   NaN   NaN ] as optFlow
%
%   Note
%       * ICS represents net flow through foramen ovale flow and any defects (e.g., ventricular
%         septal defect) present.
%       * Extension on original implementation described in Prsa2014.
%           Prsa M, Sun L, van Amerom J, Yoo S-J, Grosse-Wortmann L, Jaeggi E, Macgowan C, Seed M. 
%           Reference Ranges of Blood Flow in the Major Vessels of the Normal Human Fetal Circulation 
%           at Term by Phase Contrast Magnetic Resonance Imaging. 
%           Circ. Cardiovasc. Imaging 2014:663-671. doi: 10.1161/CIRCIMAGING.113.001859.
%       * The fmincon optimization algorithm sequential quadratic programming algorithm is used, not 
%         the active set algorithm as described in Prsa2014.
%
%   See also fmincon.

% NOTE: requires Matlab R2019b or later due to use of arguments code block for argument validation

% TODO: add note to help text re. Qref for unique cirulations

%   jfpva (joshua.vanamerom@sickkids.ca)


%% Argument Validation

arguments  

  % name                size        class       validation_functions                default_value   units / notes
  % -----------------   ---------   ---------   ---------------------------------   -------------   -------------

    % Positional Arguments

    measFlow            (1,11)      double                                                          % ml/min/kg
    circType                        char                                            = 'normal'

    % Name-Value Arguments

    NVArg.measWgt       (1,11)      double      {mustBeNonnegative}                 = inf(1,11)    
    NVArg.caFlowPerCvo  (1,1)       double      {mustBeNonnegative}                 = 0.03          % Q_ca as a percent of CVO (combined ventricular output)
    NVArg.csFlowPerCai  (1,1)       double      {mustBeNonnegative}                 = 0.03          % Q_cs as a percent of CAI (combined atrial input)
    NVArg.Qref          (1,1)       struct                                          = struct()
    NVArg.objFn                     char                                            = 'wRMSD'
    NVArg.isPctCvo      (1,1)       logical     {mustBeNonempty}                    = false         % All flow are expressed as percent of CVO (combined ventricular output)
    NVArg.isVerbose     (1,1)       logical     {mustBeNonempty}                    = true 

end

% Circulation Type
validCircType = {'normal','transposition','hlhs_aortic_atresia','tof_pa',...
                 'double_outlet_right_ventricle','ebsteins_circular_shunt'};
assert( ismember(circType,validCircType), 'circType ''%s'' invalid', circType )

% Measurement Weights
measWgt = NVArg.measWgt;

% Coronary Flows
assert( NVArg.caFlowPerCvo <= 1, 'caFlowPerCvo cannot be > 1' )
assert( NVArg.csFlowPerCai <= 1, 'csFlowPerCai cannot be > 1' )
caFlowPerCvo = NVArg.caFlowPerCvo;
csFlowPerCai = NVArg.csFlowPerCai;

% Objective Function Type
validObjFn = {'wMAD','wRMSD'};
assert( ismember(NVArg.objFn,validObjFn), 'objFn ''%s'' invalid', NVArg.objFn )


% Verbose Output
isVerbose = NVArg.isVerbose;


%% Extract Flows and Weights

vesselDesc = {'mpa','aao','svc', 'da','dao','pbf', 'uv','ics','ivc','ca','cs'};
Q = struct;
W = struct;
for iV = 1:numel(vesselDesc)
    eval( sprintf( 'Q.%s = measFlow(%i);', vesselDesc{iV}, iV ) );
    eval( sprintf( 'W.%s = measWgt(%i);', vesselDesc{iV}, iV ) );
    if isnan( Q.(vesselDesc{iV}) )  % set weight to zero if flow not measured
        W.(vesselDesc{iV}) = 0;
    end
    if ( Q.(vesselDesc{iV}) == 0)   % set weight to zero if flow was assumed to be zero
        W.(vesselDesc{iV}) = 0;
    end
end


%% Flow Validation

switch circType
    case 'ebsteins_circular_shunt'
        if isnan( Q.mpa ) && isnan( Q.pbf )
            warning( 'circulation type %s requires Q_mpa and/or Q_pbf measurements to balance flow distribution, no optimization performed', circType )
            optFlow     = nan( size( measFlow ) );
            unoptFlow   = nan( size( measFlow ) );
            wgt         = nan( size( measFlow ) );
            exitFlag    = -5;
            return
        end
end


%% Reference Flows

if all(isfield(NVArg.Qref,['cvo',vesselDesc]))
    Qref = NVArg.Qref;
else
    switch circType
        case {'normal','double_outlet_right_ventricle'}
            % Based on values for healthy fetsuses published in Prsa2014(CircCardiovascularImaging).
            Qref.cvo = 465;  
            Qref.mpa = 261;
            Qref.aao = 191;
            Qref.svc = 137;
            Qref.da  = 187;
            Qref.dao = 252;
            Qref.pbf =  74;
            Qref.uv  = 135;
            Qref.ics = 134;
            Qref.ivc = Qref.dao;
            Qref.ca  = 0.03 * Qref.cvo;
            Qref.cs  = Qref.ca;
        case 'transposition'  
            % Based on values published in Porayette2015(ISMRM).
            Qref.cvo = 482;  
            Qref.mpa = 211;
            Qref.aao = 272;
            Qref.svc = 170;
            Qref.da  = 78;
            Qref.dao = 250;
            Qref.pbf =  83;
            Qref.uv  = 133;
            Qref.ics = NaN;
            Qref.ivc = Qref.dao;
            Qref.ca  = 0.03 * Qref.cvo;
            Qref.cs  = Qref.ca;
        case 'hlhs_aortic_atresia'
             % Based on values published in AlNafisi2013(JCMR).
            Qref.cvo = 456;  
            Qref.mpa = 456;
            Qref.aao = -14;
            Qref.svc = 132;
            Qref.da  = 358;
            Qref.dao = 246;
            Qref.pbf =  81;
            Qref.uv  = 133;
            Qref.ics = -81;
            Qref.ivc = Qref.dao;
            Qref.ca  = 0.03 * Qref.cvo;
            Qref.cs  = Qref.ca;
        case 'ebsteins_circular_shunt'
            % Based on mean measured values in 5 cases as well as Q_PBF from healthy fetsuses published in Prsa2014(CircCardiovascularImaging).
            Qref.cvo = 423/.97;  
            Qref.mpa = -161;
            Qref.aao = 423;
            Qref.svc =  56;
            Qref.da  = -235;
            Qref.dao = 132;
            Qref.pbf =  74;
            Qref.uv  =  90;
            Qref.ics = 1.0309 * Qref.aao -Qref.pbf;
            Qref.ivc = Qref.svc;
            Qref.ca  = 0.03 * Qref.cvo;
            Qref.cs  = Qref.ca;
        case 'tof_pa'
            % Based on mean measured values in 4 cases.
            Qref.cvo = 407/.97;  
            Qref.mpa =   0;
            Qref.aao = 407;
            Qref.svc = 158;
            Qref.da  = -26;
            Qref.dao = 223;
            Qref.pbf =  26;
            Qref.uv  = 101;
            Qref.ivc = Qref.dao;
            Qref.ca  = 0.03 * Qref.cvo;
            Qref.cs  = Qref.ca;
            Qref.ics = Qref.svc + Qref.dao + Qref.cs;
        otherwise
            error( 'circulation type %s not recognized, no reference flows defined', circType )
    end
end

if NVArg.isPctCvo
	vesselField = fieldnames( Qref );
    for iV = 1:numel(vesselField)
        if ~strcmp( vesselField{iV}, 'cvo' )
            Qref.(vesselField{iV}) = 100 * Qref.(vesselField{iV}) / Qref.cvo;
        end
    end
end


%% Assign Missing Flows

% Use Reference Values for Major Vessels
if isnan(Q.mpa)
    Q.mpa = Qref.mpa;
end
if isnan(Q.aao)
    Q.aao = Qref.aao;
end
if isnan(Q.svc)
    Q.svc = Qref.svc;
end
if isnan(Q.dao)
    Q.dao = Qref.dao;
end
if isnan(Q.uv)
    Q.uv = Qref.uv;
end

% Derive IVC Flow
if isnan(Q.ivc)
    Q.ivc = Q.dao;
end

% Derive DA Flow
if isnan(Q.da)
    Q.da = Q.dao - Q.aao + Q.svc;
end

% Derive PBF
if isnan(Q.pbf)
    Q.pbf = Q.mpa - Q.da;
end

% Derive Coronary Flow
if isnan(Q.ca)
    switch circType
        case {'normal','transposition','double_outlet_right_ventricle'}
            Q.ca  = caFlowPerCvo/(1-caFlowPerCvo)*(Q.mpa+Q.aao);
        case {'ebsteins_circular_shunt','tof_pa'}
            Q.ca  = caFlowPerCvo/(1-caFlowPerCvo)*(Q.aao);
        case 'hlhs_aortic_atresia'
            Q.ca  = caFlowPerCvo*Q.mpa;
        otherwise 
            error( 'circulation type %s not recognized, coronary artery flow cannot be initialized', circType )
    end
end
if isnan(Q.cs)
    switch circType
        case {'normal','transposition','double_outlet_right_ventricle','hlhs_aortic_atresia'}
            Q.cs  = csFlowPerCai/(1-csFlowPerCai)*(Q.svc+Q.ivc+Q.pbf);
        case {'ebsteins_circular_shunt','tof_pa'}
            Q.cs  = csFlowPerCai/(1-csFlowPerCai)*(Q.svc+Q.ivc+Q.pbf-Q.mpa);
        otherwise 
            error( 'circulation type %s not recognized, coronary artery flow cannot be initialized', circType )
    end
end

% Derived Intra-cardiac Shunt Flow
if isnan(Q.ics)
    switch circType
        case {'normal','ebsteins_circular_shunt','tof_pa'}
            Q.ics = Q.aao + Q.ca - Q.pbf; 
        case 'transposition'
            Q.ics = Q.mpa - Q.pbf;
        case {'double_outlet_right_ventricle','hlhs_aortic_atresia'}
            Q.ics = Q.ca - Q.pbf; 
        otherwise
            error( 'circulation type %s not recognized, intra-cardiac shunt flow cannot be initialized', circType )
    end
end


%% Assign Measurement Weights

switch circType
    case {'normal', 'transposition','hlhs_aortic_atresia'}
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        refWgt = [ 0.6   0.6   0.8   0.5   1.0   0.28  1.0     0     0     0     0 ];  % per MS
    case 'ebsteins_circular_shunt'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        refWgt = [ 0.3   0.6   0.8   0.5   1.0   0.28  1.0     0     0     0     0 ];  
    case 'double_outlet_right_ventricle'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        refWgt = [ 0.6   0.3   0.8   0.5   1.0   0.28  1.0     0     0     0     0 ];  
    case 'tof_pa'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        refWgt = [ 0.0   0.6   0.8   0.5   1.0   0.28  1.0     0     0     0     0 ];  
    otherwise
        error( 'circulation type %s not recognized, intra-cardiac shunt flow cannot be initialized', circType )
end


%% Setup Optimization

% minimize objective function (f) subject to linear inequality constraints (A*x<=b) and linear
% equality constraints (Aeq*x=b) starting with initial values (x0) within lower (lb) and upper 
% bounds (ub)

% optimization values
xMeas  = [ Q.mpa, Q.aao, Q.svc, Q.da, Q.dao, Q.pbf, Q.uv, Q.ics, Q.ivc, Q.ca, Q.cs ];
wMeas  = [ W.mpa, W.aao, W.svc, W.da, W.dao, W.pbf, W.uv, W.ics, W.ivc, W.ca, W.cs ];
wMeas( isinf( wMeas ) ) = refWgt( isinf( wMeas ) );

% initial values
x0     = xMeas;

% objective function
switch NVArg.objFn
    case 'wMAD'
        f = @(x) sum( wMeas .*  abs(x-x0) ) / sum( wMeas );  % weighted mean absolute difference
    case 'wRMSD'
        f = @(x) sqrt( sum( wMeas .*  (x-x0).^2 ) / sum( wMeas ) );  % weighted root mean squared difference
        %{ WIP: penalize small flow values
        %{
        xMinPct = 5;
        if NVArg.isPctCvo
            xMin = xMinPct;
        else
            xMin = XminPct * Qref.cvo;
        end
        f = @(x) sqrt( sum( wMeas .*  (x-x0).^2 ) / sum( wMeas ) ) + sqrt( mean( (abs(x0)>xMin) .* ((max(0,(xMin-abs(x)))).^2) ) );
        %}
    case 'wRMSD'
        f = @(x) sqrt( sum( wMeas .*  (x-x0).^2 ) / sum( wMeas ) );
    otherwise
        error( 'objective function type %s not recognized', NVArg.ojbFn )
end

% linear inequality constraint: A*x <= b
%        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
A      = [   0     0     0     0    -1     0     1     0     0     0     0 ]; % placenta
b      = [ 0 ];

% nonlinear constraints
nonlcon = [];

% options
options = optimoptions('fmincon','Display','none','Algorithm','sqp');  % NOTE: to show intermediate results use Display 'iter-detailed'


%% Setup Optimization Equality Constraints

% linear equality constraint: Aeq*x = beq

Rca    = (1-caFlowPerCvo)/caFlowPerCvo;  % ratio of non-coronoary ventricular output to coronary artery flow
Rcs    = (1-csFlowPerCai)/csFlowPerCai;  % ratio of non-coronoary atrial input to coronary sinus flow

%                                    [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
EqNormal.heart      = struct( 'Aeq', [  -1    -1     1     0     0     1     0     0     1    -1     1 ], 'beq', 0 );  % heart chambers
EqNormal.pa         = struct( 'Aeq', [   1     0     0    -1     0    -1     0     0     0     0     0 ], 'beq', 0 );  % pulmonary artery branch junction
EqNormal.arch       = struct( 'Aeq', [   0     1    -1     1    -1     0     0     0     0     0     0 ], 'beq', 0 );  % arch junction
EqNormal.rva        = struct( 'Aeq', [  -1     0     1     0     0     0     0    -1     1     0     1 ], 'beq', 0 );  % right ventriculo-arterial connection
EqNormal.lva        = struct( 'Aeq', [   0    -1     0     0     0     1     0     1     0    -1     0 ], 'beq', 0 );  % left ventriculo-arterial connection
EqNormal.placenta   = struct( 'Aeq', [   0     0     0     0     1     0     0     0    -1     0     0 ], 'beq', 0 );  % lower body and placenta circulation
EqNormal.coronary   = struct( 'Aeq', [   0     0     0     0     0     0     0     0     0     1    -1 ], 'beq', 0 );  % coronary circulation
EqNormal.output     = struct( 'Aeq', [   1     1     0     0     0     0     0     0     0  -Rca     0 ], 'beq', 0 );  % ventricular output
EqNormal.input      = struct( 'Aeq', [   0     0     1     0     0     1     0     0     1     0  -Rcs ], 'beq', 0 );  % atrial input
EqNormal.cvo        = struct( 'Aeq', [   1     1     0     0     0     0     0     0     0     1     0 ], 'beq', 100 );% combined ventricular output
EqNormal.cai        = struct( 'Aeq', [   0     0     1     0     0     1     0     0     1     0     1 ], 'beq', 100 );% combined atrial input

switch circType
    case 'normal'
        Aeq    = [ 
                     EqNormal.heart.Aeq
                     EqNormal.pa.Aeq
                     EqNormal.arch.Aeq
                     EqNormal.rva.Aeq
                     EqNormal.lva.Aeq
                     EqNormal.placenta.Aeq
                     EqNormal.coronary.Aeq
                     EqNormal.output.Aeq
                     EqNormal.input.Aeq
                 ];
        beq    = [ 
                     EqNormal.heart.beq
                     EqNormal.pa.beq
                     EqNormal.arch.beq
                     EqNormal.rva.beq
                     EqNormal.lva.beq
                     EqNormal.placenta.beq
                     EqNormal.coronary.beq
                     EqNormal.output.beq
                     EqNormal.input.beq
                 ];
        if NVArg.isPctCvo
           Aeq = [  Aeq; EqNormal.cvo.Aeq; EqNormal.cai.Aeq ];
           beq = [  beq; EqNormal.cvo.beq; EqNormal.cai.beq ];
        end
    case 'transposition'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        Aeq    = [   EqNormal.heart.Aeq
                     EqNormal.pa.Aeq
                     EqNormal.arch.Aeq
                     0    -1     1     0     0     0     0    -1     1    -1     1  % right ventriculo-arterial connection
                    -1     0     0     0     0     1     0     1     0     0     0  % left ventriculo-arterial connection
                     EqNormal.placenta.Aeq
                     EqNormal.coronary.Aeq
                     EqNormal.output.Aeq
                     EqNormal.input.Aeq
                 ];
        beq    = [ 
                     EqNormal.heart.beq
                     EqNormal.pa.beq
                     EqNormal.arch.beq
                     0
                     0
                     EqNormal.placenta.beq
                     EqNormal.coronary.beq
                     EqNormal.output.beq
                     EqNormal.input.beq
                 ];
        if NVArg.isPctCvo
           Aeq = [  Aeq; EqNormal.cvo.Aeq; EqNormal.cai.Aeq ];
           beq = [  beq; EqNormal.cvo.beq; EqNormal.cai.beq ];
        end
    case 'ebsteins_circular_shunt'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        Aeq    = [ 
                     EqNormal.heart.Aeq
                     EqNormal.pa.Aeq
                     EqNormal.arch.Aeq
                     EqNormal.rva.Aeq
                     EqNormal.lva.Aeq
                     EqNormal.placenta.Aeq
                     EqNormal.coronary.Aeq
                     0     1     0     0     0     0     0     0     0  -Rca     0  % combined ventricular output
                    -1     0     1     0     0     1     0     0     1     0  -Rcs  % combined atrial input
                 ];
        beq    = [ 
                     EqNormal.heart.beq
                     EqNormal.pa.beq
                     EqNormal.arch.beq
                     EqNormal.rva.beq
                     EqNormal.lva.beq
                     EqNormal.placenta.beq
                     EqNormal.coronary.beq
                     0
                     0
                 ];
    case 'double_outlet_right_ventricle'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        Aeq    = [ 
                     EqNormal.heart.Aeq
                     EqNormal.pa.Aeq
                     EqNormal.arch.Aeq
                    -1    -1     1     0     0     0     0    -1     1    -1     1  % right ventriculo-arterial connection
                     0     0     0     0     0     1     0     1     0     0     0  % left ventriculo-arterial connection
                     EqNormal.placenta.Aeq
                     EqNormal.coronary.Aeq
                     EqNormal.output.Aeq
                     EqNormal.input.Aeq
                 ];
        beq    = [ 
                     EqNormal.heart.beq
                     EqNormal.pa.beq
                     EqNormal.arch.beq
                     0
                     0
                     EqNormal.placenta.beq
                     EqNormal.coronary.beq
                     EqNormal.output.beq
                     EqNormal.input.beq
                 ];
        if NVArg.isPctCvo
           Aeq = [  Aeq; EqNormal.cvo.Aeq; EqNormal.cai.Aeq ];
           beq = [  beq; EqNormal.cvo.beq; EqNormal.cai.beq ];
        end
    case 'hlhs_aortic_atresia'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        Aeq    = [ 
                    -1     0     1     0     0     1     0     0     1     0     1  % heart chambers
                     EqNormal.pa.Aeq
                     EqNormal.arch.Aeq
                     EqNormal.rva.Aeq
                     0     0     0     0     0     1     0     1     0     0     0 	% left ventriculo-arterial connection
                     EqNormal.placenta.Aeq
                     EqNormal.coronary.Aeq
                     EqNormal.input.Aeq
                 ];
        beq    = [ 
                     0
                     EqNormal.pa.beq
                     EqNormal.arch.beq
                     EqNormal.rva.beq
                     0
                     EqNormal.placenta.beq
                     EqNormal.coronary.beq
                     EqNormal.input.beq
                 ];
        if NVArg.isPctCvo
           Aeq = [  
                    Aeq 
                    1     0     0     0     0     0     0     0     0     0     0  % combined ventricular output (only one flow contributing to CVO, so remove constraint) 
                    EqNormal.cai.Aeq 
                 ];
           beq = [  
                    beq 
                    100
                    EqNormal.cai.beq 
                 ];
        end
    case 'tof_pa'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        Aeq    = [ 
                     0    -1     1     0     0     1     0     0     1    -1     1  % heart chambers
                     0     0     0    -1     0    -1     0     0     0     0     0  % pulmonary artery branch junction
                     EqNormal.arch.Aeq
                     0     0     1     0     0     0     0    -1     1     0     1  % right ventriculo-arterial connection
                     EqNormal.lva.Aeq
                     EqNormal.placenta.Aeq
                     EqNormal.coronary.Aeq
                     0     1     0     0     0     0     0     0     0  -Rca     0  % ventricular output
                     EqNormal.input.Aeq
                 ];
        beq    = [ 
                     0
                     0
                     EqNormal.arch.beq
                     0
                     EqNormal.lva.beq
                     EqNormal.placenta.beq
                     EqNormal.coronary.beq
                     0
                     EqNormal.input.beq
                 ];
        if NVArg.isPctCvo
           Aeq = [  
                    Aeq 
                    0     1     0     0     0     0     0     0     0     1     0  % combined ventricular output
                    EqNormal.cai.Aeq 
                 ];
           beq = [  
                    beq 
                    100
                    EqNormal.cai.beq 
                 ];
        end
    otherwise 
        error( 'circulation type %s not recognized, equality constraints not defined', circType )
end


%% Setup Optimization Bounds

% lower and upper bounds
if NVArg.isPctCvo
    maxBound = 100;
else
    maxBound = 1e9 * Qref.cvo;
end
switch circType 
    case {'normal','transposition','double_outlet_right_ventricle'}
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        lb =     [   0     0     0     1     0     0     0     1     0     0     0 ] * -maxBound;
        ub =     [   1     1     1     1     1     1     1     1     1     1     1 ] * +maxBound;
    case 'ebsteins_circular_shunt'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        lb =     [   1     0     0     1     0     0     0     1     0     0     0 ] * -maxBound;
        ub =     [   0     1     1     0     1     1     1     1     1     1     1 ] * +maxBound;
    case 'hlhs_aortic_atresia'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        lb =     [   0     1     0     1     0     0     0     1     0     0     0 ] * -maxBound;
        ub =     [   1     0     1     1     1     1     1     1     1     1     1 ] * +maxBound;
    case 'tof_pa'
        %        [ mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs ] 
        lb =     [   0     0     0     1     0     0     0     1     0     0     0 ] * -maxBound;
        ub =     [   0     1     1     0     1     1     1     1     1     1     1 ] * +maxBound;
    otherwise 
        error( 'circulation type %s not recognized, bounds not defined', circType )
end

           
%% Run Optimization

[x,~,exitFlag,output] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

if ( exitFlag < 1 )
    warning( output.message )
end


%% Assign Output

optFlow     = x;
unoptFlow   = x0;
wgt         = wMeas;


%% Verbose Output

if ( isVerbose )
    
    fprintf( '\n%s\n%s\n\n', mfilename, repmat('=',size(mfilename)) )
    
    fprintf( 'circulation type: %s\n\n', circType )

    % Check Constraints
    fprintf( '%-98s %s\n', 'Equality Constraints:', '(should be Aeq*x-beq == 0)' )
    for iC = 1:numel(beq) 
        AeqStr = '';
        for iV = 1:numel(vesselDesc)
            if Aeq(iC,iV) ~= 0
                if sign(Aeq(iC,iV)) == 1
                    signStr = '+';
                else
                    signStr = '-';
                end
                if abs(Aeq(iC,iV)) == 1
                    AeqStr = sprintf( '%s %s Q_%-3s', AeqStr, signStr, vesselDesc{iV} );
                else
                    AeqStr = sprintf( '%s %s %g Q_%-3s', AeqStr, signStr, abs(Aeq(iC,iV)), vesselDesc{iV} );
                end
            end
        end
        fprintf( '%-59s = %12g,     ', AeqStr, beq(iC) )
        fprintf( 'Aeq*x0-beq = %12g,     ', Aeq(iC,:)*x0.'-beq(iC) )
        fprintf( 'Aeq*x-beq = %12g     ', Aeq(iC,:)*x.'-beq(iC) )
        fprintf( '\n' )
    end
    fprintf('\n')

    fprintf( '%-98s %s\n', 'Inequality Constraints:', '(should be A*x-b    <= 0)' )
    for iC = 1:numel(b) 
        AStr = '';
        for iV = 1:numel(vesselDesc)
            if A(iC,iV) ~= 0
                if sign(A(iC,iV)) == 1
                    signStr = '+';
                else
                    signStr = '-';
                end
                if abs(A(iC,iV)) == 1
                    AStr = sprintf( '%s %s Q_%-3s', AStr, signStr, vesselDesc{iV} );
                else
                    AStr = sprintf( '%s %s %g Q_%-3s', AStr, signStr, abs(A(iC,iV)), vesselDesc{iV} );
                end
            end
        end
        fprintf( '%-58s <= %12g,     ', AStr, b(iC) )
        fprintf( 'A*x0-b     = %12g,     ', A(iC,:)*x.'-b(iC) )
        fprintf( 'A*x-b     = %12g', A(iC,:)*x0.'-b(iC) ) 
        fprintf( '\n' )
    end
    fprintf('\n')
    
    fprintf( '%-98s \n', 'Bounds:' )
    fprintf( '    %-7s ', '' )
    for iV = 1:numel(vesselDesc)
        fprintf( '%5s ', vesselDesc{iV} )
    end
    fprintf( '\b\n' )
    fprintf( '    %-7s ', 'lower' )
    fprintf( '%5i ', lb )
    fprintf( '\b\n' )
    fprintf( '    %-7s ', 'upper' )
    fprintf( '%5i ', ub )
    fprintf( '\b\n\n' )

    fprintf( '%-98s \n', 'Weights:' )
    fprintf( '    %-7s ', '' )
    for iV = 1:numel(vesselDesc)
        fprintf( '%5s ', vesselDesc{iV} )
    end
    fprintf( '\b\n' )
    fprintf( '    %-7s ', '' )
    fprintf( '%5.2f ', wMeas )
    fprintf( '\b\n\n' )
    
    % Print Final Results

    fprintf( 'Results:\n' )
    fprintf( '    %-7s %-6s %-6s %-6s %-6s %-6s\n', 'vessel', 'wgt', 'meas', 'x0', 'x', 'diff' ) 
    fprintf( '    %-7s %-6s %6s %-6s %-6s %-6s\n', '-------', '------', '------', '------', '------', '------' ) 
    for iV = 1:numel(vesselDesc) 
        fprintf( '    %-7s %6.2f %6.1f %6.1f %6.1f %6.1f\n', vesselDesc{iV}, wMeas(iV), measFlow(iV), x0(iV), x(iV), x(iV) - x0(iV) )
    end
    fprintf( '\n' )

end  % if ( isVerbose )

    
end  % balance_fetal_flow_distribution(...)
