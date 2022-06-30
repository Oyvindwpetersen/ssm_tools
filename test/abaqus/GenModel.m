%%

clc
clear all
close all

%%

% clear all
clc

% addPathTemp('C:\Work\Abaqus\Abaqusgen');

folderNameModel=[cd '\']
fileNameMainInp='ssbeam_twospan.inp';

dirMainInp=[folderNameModel fileNameMainInp];

fid = fopen(dirMainInp,'wt');

%% Gen part

genMaterial(fid,'STEEL',210e9,0.3,7850);

genPart(fid,'BEAM');

%% COL

N_node=81
xcoord1=linspace(0,40,N_node);
ycoord1=zeros(N_node,1);
zcoord1=zeros(N_node,1);

[node_no1 el_no1 matrix_node1 matrix_el1]=lineBeamMesh(xcoord1,ycoord1,zcoord1,2,1000,1,1000,1);
genNode(fid,matrix_node1,'COL1');
genElement(fid,matrix_el1,'B31','COL1');

genNset2(fid,'SUPPORT',[node_no1(1) node_no1(end)],'');
genNset2(fid,'ALL',[node_no1(1:end) ],'');


genLine(fid,'*ELEMENT, TYPE=SPRING1, ELSET=SPRING_TO_GROUND')
genLine(fid,'3000, 1030')

genLine(fid,'*SPRING, ELSET=SPRING_TO_GROUND')
genLine(fid,'2')
genLine(fid,'1e7')


genSectionBeam(fid,'COL1','STEEL','BOX',[0.1 0.1 0.01 0.01 0.01 0.01],[0 1 0]);

    
    

%% End part, assembly

genPartEnd(fid);
genAssembly(fid,'ASSEMBLY_BEAM');
genInstance(fid,'BEAM','BEAM');
genInstanceEnd(fid);



genAssemblyEnd(fid);

fclose(fid);

%% Step
fid = fopen(dirMainInp,'a+');

genBoundary2(fid,'BEAM.SUPPORT',[1 4 0],'MOD');
genBoundary2(fid,'BEAM.ALL',[3 3 0],'MOD');
genBoundary2(fid,'BEAM.ALL',[4 4 0],'MOD');

genComment(fid,'Modal',true);
genStep(fid,'','')
genFrequency(fid,20,'mass')
genFieldOutput(fid,'node',{'U' 'COORD'},' ','')
genFieldOutput(fid,'element','SF',' ','')

genStepEnd(fid)


fclose(fid);

% return
%%  Run

% cd('abaqus');
clc
tic
JobName=fileNameMainInp(1:end-4)
InpName=fileNameMainInp(1:end-4)
system(['abaqus job=' JobName ' input=' InpName ' interactive'])
toc
pause(1)

% cd ..

printAbaqusError([ folderNameModel JobName ]);


%%


dir_odb=[cd]
dir_export=dir_odb; % folder to export mat file
dir_python=[''] % folder of the pythonscript AbaqusExportModal

JobName='ssbeam_twospan.odb' % odb file name

FrequencyStepNumber=-1; % -1 gives last step in odb
ExportFileName='ssbeam_twospan' % file name of export

% path(path,dir_python);
tic
AbaqusExportModal(dir_odb,dir_export,dir_python,JobName,FrequencyStepNumber,ExportFileName)
toc


