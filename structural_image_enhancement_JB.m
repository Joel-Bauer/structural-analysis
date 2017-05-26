
function [dataout,  levels, channels, frames, dx, dy, info] = structural_image_enhancement_NormCorre_JB(yearmonthday, mouse, exp, datapath, redo_reg, bscope2)
% Uses toolboxes:
% -Distributed computing toolbox
%
%% Find files
% If date is not specified the the program will look for the experiment in
% all of the imaging session folders.
% For every dataset this should produce 3 tiff files and one mat file in
% the ProcessedStacks folder and one Backup analysis script in the StorVars
% folder
% The three files which are produced are exp00000_dataout.mat (processed stack
% in matlab variable format), exp00000_data.tif (processed stack in tiff
% format) and exp00000_data_reg.tif  (the same as exp00000_data.tif but
% with between level registration) 
clearvars -global data
global data  

%%%%
do_NormCorre_reg = 1;
%%%%

if isempty(yearmonthday)
    yearmonthday = {'*'};
    warning('Experiment date not specified.')
end

if isempty(exp)
    startexp = 1;
    endexp = 99999;
    warning('Experiment ID not specified.')
else
    %startexp = exp-1;
    %endexp = exp+1;
end

PC = getenv('computername');
if strcmp('P1-437',PC)==0
    warning('PC name is %s. The function structural_imaging_enhancement_JB may not run on this computer.',PC)
elseif strcmp('P1-437',PC)==1 
    
end

cd(datapath)

% produce a list (idxs) of all experiments of the chosen mouse ID
s = rdir([cd '\' mouse '\ImagingData\' yearmonthday{1} '**\*exp*_' '*.tif*']); %Lists all tif files
expnumbs = regexp({s.name},'(?<=exp)[0-9]*', 'match', 'once')'; %note the 'once option to get rid of cell arrays of cells and make cellfun work
[idxs b] = unique(expnumbs);
for i=1:size(idxs,1)
    if str2double(idxs(i))== exp;
        experiment_directory = fileparts(s(b(i)).name);
    end
end

% If date of exp was not specified it will be specified here according to
% the location of the exp directory. 
if cell2mat(yearmonthday) == ('*')
    yearmonthday{1} = experiment_directory(end-9:end);
end

cd(experiment_directory) %move to directory containing data

%% Load experiment file and image enhancement parameters
if redo_reg == 1 
    try
        delete([cd '\exp' num2str(exp) '_regParameters.mat']);
    end
    [data, info, levels, channels, frames] = loadSI_local(num2str(exp), cd, 0); %variable frames indicates frames per level per channel
else
    try
        load([cd '\exp' num2str(exp) '_regParameters.mat']); % This folder contains the parameters from previouse image enhancement runs
        [data, ~, levels, channels, frames] = loadSI_local(num2str(exp), cd, 0);
    catch
        disp('No previouse image enhancement parameters found')
        [data, info, levels, channels, frames] = loadSI_local(num2str(exp), cd, 0);
        redo_reg = 1;
    end
end

if levels == 4 && isempty(info)
    return
end


clear b startexp endexp expnumbs i idxs % these varialbes are no longer required


%% Image enhancement 
% What parts of the work-flow should be done or redone
% Alows easier troubleshooting to skip specific enhancement methods

if redo_reg == 1
    do_correct_line_shift = 1; redo_correct_line_shift = 1; 
elseif redo_reg == 0
    do_correct_line_shift = 1; redo_correct_line_shift = 0;
end

% If that stack has a zoom higher than 5 and was not frame averaged then
% within level registration and between level registration is applied.
if info.zoom > 5 && info.frame_averaging == 1
    if redo_reg == 1
        do_xcorr_frame_reg    = 1; redo_xcorr_frame_reg    = 1;
        do_LucasKanade_reg    = 1; redo_LucasKanade_reg    = 1;
        do_ImgJ_stack_reg     = 0;
        if do_NormCorre_reg == 1
        	do_xcorr_frame_reg    = 0; redo_xcorr_frame_reg    = 0;
            do_LucasKanade_reg    = 0; redo_LucasKanade_reg    = 0;
            do_NormCorre_reg      = 1; redo_NormCorre_reg      = 1;
        else
            do_NormCorre_reg      = 0; redo_NormCorre_reg      = 0;
        end
    elseif redo_reg == 0
        do_xcorr_frame_reg    = 1; redo_xcorr_frame_reg    = 0;
        do_LucasKanade_reg    = 1; redo_LucasKanade_reg    = 0;
        do_NormCorre_reg      = 0; redo_NormCorre_reg      = 0;
        do_ImgJ_stack_reg     = 0;
    end
    disp('Single frames will be registered both within and between levels')
elseif info.zoom < 5 && info.frame_averaging == 1
    do_xcorr_frame_reg    = 0; redo_xcorr_frame_reg    = 0;
    do_LucasKanade_reg    = 0; redo_LucasKanade_reg    = 0;
    do_NormCorre_reg      = 0; redo_NormCorre_reg      = 0;
    do_ImgJ_stack_reg     = 0;
    disp('No frame registration will be done. Averaging all frames per level!')
    for i = 1:levels
        frame_index = (i-1)*frames+1:i*frames;
        data_temp(:,:,i) = uint16(mean(data(:,:,frame_index),3));
    end
    %clear data; %dont do this as data will no longer refer to the global variable for this work space.
    data = data_temp; clear data_temp
else
    do_xcorr_frame_reg    = 0; redo_xcorr_frame_reg    = 0;
    do_LucasKanade_reg    = 0; redo_LucasKanade_reg    = 0;
    do_NormCorre_reg      = 0; redo_NormCorre_reg      = 0;
    do_ImgJ_stack_reg     = 0;
    disp('No frame registration will be done')
end

% Manualy enable CANDLE denoising and pixelvalue thresholding
do_CANDLE_denoising = 0;
threshold_noise = 1;
do_median_frame_avg =0; %Does not work at all!!

%% Threshold pixelvalues to get rid of cohirent noise
% Threshold pixel values at the mode pixel value to remove piezo induced
% noise
% DOES NOT WORK WITH LOW AMPLITUDE DATA!!!
% 
% 
% if size(data,3) == info.level*info.frames_per_level_aquired; %info.frame_averaging == 0 % change back if it failes here using bscope1
%     total_frames = frames*levels;
% elseif size(data,3) == info.level; %info.frame_averaging == 1
%     total_frames = levels;
% else
%     total_frames = levels;
% end

if threshold_noise == 1
    tic
    for i = 1:size(data,3); %(total_frames)
        one_frame = data(:,:,i);
        frame_mean(i) = mean(one_frame(:));
        frame_std(i) = std(double(one_frame(:)));
        frame_threshold(i) = frame_mean(i)+frame_std(i);
    end
    threshold = mean(frame_threshold);
    for i = 1:size(data,3); %(total_frames)
        for k = 1:size(data,1)
            line_logical=data(k,:,i)>threshold;
            data(k,:,i) = double(data(k,:,i)).*line_logical; %
        end 
    end
end; 
clear frame_mean frame_std threshold one_frame

%% Line shift correction 
if do_correct_line_shift == 1
    if ~exist('line_shift','var') || isempty(line_shift)
        tic;
        if redo_correct_line_shift == 0
            disp('Unable to reuse prevouse line shift estimation.')
        end
        disp('Estimating and correcting line shift.')
        [line_shift]=correct_line_shift_local([],[],bscope2); %variable data is passed into function as a global variable to save RAM
        disp(['Estimating and correcting line shift took ' num2str(toc) ' sec.'])  
    elseif exist('line_shift','var') && redo_correct_line_shift == 1
        tic;
        disp('Recalculating and correcting line shift. (previouse line shift estimatino not reused)')
        [line_shift]=correct_line_shift_local([],[],bscope2); 
        disp(['Estimating and correcting line shift took ' num2str(toc) ' sec.'])
    elseif exist('line_shift','var') && redo_correct_line_shift == 0
        tic;
        disp(['Using line shift estimation from previouse run (' num2str(line_shift) ').'])
        correct_line_shift_local([], line_shift,bscope2); 
        disp(['Correcting line shift took ' num2str(toc) ' sec.'])
    end
else
    line_shift = [];
    disp('No line shift correction done.')
end

%% Xcorr frame registration
if do_xcorr_frame_reg == 1 && frames>1
    if redo_xcorr_frame_reg == 0 && (exist('dx','var') || exist('dy','var') || ~isempty(dx) || ~isempty(dy))
        disp('Using within level Xcorr registration parameters from previouse run.')
    elseif redo_xcorr_frame_reg == 1 
        disp('Redoing Xcorr within level registration.')
        dx=[]; dy=[];
        for i = 1:levels
            disp(['Layer ' num2str(i) '/' num2str(levels)])
%             % a progress update using this template 
%             exp =90:102;
%             fprintf(1,'For loop progression (of %s):       ',num2str(exp(length(exp))));
%             for i = 1:length(exp)
%                 fprintf(1,'\b\b\b\b\b\b%5.5s\n',num2str(exp(i)));
%                 for k = 1:10
%                     fprintf('For loop 2: %5.5s', num2str(k));pause(0.05);
%                     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
%                 end
%             end
            frame_index = (i-1)*frames+1:i*frames;
            [dx_level, dy_level] = cpu_register_local(data(:,:,frame_index)); % Check this function does the same as gpu_register
            dx = cat(1, dx, dx_level);
            dy = cat(1, dy, dy_level);
        end 
        % dx=round(dx-mean(dx)/2); % not sure what this does but its in the register_multilayer_gpu function ask tobias.
        % dy=round(dy-mean(dy)/2); % not sure what this does but its in the register_multilayer_gpu function
    elseif redo_xcorr_frame_reg == 0 && (~exist('dx','var') || ~exist('dy','var') || isempty(dx) || isempty(dy))
        disp('Unable to find prevouse Xcorr registration parameters. Redoing within level registration.')
        dx=[]; dy=[];
        for i = 1:levels
            disp(['Layer ' num2str(i)])
            frame_index = (i-1)*frames+1:i*frames;
            [dx_level, dy_level] = cpu_register_local(data(:,:,frame_index)); % Check this function does the same as gpu_register
            dx = cat(1, dx, dx_level);
            dy = cat(1, dy, dy_level);
        end 
    end
    if (range(dx)+range(dy)) == 0;
        disp('Xcorr registration detected no movement. Xcorr based frame shifting was therefore skipped');
    else
        disp('Shifting frames based on Xcorr registration parameters');
        data = shift_data_TR(data,dx,dy);
    end
else
    disp('No within level frame registration done.')
    dx=[]; dy=[];
end

%% Lucas Kanade frame registration (addapted from the function
% register_multilayer_gpu.m)
if do_LucasKanade_reg == 1&& frames>1
    if redo_LucasKanade_reg == 0 && (exist('LKdx','var') || exist('LKdy','var') || ~isempty(dxLK) || ~isempty(dyLK))
        disp('Using within level Lucas Kanad registration parameters from previouse run.')
    elseif redo_LucasKanade_reg == 1 
        disp('Redoing Lucas Kanad registration.')
        for i = 1:levels
            disp(['Layer ' num2str(i) '/' num2str(levels)]) 
            frame_index = (i-1)*frames+1:i*frames;
            template = mean(data(:,:,frame_index),3); % Makes new template for every level
            [LKdx_level, LKdy_level] = LucasKanade_gpu(data(:,:,frame_index), template); % Check this function does the same as gpu_register % anoying output variable format!!!!
            LKdx{i} = LKdx_level; % Yes im making a cell of cells of arrays. DEAL WITH IT!!!
            LKdy{i} = LKdy_level;
        end 
        % dx=round(dx-mean(dx)/2); % not sure what this does but its in the register_multilayer_gpu function so iv left it in here as a note.
        % dy=round(dy-mean(dy)/2); % same with this
    elseif redo_LucasKanade_reg == 0 && (~exist('LKdx','var') || ~exist('LKdy','var') || isempty(LKdx) || isempty(LKdy))
        disp('Unable to find prevouse Lucas Kanad registration parameters. Redoing within level registration.')
        for i = 1:levels
            disp(['Layer ' num2str(i) '/' num2str(levels)]) 
            frame_index = (i-1)*frames+1:i*frames;
            template = mean(data(:,:,frame_index),3); % Makes new template for every level
            [LKdx_level, LKdy_level] = LucasKanade_gpu(data(:,:,frame_index), template); % Check this function does the same as gpu_register % anoying output variable format!!!!
            LKdx{i} = LKdx_level; % Yes im making a cell of cells of arrays. DEAL WITH IT!!!
            LKdy{i} = LKdy_level;
        end 
    end
    for i = 1:levels
        LKdx_level = LKdx{i}; % Unpacking the cell which is a cell of cells of arrays, giving a single cell of arrays.
        LKdy_level = LKdx{i};
        disp(['Layer ' num2str(i) '/' num2str(levels)]) 
        frame_index = (i-1)*frames+1:i*frames;
        data(:,:,frame_index) = shift_LucasKanade_gpu_JB(data(:,:,frame_index),LKdx_level,LKdy_level); % _JB is a version of the function which does not use the GPU
        % Need to test if creating a new data set and then overwriting data is faster than changing each line
    end
else
    LKdx = []; LKdy = [];
    disp('No within level Lucas Kanade frame registration done.')
end

%% NormCorre non-rigid registration
if do_NormCorre_reg == 1&& frames>1
    if redo_NormCorre_reg == 1
        disp('doing NormCorre registration.')
        d1 = size(data,1);
        d2 = size(data,2);
        grid_size= [32,32];
        mot_uf = 4;
        bin_width = 50;
        max_shift = 15;
        max_dev = 3;
        us_fac = 50 ;       
        
        for i = 1:levels
            disp(['Layer ' num2str(i) '/' num2str(levels)]) 
            frame_index = (i-1)*frames+1:i*frames;
            options_nonrigid = NoRMCorreSetParms('d1',d1,'d2',d2,'grid_size',grid_size,'mot_uf',mot_uf,'bin_width',bin_width,'max_shift',max_shift,'max_dev',max_dev,'us_fac',us_fac);
            tic; [data(:,:,frame_index),shifts2,~] = normcorre_batch(data(:,:,frame_index),options_nonrigid); toc
        end
    end
end
        
%% Average Within level frames

if size(data,3)~= levels && ~do_median_frame_avg;
    disp('Doing mean frame averaging instead of median frame averaging')
    for i = 1:levels
        frame_index = (i-1)*frames+1:i*frames;
        dataout(:,:,i) = uint16(mean(data(:,:,frame_index),3));
    end
elseif do_median_frame_avg && size(data,3)~= levels; % not tested yet but try and compare difference between mean and median averaging 
    disp('Doing median frame averaging instead of mean frame averaging')
    for i = 1:levels
        frame_index = (i-1)*frames+1:i*frames;
        dataout(:,:,i) = uint16(median(data(:,:,frame_index),3));
    end
else
    dataout = data;
end
clear data

%% ImageJ Z stack registration 
%(addapted from the function register_multilayer_gpu.m)

if do_ImgJ_stack_reg  == 1
    try
        MIJpath = 'C:\Users\jbauer\Google Drive\ImageJ\mij.jar';
        IJpath = 'C:\Users\jbauer\Google Drive\ImageJ\';
        javaaddpath 'C:\MATLAB\R2013b\64bit\java\mij.jar'; 
        javaaddpath 'C:\MATLAB\R2013b\64bit\java\ij.jar'; 
        MIJ.start(IJpath);
        MIJ.run('Close All');
    catch
        disp('File path for ImageJ not found. Z-Stack registration skipped')
    end
    if ~isempty(MIJpath) || ~isempty(IJpath) || ~isempty(MIJpath)
        try
            MIJ.createImage(uint16(dataout)); 
            MIJ.run('StackReg ', ['transformation=Translation']);
            dataout_reg = MIJ.getCurrentImage;
            dataout_reg = int16(dataout_reg); % image j aligned stack variable
            MIJ.run('Close');
            MIJ.exit
        catch
            disp('ImageJ crashed - run again manually on registered stack')
        end
    end
    
else 
    dataout_reg = dataout;
end

% Reload header

img_chuncks=dir([cd '\exp' num2str(exp) '**.tif']);
temp=imfinfo(img_chuncks(1).name);
img_header = temp(1);


%% Make structural analysis data folder
mkdir('StructuralAnalysis'); 
cd 'StructuralAnalysis'
mkdir('ProcessedStacks')
mkdir('StoreVars')

%% CANDLE denoising (saves output Tiff in StructuralAnalysis folder
%replace with 3D median filter
if do_CANDLE_denoising  == 1; 
    cd 'ProcessedStacks';
    try; delete([cd '\exp' num2str(exp) '_data_denoised.tif']); end;
    FileNameDenoised = [cd '\exp' num2str(exp) '_data_denoised.tif']; % Def .tif dataout file name
    searchradius = 3;% 
    patchradius = 1;    %A radius of 1 voxel (i.e., patch of 3x3x3 voxels) is usually appropriate for a greater view of fine structures. 
                        %A radius of 2 voxels (i.e., patch of 5x5x5 voxels) should be favoured for a restricted field of view of large structures.
    beta = 0.1;%
    background = 1;%
    if do_CANDLE_denoising == 1;
        [dataout_denoised]=CANDLEfilter_local(dataout_reg, FileNameDenoised, searchradius, patchradius, beta, background);
    end
    cd ..
else
    disp('CANDLE denoising skipped')
end



%% Save output data as .mat and .tiff
% Save shift parameters
% Return to ImagingData folder
dbstop if error

cd ..
disp('Saving dataout'); 
if do_NormCorre_reg == 0
    save([cd '\exp' num2str(exp) '_regParameters.mat'], 'info', 'line_shift', 'dx', 'dy', 'LKdx', 'LKdy', '-mat') % Saving registration parameters
end
cd 'StructuralAnalysis\ProcessedStacks';

% Save data output as tiff and mat files in StructuralAnalysis folder
FileNameMat = [cd '\exp' num2str(exp) '_dataout_NormCorre.mat']; % Def .mat dataout file name
FileNameTif = [cd '\exp' num2str(exp) '_dataout_NormCorre.tif']; % Def .tif dataout file name without ImageJ reg
%FileNameRegTif = [cd '\exp' num2str(exp) '_dataout_reg.tif']; % Def .tif dataout file name with ImageJ reg

try; delete(FileNameMat); end; %delete pre-existing dataout files!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Its not deleting old files
try; delete([cd '\exp' num2str(exp) '_data_NormCorre.mat']); end;
try; delete(FileNameTif); end;
%try; delete(FileNameRegTif); end;
try; delete([cd '\exp' num2str(exp) '_data_NormCorre.tif']); end;
%try; delete([cd '\exp' num2str(exp) '_data_reg.tif']); end;
pause(5)

% Save averaged data (no ImageJ registration). Tiff
% add saving of header
try
    dataoutINT=uint16(dataout); % convert output data to integer (to make sure it is not a double when converted to TIFF 
end


%readload header!! to save it to new fun
try
    % imwrite(dataoutINT(:,:,1),FileNameTif,'tiff','writemode','overwrite','compression','none'); %save first frame as *dataout.tif
    imwrite2tif(dataoutINT(:,:,1), img_header(1), FileNameTif, 'uint16') % This function creates a new Tif file INCLUDING the header
    for K=2:length(dataoutINT(1, 1, :))
       imwrite(dataoutINT(:, :, K), FileNameTif, 'tiff', 'writemode', 'append',  'compression','none'); % save subsequent frames in *dataout.tif
    end
catch
    warning('Tiff files may not have been saved due to stupid "You may not have write permission" ')
end

% Save averaged data with ImageJ registration. Tiff
try
    dataout_regINT=uint16(dataout_reg); % convert output data to integer (to make sure it is not a double when converted to TIFF 
end
% try
%     % imwrite(dataout_regINT(:,:,1),FileNameRegTif,'tiff','writemode','overwrite','compression','none'); %save first frame as *dataout.tif
%     imwrite2tif(dataout_regINT(:,:,1), img_header(1), FileNameRegTif, 'uint16') % This function creates a new Tif file INCLUDING the header
%     for K=2:length(dataout_regINT(1, 1, :))
%        imwrite(dataout_regINT(:, :, K), FileNameRegTif, 'tiff', 'writemode', 'append',  'compression','none'); % save subsequent frames in *dataout.tif
%     end
% catch
%     warning('Check writing permission in current directory. Try using imwrite instead of imwrite2tif. ')
% end


% Make a file in the Structural Analysis folder indicating how the analysed
% stacks where aquired.
cd ..
FileNameAndLocation=[mfilename('fullpath')];
newbackup=[cd '\StoreVars\' mfilename '_backup_exp' num2str(exp) '.m'];
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);


% For every dataset this should produce 3 tiff files and one mat file in
% the ProcessedStacks folder and one Backup analysis script in the StorVars
% folder

dbclear if error
end

function [data, info, levels, channels, frames] = loadSI_local(exp, diri, fct)


%% get infromation from Tif file header
files = dir([diri '\*' exp '*.tif*']);

tio = Tiff([diri '\' files(1).name], 'r');
info.ImageDescription = tio.getTag('ImageDescription');

info.Width = tio.getTag('ImageWidth');
info.Height = tio.getTag('ImageLength');
info.FileID = exp;

info.tnumber_of_frames =  regexp(info(1).ImageDescription,'(?<=loggingFramesPerFile = )\d+\.?\d*', 'match');
info.tnumber_of_frames = str2num(info.tnumber_of_frames{1});

info.tchannels = regexp(info(1).ImageDescription,'(?<=channelsSave = )\S*', 'match');
info.tchannels = str2num(info.tchannels{1});
info.tchannels = length(info.tchannels);


try %try to get the zoom factor of the Tif
    info.objective = str2num(cell2mat(regexp(info.ImageDescription,'(?<=userFunctionsCfg__4.Arguments = {)[0-9]*', 'match')));
    [info.fieldsize, info.zoom,~, info.zoom1_field, info.pixelsize] = getzoom([diri '\' files(1).name], info.objective, info);
catch
    warning('Unable to aquire zoomfactor with function getzoom.') 
end

info.frames_per_level_aquired = str2num(cell2mat(regexp(info(1).ImageDescription,'(?<=acqNumFrames = )[0-9]*', 'match')));
info.frame_averaging = str2num(cell2mat(regexp(info(1).ImageDescription,'(?<=acqNumAveragedFrames = )[0-9]*', 'match')));
frames = info.frames_per_level_aquired/info.frame_averaging;
info.frames_per_level_saved = frames;



%get number of saved channels in the Tif file
channelsSave_str= info(1).ImageDescription(...
  findstr('channelsSave',info(1).ImageDescription)+15:...
  findstr('channelsSave',info(1).ImageDescription)+20);
if strcmp(channelsSave_str(1),'2')
    channels=1;
elseif strcmp(channelsSave_str(1:5),'[1;2]')
    channels=2;
else
    warning('Number of channels in Tif can not be detected, and is therefore assumed to be 1')
    channels=1;
end
info.channels = channels;


% get number of levels in the Tif file
levels = str2num(cell2mat(regexp(info(1).ImageDescription,'(?<=stackNumSlices = )[0-9]*', 'match')));
info.level = levels;

if levels == 4
    disp('Number of levels is 4.')
    disp(['Assuming recording was functional and canceling Structural image enhancement for exp' num2str(exp)])
    fct = 1;
    data = [];
    info= [];
    levels= 4;
    channels= [];
    frames= [];
    return
end
    
if fct == 0%%%%%%%%%%%%%%%%%
    PiezoMotor_stack = str2num(cell2mat(regexp(info(1).ImageDescription,'(?<=motorSecondMotorZEnable = ).', 'match')));
    info.PiezoMotor_stack = PiezoMotor_stack;
    if PiezoMotor_stack == 1
        string_temp = regexp(info(1).ImageDescription,'(?<=stackZStepSize = ).........', 'match');
        string_temp = strsplit(string_temp {1}, '\n');
        Zstack_step_size=str2num(string_temp {1});
        info.Zstack_step_size = Zstack_step_size;
        info.Zstack_volume = info.level*Zstack_step_size; 
    elseif PiezoMotor_stack == 0
        string_temp = regexp(info(1).ImageDescription,'(?<=stackZStepSize = -).........', 'match');
        string_temp = strsplit(string_temp {1}, '\n');
        Zstack_step_size=str2num(string_temp {1});
        info.Zstack_step_size = Zstack_step_size;
        info.Zstack_volume = info.level*Zstack_step_size;
    end
end

for i = 1:length(files);
    disp(['loading Header ' diri '\' files(i).name]);
    al{i} = imfinfo([diri '\' files(i).name]);
    chunksize(i) = al{i}(1).FileSize;
    b(i) = length(al{i});
    if chunksize(i) > 3.5*1024*1024*1024;
        disp('- - - - - - - -')
        disp('Tiff most likely corrupted (>3.5GB)')
        disp('- - - - - - - -')
        return
    end
end
img_header = al{1};

info.totalframes =sum(b);


%% Memory check
[~,sysv] = memory;
available_RAM = sysv.PhysicalMemory.Available;
gigasize = sum(chunksize);

if gigasize > available_RAM*0.95 %Error message if file size is over 95% of available RAM
    disp('- - - - - - - - - - - - - - - - - - -')
    disp('FILE MAY BE TOO BIG FOR RAM LOADING! ')
    disp(['Total file size is:  ' num2str(gigasize/1e+9), ' GB'])
    disp(['Available RAM:       ', num2str(available_RAM/1e+9) , ' GB'])
    disp('- - - - - - - - - - - - - - - - - - -')
end



%% load Tif files
data = zeros([info(1).Height info(1).Width sum(b)], 'int16');% preallocate data m?????????????????????????????????

filecat = 0;
% for i = 1:length(files);
%     disp(['loading Tiff ' diri '\' files(i).name]);
%     data(:,:,filecat+1:filecat+b(i)) = load_tiff_TR_TIFF([diri '\' files(i).name],al{i}, 0, 1);
%     filecat = filecat + b(i);
% end
for i = 1:length(files);
    disp(['loading Tiff ' diri '\' files(i).name]);
    data_range=(filecat+1:filecat+b(i));
    if strcmp(diri(1),'D')
        data(:,:,data_range) = load_tiff_TR_JB([cd '\' files(i).name],1,al{i},0,1); % also loads the negative values wwithout cutting them off
    elseif strcmp(diri(1),'I')
        tic; 
        data(:,:,data_range) = ReadTiffStack([cd '\' files(i).name],al{i});
        disp(['Elapsed time is' num2str(toc) ' seconds.'])
    end
    filecat = filecat + b(i);
end
clear al; % al can be quite big and is not used again so it is deleted here to save RAM

% if min(data(:))<0
%     data = data - min(data(:));  % rescales data to get rid of netagive values.
% end

end

function [dx, dy] = cpu_register_local(data)
clear template
dx = [];
dy = [];
for i = 1:round(size(data,3)/10);
    if size(data,3)<10
       template=(mean(data,3)); 
    else
        try
            template = mean(data(:,:,(i-1)*10+1:i*10),3); 
        catch
            template = mean(data(:,:,(i-1)*10+1:end),3);
        end 
    end
    
    bound = 0.05;
    template=template-mean(mean(template)); %?

    % determine how much of the images to use for alignment, the larger the
    % boundary the less pixels are used for alignment and the faster the
    % algorithm runs.

    boundary=round(bound*max(size(data(:,:,1))));
    template=template(boundary+1:end-boundary,boundary+1:end-boundary);
    F_template_size = [size(template,1) size(template,2)];
    F_in_size = [size(data,1)-2*boundary size(data,2)-2*boundary];
    outsize = F_template_size + F_in_size - 1;
    low_pass_thresh=round(min(F_template_size)/4);
    high_pass_thresh=round(min(F_template_size)/40);
    template = im2double(template);

    % fourier transform and band pass filter the template
    F_template = fft2(rot90(template,2),outsize(1),outsize(2));

    F_template([1:high_pass_thresh size(F_template,1)-high_pass_thresh+2:size(F_template,1)],:) = 0;
    F_template(:,[1:high_pass_thresh size(F_template,2)-high_pass_thresh+2:size(F_template,2)]) = 0;

    F_template(low_pass_thresh+2:end-low_pass_thresh,:)=0;
    F_template(:,low_pass_thresh+2:end-low_pass_thresh)=0;


    try
        [dxt, dyt] = register_frames_CPU_gfor(data(:,:,(i-1)*10+1:i*10),F_template, boundary,low_pass_thresh,high_pass_thresh,outsize); %, template); 
    catch
        [dxt, dyt] = register_frames_CPU_gfor(data(:,:,(i-1)*10+1:end),F_template, boundary,low_pass_thresh,high_pass_thresh,outsize); %, template);
    end
    dx = cat(1, dx, dxt);
    dy = cat(1, dy, dyt);
end

fprintf('\n') % print empty line?
end

function [shift]=correct_line_shift_local(template,shift,bscope2)
% estimates and corrects the shift between successive scan lines
global data 
% Define boundry to remove sides of image (center of the image can have a different line shift to the sides)
if ~exist('template','var') || isempty(template)
    if size(data,1)== 512
        template_boundries = round(size(data,1)*0.2):round(size(data,1)-size(data,1)*0.2); %boudrie size is too big for 1024*1024 images
        template=mean(data(4:end,template_boundries,:),3); %to ignore first few frames i changed this without testing from template=mean(data(:,template_boundries,:),3);
    elseif size(data,1) == 1024
        template_boundries = round(size(data,1)*0.3):round(size(data,1)-size(data,1)*0.3); %boudrie size is too big for 1024*1024 images
        template=mean(data(template_boundries,template_boundries,round(size(data,3)*0.3):round(size(data,3)*0.4)),3);
    end
end

twoDdata=0;
if size(data,3)==1
    tdata(:,:,1)=data;
    tdata(:,:,2)=data;
    twoDdata=1;
    data=tdata;
end

if nargin<3 || isempty(shift)
    shift=estimate_line_shift_local(template);
    if bscope2 == 1 && size(data,1)==512 
        shift=shift*-1; % for some wird reason the line shift screws up otherwise. dont understand why. Dont think this was necessary for bscope1 data
    end
end

function [shift]=estimate_line_shift_local(template) %estimate_line_shift function as nested function
% this function attempts to estimate the shift between alternating lines 
disp('Estimating line shift');

for ind=1:size(template,1)-1
    if rem(ind,2)==1
        line_xcorr(ind,:) = xcorr(template(ind,:)-mean(template(ind,:)),template(ind+1,:)-mean(template(ind+1,:)),'unbiased');
    else
        line_xcorr(ind,:) = flipdim(xcorr(template(ind,:)-mean(template(ind,:)),template(ind+1,:)-mean(template(ind+1,:)),'unbiased'),2);
    end
end

boundary=floor(size(template,2)/2);
[~,b]=max(mean(line_xcorr(:,boundary+1:end-boundary)));
shift=round(size(line_xcorr,2)/2)-b-boundary;
disp(['Estimated line shift is ' num2str(shift) '.  Correcting shift.']);
end

clear template % To save RAM



if shift~=0
    for i=1:size(data,1)
        if rem(i,2)==1
            data(i,:,:)=circshift(squeeze(data(i,:,:)),[shift 0]);
        end
    end
end

if twoDdata
    data=data(:,:,1);
end

end

function [imgOut] = CANDLEfilter_local(data, filenameOUT, searchradius, patchradius, beta, background)
imgIn = data;
nbfile = 1;
if nargin <2 || isempty(filenameOUT)
    filenameOUT= [cd '\expUNKOWN' '_data_denoised.tif']';%
end
if nargin < 3 % defalt variables
    searchradius = 3;% 
    patchradius = 2;%
    beta = 0.1;%
    background = 1;%
end
dim = size(data,3);%


for f = 1 : nbfile

%         figure
% 
%         if(nbfile>1)
%             filenamein = namein{f};
%         else
%             filenamein = namein;
%         end
% 
%         disp(['Input file : ', fullfile(pathin, filenamein)])
%         [pathstr, name_s, ext]=fileparts(fullfile(pathin, filenamein));
%         nout=[name_s suffix ext];
%         pathout = pathin;
% 
%         disp(['Output file: ', fullfile(pathout, nout)])
% 
%         info = imfinfo([pathin filenamein]);
%         dim = numel(info);
% 
%         com = sprintf('\nNumber of slices: %d \n', dim);
%         disp(com)
% 
%         % Image reading
%         for i = 1:dim
%             imgIn(:,:,i) = imread([pathin filenamein],i);
%         end



        s=size(imgIn);
        if (size(s)~=3)
            error('Input image must be a 3-D array.')
        end

        fact = 64-searchradius;

        substack = ceil(dim/fact);

        if (substack==2)
            fact=fact*2;
            substack=1;
        end

        if (substack==1)
            fact=dim;
        end



        for ind=1 : substack

            % To avoid out of memory the stack is processed by substacks of 64 images
            s_init = 1 + (ind-1)*fact; % beguinning of the substack
            s_end =  ind*fact; % end of the substack


            % Limits of the substack
            if (s_end>dim) s_end = dim; end
            if (s_end-s_init<12) s_init = s_end-12; end
            if (s_init<1) s_init = 1; end

            % Padding between substacks
            if ((s_init-searchradius)<=1)
                pad_init = 0;
            else
                pad_init = searchradius;
            end

            if((s_end+searchradius)>dim)
                pad_end = 0;
            else
                pad_end = searchradius;
            end



            com = sprintf('Processing of the slices: %d to %d', s_init,s_end);
            disp(com)

            % Subimg to denoise
            img = single(imgIn(:,:,s_init-pad_init:s_end+pad_end));


            com = sprintf('\n\t Preprocessing: 3D Median filter');
            disp(com)
            tic
            fimgMed = median3D(img,1);
            t = toc;
            com = sprintf('\t Elapsed time: %0.1f s', t);
            disp(com)

            % Background detection
            if (background==1)
                mini = min(fimgMed(:));%?
                maxi = max(fimgMed(:));%?
                average_val = mean(fimgMed(:));
                delta = (maxi-mini)/10000;
                k=3;
                N=k^3;
                ConvMed=(convn(fimgMed,ones(k,k,k),'same')/N);
                bins = double(mini):delta:double(maxi);
                [nb,histb] = histc(ConvMed(:),bins);
                [val loc] = max(nb);
                Threshold = bins(loc);
                smask = size(fimgMed);
                mask = ones(smask);
                map = find(ConvMed<(Threshold+5*delta));
                mask(map) = 0;
            else
                smask = size(fimgMed);
                mask = ones(smask);
            end



            % Anscombe transform to convert Poisson noise into Gaussian noise
            img = 2 * sqrt(img + (3/8));
            fimgMed = 2 * sqrt(fimgMed + (3/8));

            % Estimation of std for Gaussian noise
            com = sprintf('\n\t Noise estimation: Wavelet-based local estimation');
            disp(com)
            tic
            [MAP,h] = GaussianSTD_MAP(img,2*searchradius);
            t = toc;
            com = sprintf('\t Elapsed time: %0.1f s', t);
            disp(com)



            % Denoising
            com = sprintf('\n\t Denoising: 3D Optimized Nonlocal means filter');
            disp(com)
            tic
            fimg=adapt_onlm_collaborative(single(img),searchradius, patchradius,single(MAP),beta,single(fimgMed),single(mask));
            t = toc;
            com = sprintf('\t Elapsed time: %0.1f s \n', t);
            disp(com)



            % Optimal Inverse Anscombe Transform
            fimg = OVST(fimg);
            img =  OVST(img);
            fimgMed = OVST(fimgMed);


            % Convertion  
            % Since I do not want libTIFF dependencies (version issues, OS issues...), I cannot write in 24bit.
            % So I decided to write the output image in 16bit.   
            imgOut(:,:,s_init:s_end) = uint16(fimg(:,:,pad_init+1:end-pad_end));



            % Substack Display
            slice = floor((s_end-s_init)/2);
            mini = min(imgOut(:));
            maxi = max(imgOut(:));
            subplot(1,2,1)
            imagesc(img(:,:,slice)); colormap('bone');
            axis image
            axis off
            title('Input image');
            subplot(1,2,2)
            imagesc(fimg(:,:,slice),[mini maxi-0.1*maxi]); colormap('bone');
            axis image
            axis off
            title('Denoised with CANDLE');
            drawnow;

        end 

        % Display
        mini = min(imgOut(:));
        maxi = max(imgOut(:));
        MIPnoisy = max(imgIn, [], 3);
        MIPdenoised = max(imgOut, [], 3);
        figure;
        subplot(1,2,1)
        imagesc(MIPnoisy); colormap('bone');
        axis image
        axis off
        title('MIP of Input image');
        subplot(1,2,2)
        imagesc(MIPdenoised,[mini maxi-0.1*maxi]); colormap('bone');
        axis image
        axis off
        title('MIP of Denoised image');
        % Save Tiff file
        imwrite(imgOut(:,:,1),filenameOUT,'Compression','none');
        for k = 2:dim
            imwrite(imgOut(:,:,k),filenameOUT,'Compression','none','writemode', 'append');
        end
end
end

