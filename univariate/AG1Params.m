function par = AG1Params(subject)

% sets up parameters for batching of analyses
% written by amg
% based on code by jbh


fprintf('\nEstablishing Parameters for Subject %s... \n', subject);

spm_defaults; %set-up defaults 
par.substr = subject;

par.Tasks = {'Ret'}; %either 'Enc' or 'Ret'

%% experiment-wide directory names
par.pardir = '/Users/gordonam/Studies/AG1/';
par.scriptsdir = fullfile(par.pardir,'scripts');
par.groupbehavdir = fullfile(par.pardir, 'behav_data');
par.fmridir = fullfile(par.pardir, 'fmri_data2');

%% subject-specific directory names
par.subdir = fullfile(par.fmridir, subject);
par.anatdir = fullfile(par.subdir, 'anat');
par.funcdir = fullfile(par.subdir, 'functional');
par.behavdir = fullfile(par.subdir, 'behav');
par.rawdir = fullfile(par.subdir, 'raw');
par.logdir = fullfile(par.subdir, 'logfiles');
par.artrepdir = fullfile(par.subdir, 'artRepair');

%% Encoding-specific variables
par.Enc.dir = fullfile(par.subdir, 'EncAnalysis');
par.Enc.dur = 420;
par.Enc.nTRs = 210;
par.Enc.nTrials = 180;
par.Enc.Scans = [1 2 3 4 5 6 ];
par.Enc.nSess = 6;

%% Retrieval-specific variables
par.Ret.dir = fullfile(par.subdir, 'RetAnalysis');
par.Ret.dur = 280;
par.Ret.nTRs = 140;
par.Ret.nTrials = 168;
par.Ret.Scans = [7 8 9 10 11 12];
par.Ret.nSess = 6;

%% which scans were collected
par.scans_to_include = [1 2 3 4 5 6 7 8 9 10 11 12];

if strcmp(subject, 'ag1_071409') %missing fMRI data for scan04
    par.scans_to_include = [1 2 3 5 6 7 8 9 10 11 12];
    par.EncScans = [1 2 3 4 5];
    par.RetScans = [6 7 8 9 10 11];
end

%% artifact detection parameters
par.art.motThresh = .5;
par.art.sigThresh = 4;

%% scan acquisition info
par.numscans = length(par.scans_to_include); %number of scans
par.TR = 2; %TR in seconds
par.numslice = 30;
par.TA = par.TR-par.TR/par.numslice; %TA
par.sliceorder = 1:1:par.numslice; %slice order (assumes ascending)
par.refslice = floor(par.numslice/2); %reference slice (assumes middle)

par.maxvol = [216 216 216 216 216 216 146 146 146 146 146 146]; %highest volume number FOR EACH RUN
par.dropvol = 6; %dropped volumes (Assumed at start of scan)
par.minvol = par.dropvol+1; %first non-dropped volume
par.numvols = par.maxvol-par.dropvol; %number of volumes per scan
par.slicetiming(1) = par.TA / (par.numslice -1);%timing var for slice timing
par.slicetiming(2) = par.TR - par.TA; %timing var for slice timing

%% anatomy info
par.inplaneimg = [par.anatdir '/In001.img'];
par.hiresimg = [par.anatdir '/V001.img'];

%% variables for realigning/reslicing

% realign flags...
par.realflag.quality = 0.9;
par.realflag.fwhm = 5;
par.realflag.sep = 4;
par.realflag.rtm = 1;
par.realflag.interp = 2;

% reslice flags
par.reslflag.mask = 1;
par.reslflag.mean = 1;
par.reslflag.interp = 4;
par.reslflag.which = 0;

%% coregistration info
par.cor_meanfunc = [par.funcdir '/scan01/meanascan01.V007.img'];
par.cor_inimg = [par.anatdir '/In001.img'];
par.cor_hiresimg = [par.anatdir '/V001.img'];

%% segmenting info
par.img2bSeg = par.cor_hiresimg;
par.segopts = struct('biascor',1,'GM',[0 0 1],'WM',[0 0 1],'CSF',[0 0 0],'cleanup',0);
par.segbiasfwhm = 60; % 60 is default in gui, 75 is default in command line for reasons unexplained

%% normalization:

% gray matter:
par.graytemp = '/Applications/MATLAB_R2007b/spm5/apriori/grey.nii';
par.grsegs(1,:) = [par.anatdir '/c1' 'V001.img']; 
par.grsegs(2,:) = [par.anatdir '/c2' 'V001.img']; 
par.graysrcimg = [par.anatdir '/c1' 'V001.img']; 
par.graywrimg = [par.anatdir '/c1' 'V001.img']; 
par.grflags.smosrc = 8;
par.grflags.smoref = 0;
par.grflags.regtype = 'mni';
par.grflags.cutoff = 25;
par.grflags.nits = 16;
par.grflags.reg = 1;
par.grwrflags.preserve = 0;
par.grwrflags.bb = [[-78 -112 -50];[78 76 85]];
par.grwrflags.vox        = [2 2 2];
par.grwrflags.interp     = 1;
par.grwrflags.wrap       = [0 0 0];

% spgm normalization:
par.spgrnrmflags.preserve = 0;
par.spgrnrmflags.bb = [[-78 -112 -50];[78 76 85]];
par.spgrnrmflags.vox = [.5 .5 .5];
par.spgrnrmflags.interp     = 1;
par.spgrnrmflags.wrap       = [0 0 0];

%% smoothing funcs
par.smoothkernel = [6 6 6];

%% mask specification
par.specwrflags.preserve = 0;
par.specwrflags.bb = [[-78 -112 -50];[78 76 85]];
par.specwrflags.vox        = [1 1 1];
par.specwrflags.interp     = 1;
par.specwrflags.wrap       = [0 0 0];

par.specsmooth = [8 8 8];
par.maskimg = [par.anatdir '/mask.img']; 
par.smaskimg = [par.anatdir '/smask.img'];
par.tsmaskimg = [par.anatdir '/tsmask.img'];
par.wmaskimg = [par.anatdir '/wmask.img'];
par.swmaskimg = [par.anatdir '/swmask.img'];
par.tswmaskimg = [par.anatdir '/tswmask.img'];
par.twmaskimg = [par.anatdir '/twmask.img'];
par.addsegs = 'i1 + i2';
par.maskthresh = 'i1 > .2';


%% model variables
par.timing.fmri_t0 = par.refslice;  
par.timing.fmri_t = par.numslice; 
par.timing.units = 'secs';
par.bases.hrf.derivs = [0 0]; 
par.volt = 1;

par.sess.multi = {fullfile(par.behavdir, 'ons.mat')};
par.sess.multi_reg = {fullfile(par.behavdir, 'regs.mat')};
par.sess.hpf = 128;

par.cvi = 'AR(1)'; 
par.global = 'None';

% explicit mask
par.mask = {par.tswmaskimg};

% contrast variables
par.constat = 'T';

%% list and store all nifti files at various levels of processing
ac = 1;

for H = 1:par.numscans
    I = par.scans_to_include(H);
    
    scanNum = num2str(I);
    scanNumPrepended = prepend(scanNum,2);
    
    for J = par.minvol:par.maxvol(I)
        JFormatted = prepend(num2str(J),3);
        
        par.scanfiles{H}(J-par.dropvol,:) = [par.funcdir '/scan' scanNumPrepended '/scan' scanNumPrepended '.V' JFormatted '.img']; %scan files
        par.ascanfiles{H}(J-par.dropvol,:) = [par.funcdir '/scan' scanNumPrepended '/ascan' scanNumPrepended '.V' JFormatted '.img']; %ascan files
        par.wascanfiles{H}(J-par.dropvol,:) = [par.funcdir '/scan' scanNumPrepended '/wascan' scanNumPrepended '.V' JFormatted '.img']; %wascan files
        par.swascanfiles{H}(J-par.dropvol,:) = [par.funcdir '/scan' scanNumPrepended '/swascan' scanNumPrepended '.V' JFormatted '.img']; %swascan files
        
        par.allscanfiles(ac,:) = [par.funcdir '/scan' scanNumPrepended '/scan' scanNumPrepended '.V' JFormatted '.img']; %scan files
        par.allascanfiles(ac,:) = [par.funcdir '/scan' scanNumPrepended '/ascan' scanNumPrepended '.V' JFormatted '.img']; %ascan files
        par.allwascanfiles(ac,:) = [par.funcdir '/scan' scanNumPrepended '/wascan' scanNumPrepended '.V' JFormatted '.img']; %wascan files
        par.allswascanfiles(ac,:) = [par.funcdir '/scan' scanNumPrepended '/swascan' scanNumPrepended '.V' JFormatted '.img']; %swascan files
        ac = ac + 1;
    end
end

par.Enc.swafiles{1} = vertcat(par.swascanfiles{par.EncScans});
par.Ret.swafiles{1} = vertcat(par.swascanfiles{par.RetScans});

