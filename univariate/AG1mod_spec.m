function AG1mod_spec(subpar,task)

% Set up the design matrix and run a design.
% written by jbh
% modified by amg


origdir = pwd;
% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = par_params(subpar);
end


spm_defaults;
global defaults
defaults.modality='FMRI';

% make results dir if not there:
if ~exist(par.(task).dir,'dir')
    mkdir(par.(task).dir);
end

cd(par.(task).dir);

%- rename files from previous analyses...
%-------------------------------------------------------------------
if exist(fullfile(par.(task).dir,'SPM.mat'),'file')
    movefile('SPM.mat',['SPM.mat-' date]);
end

% If we've gotten to this point we're committed to overwriting files.
% Delete them so we don't get stuck in spm_spm
%------------------------------------------------------------------------
files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
    '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
    '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$'};

for i=1:length(files)
    j = spm_select('List',pwd,files{i});
    for k=1:size(j,1)
        spm_unlink(deblank(j(k,:)));
    end
end

% Variables
%-------------------------------------------------------------
SPM.xY.RT = par.TR;
SPM.xY.P = [];

% Slice timing
defaults.stats.fmri.t=par.timing.fmri_t;
defaults.stats.fmri.t0=par.timing.fmri_t0;

% Basis function variables
%-------------------------------------------------------------
SPM.xBF.UNITS = par.timing.units;
SPM.xBF.dt    = par.TR/defaults.stats.fmri.t;
SPM.xBF.T     = defaults.stats.fmri.t;
SPM.xBF.T0    = defaults.stats.fmri.t0;

% Basis functions
%-------------------------------------------------------------
if strcmp(fieldnames(par.bases),'hrf')
    if all(par.bases.hrf.derivs == [0 0])
        SPM.xBF.name = 'hrf';
    elseif all(par.bases.hrf.derivs == [1 0])
        SPM.xBF.name = 'hrf (with time derivative)';
    elseif all(par.bases.hrf.derivs == [1 1])
        SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
    else
        error('Unrecognized hrf derivative choices.')
    end
else
    error('non-HRF selected! Not supported by this batchmode!');
end

SPM.xBF = spm_get_bf(SPM.xBF);
if isempty(par.sess),
    SPM.xBF.Volterra = false;
else
    SPM.xBF.Volterra = par.volt;
end;



for i = 1:numel(par.sess),
    if par.SPMSess
        sess = par.sess(i);
    else
        sess.multi = fullfile(par.(task).dir, 'ons.mat');
        sess.multi_reg = fullfile(par.(task).dir, 'regs.mat');
    end
    % Image filenames
    %-------------------------------------------------------------
    SPM.nscan(i) = size(par.(task).swafiles{i},1);
    SPM.xY.P     = strvcat(SPM.xY.P, par.(task).swafiles{i});
    U = [];

    % Augment the singly-specified conditions with the multiple
    % conditions specified in a .mat file provided by the user
    %------------------------------------------------------------
    if ~isempty(sess.multi)
        try
            multicond = load(sess.multi);
        catch
            error('Cannot load %s',sess.multi);
        end
        if ~(isfield(multicond,'names')&&isfield(multicond,'onsets')&&...
	     isfield(multicond,'durations')) || ...
		~all([numel(multicond.names),numel(multicond.onsets), ...
		      numel(multicond.durations)]==numel(multicond.names))
            error(['Multiple conditions MAT-file ''%s'' is invalid.\n',...
		   'File must contain names, onsets, and durations '...
		   'cell arrays of equal length.\n'],sess.multi{1});
        end
	
        %-contains three cell arrays: names, onsets and durations
        for j=1:length(multicond.onsets)
            cond.name     = multicond.names{j};
            cond.onset    = multicond.onsets{j};
            cond.duration = multicond.durations{j};
            
            
            % Mutiple Conditions Time Modulation
            %------------------------------------------------------
            % initialize the variable.
            cond.tmod = 0;
            if isfield(multicond,'tmod');
                try
                    cond.tmod = multicond.tmod{j};
                catch
                    error('Error specifying time modulation.');
                end
            end

            % Mutiple Conditions Parametric Modulation
            %------------------------------------------------------
            % initialize the parametric modulation variable.
            cond.pmod = [];
            if isfield(multicond,'pmod')
                % only access existing modulators
                try
                    % check if there is a parametric modulator. this allows
                    % pmod structures with fewer entries than conditions.
                    % then check whether any cells are filled in.
                    if (j <= numel(multicond.pmod)) && ...
                            ~isempty(multicond.pmod(j).name)

                        % we assume that the number of cells in each
                        % field of pmod is the same (or should be).
                        for ii = 1:numel(multicond.pmod(j).name)
                            cond.pmod(ii).name  = multicond.pmod(j).name{ii};
                            cond.pmod(ii).param = multicond.pmod(j).param{ii};
                            cond.pmod(ii).poly  = multicond.pmod(j).poly{ii};
                        end
                    end;
                catch
                    error('Error specifying parametric modulation.');
                end
            end
        sess.cond(j) = cond;
       end

    end

    % Configure the input structure array
    %-------------------------------------------------------------
    for j = 1:length(sess.cond),
        cond      = sess.cond(j);
        U(j).name = {cond.name};
        U(j).ons  = cond.onset(:);
        U(j).dur  = cond.duration(:);
        if length(U(j).dur) == 1
            U(j).dur    = U(j).dur*ones(size(U(j).ons));
        elseif length(U(j).dur) ~= length(U(j).ons)
            error('Mismatch between number of onset and number of durations.')
        end

        P  = [];
        q1 = 0;
        if cond.tmod>0,
            % time effects
            P(1).name = 'time';
            P(1).P    = U(j).ons*par.TR;
            P(1).h    = cond.tmod;
            q1        = 1;
        end;
        if ~isempty(cond.pmod)
            for q = 1:numel(cond.pmod),
                % Parametric effects
                q1 = q1 + 1;
                P(q1).name = cond.pmod(q).name;
                P(q1).P    = cond.pmod(q).param(:);
                P(q1).h    = cond.pmod(q).poly;
            end;
        end
        if isempty(P)
            P.name = 'none';
            P.h    = 0;
        end
        U(j).P = P;

    end

    SPM.Sess(i).U = U;


    % User specified regressors  ASSUMING NONE!!!!
    %-------------------------------------------------------------
    C = [];
    Cname = {}; 

    % Augment the singly-specified regressors with the multiple regressors
    % specified in the regressors.mat file
    %------------------------------------------------------------
    if ~strcmp(sess.multi_reg,'')
        tmp=load(char(sess.multi_reg));
        if isstruct(tmp) && isfield(tmp,'R')
            R = tmp.R;
        elseif isnumeric(tmp)
            % load from e.g. text file
            R = tmp;
        else
            warning('Can''t load user specified regressors in %s', ...
                char(sess.multi_reg{:}));
            R = [];
        end

        C=[C, R];
        nr=size(R,2);
        nq=length(Cname);
        for inr=1:nr,
            Cname{inr+nq}=['R',int2str(inr)];
        end
    end
    SPM.Sess(i).C.C    = C;
    SPM.Sess(i).C.name = Cname;

end

% Factorial design
%-------------------------------------------------------------
if isfield(par,'fact')
    if ~isempty(par.fact)
        NC=length(SPM.Sess(1).U); % Number of conditions
        CheckNC=1;
        for i=1:length(par.fact)
            SPM.factor(i).name=par.fact(i).name;
            SPM.factor(i).levels=par.fact(i).levels;
            CheckNC=CheckNC*SPM.factor(i).levels;
        end
        if ~(CheckNC==NC)
            disp('Error in fmri_spec job: factors do not match conditions');
            return
        end
    end
else
    SPM.factor=[];
end

% Globals
%-------------------------------------------------------------
SPM.xGX.iGXcalc = par.global;
SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca  = 'session specific';

% High Pass filter
%-------------------------------------------------------------
for i = 1:numel(par.sess),
    SPM.xX.K(i).HParam = par.sess(i).hpf;
end

% Autocorrelation
%-------------------------------------------------------------
SPM.xVi.form = par.cvi;

% Let SPM configure the design
%-------------------------------------------------------------
SPM = spm_fmri_spm_ui(SPM);

if ~isempty(par.mask)&&~isempty(par.mask{1})
    SPM.xM.VM         = spm_vol(par.mask{:});
    SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
end


%-Save SPM.mat
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')   %-#
if spm_matlab_version_chk('7') >= 0
    save('SPM','-V6','SPM'); %uncommented 071609
else
    save('SPM','SPM'); %uncommented 071609
end;

fprintf('%30s\n','...SPM.mat saved')                     %-#


cd(origdir); % Change back dir

%print out ps
par_print(par);


fprintf('Done\n')
return
%-------------------------------------------------------------------------