
function AG1wholeshebang(subpar, flags)

% runs preprocessing and univariate analyses for a subject
% written by jbh
% modified by amg
%
% flags defaults to all stages, or individual stages can be chosen from:
% m = makevols
% a = make anatomicals
% s = slice timing
% l = realignment
% c = coregister inplane anat to mean func, and hires to inplane
% g = segment and normalize gray matter
% n = normalize functionals
% h = smooth functionals
% w = normalize anatomy
% k = make specmask
% z = make artifact indices
% f = make mean functinoal
% d = make regressors
% p = specify model
% e = run model
% t = make contrasts
%

origdir = pwd;

if ~exist('subpar', 'var')
    error('Must specify subject!');
end

par = subpar;


% Preprocessing
if ismember('m',flags); par_makevols(par); end

if ismember('a',flags); par_makemoveanat(par); end

if ismember('s',flags); par_slicetime(par); end

if ismember('l',flags); par_realign(par); end

if ismember('c',flags); par_coregwrapper(par); end

if ismember('g',flags); par_segnorm(par); end

if ismember('n',flags); par_normfuncs(par); end

if ismember('h',flags); par_smoothfuncs(par); end

if ismember('w',flags); AG1_normSPGR(par); end

if ismember('k',flags); par_makespec(par); end

if ismember('z',flags); AG1_ArtScansToOns(par); end

if ismember('f',flags); AG1_meanFuncs(par); end

% Modeling, etc
    if ismember('d',flags); AG1_make_regs(par, 1); end
for t = 1:length(par.Tasks)
    if ismember('p',flags); AG1_mod_spec(par, par.Tasks{t}); end
    if ismember('e',flags); AG1_mod_est(par, par.Tasks{t}); end
    if ismember('t',flags); AG1_setcontrasts(par, par.Tasks{t}); end
end

cd(origdir);
