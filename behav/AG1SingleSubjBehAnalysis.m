function [dat idx] = AG1SingleSubjBehAnalysis(subject)
% analyzes behavioral results of AG1 experiment
%
% wrtten by amg

%% setup

% directories of interest
expDir = '/Users/gordonam/Studies/AG1/fmri_data2';
myDir.subj = fullfile(expDir, subject);
myDir.Beh= (fullfile(expDir, subject, 'behav'));
cd (myDir.Beh);

%data file names for the tasks
encDirH = dir ('AG1_enc*');
retDirH = dir ('AG1_ret*');

%load the raw data files
rawDat.Enc = load(encDirH.name);
rawDat.Ret = load(retDirH.name);

% split data into tasks
dat.sub.rawDat = rawDat;
dat.thisHand.Enc = rawDat.Enc.S.encHandNum;
dat.thisHand.Ret = rawDat.Ret.S.retHandNum;


%% Encoding trial indices

% mapping hands to responses
HandNum.enc = rawDat.Enc.S.encHandNum;
Hands.enc = {'1!' '5%'};
thisHand.enc = Hands.enc{HandNum.enc};
thatHand.enc = Hands.enc{3-HandNum.enc};

% raw stored behavioral variables
enc.Cond = vertcat(rawDat.Enc.EncData(:).cond) ;
enc.Resp = horzcat(rawDat.Enc.EncData(:).conf)';
enc.RT = cell2mat(horzcat(rawDat.Enc.EncData(:).confRT))';
enc.Item = vertcat(rawDat.Enc.EncData(:).item);
enc.Onset = horzcat(rawDat.Enc.EncData(:).onset);

% session index
for j = 1:length(rawDat.Enc.EncData)
    sess_h{j} = j*ones(size(rawDat.Enc.EncData(j).cond));
end
enc.Sess = vertcat(sess_h{:});

%index of person or scene presentations
idx.enc.person = (enc.Cond==1);
idx.enc.scene = (enc.Cond==2);

%index of confident, nonconfident, other, and clean responses
idx.enc.conf = strcmp(enc.Resp, thisHand.enc);
idx.enc.nconf = strcmp(enc.Resp, thatHand.enc);
idx.enc.other = ones(180,1)-idx.enc.conf - idx.enc.nconf;
idx.enc.cleanResp = idx.enc.conf + idx.enc.nconf;

% index of confident/nonconfident X person/scene responses
idx.enc.personConf = idx.enc.person.*idx.enc.conf;
idx.enc.sceneConf = idx.enc.scene.*idx.enc.conf;
idx.enc.personNConf = idx.enc.person.*idx.enc.nconf;
idx.enc.sceneNConf = idx.enc.scene.*idx.enc.nconf;

% onsets and session
idx.enc.onsets = enc.Onset';
idx.enc.sessToSPM = enc.Sess;


%% Encoding results

%percent confident for person and scene trials
dat.sub.Enc.pctConf.person = sum(idx.enc.personConf)/sum(idx.enc.person.*idx.enc.cleanResp);
dat.sub.Enc.pctConf.scene = sum(idx.enc.sceneConf)/sum(idx.enc.scene.*idx.enc.cleanResp);

%RT for person and scene trials
dat.sub.Enc.RT.person = median(enc.RT(find(idx.enc.person.*idx.enc.conf)));
dat.sub.Enc.RT.scene = median(enc.RT(find(idx.enc.scene.*idx.enc.conf)));


%% Ret trial indices

% mapping hands to responses
HandNum.ret = rawDat.Ret.S.retHandNum;
Hands.ret = {'1!' '5%'};
thisHand.ret = Hands.ret(HandNum.ret);
thatHand.ret = Hands.ret(3-HandNum.ret);

% words
ret.item = vertcat(rawDat.Ret.retData(:).item);



% retrieval responses and reaction times (use the first recorded response)
ret.resp = {};
ret.respRT = [];
ret.sess = [];
for i = 1:length(rawDat.Ret.retData)
    for j = 1:length(rawDat.Ret.retData(i).stimRT)
        ret.respRT = [ret.respRT; rawDat.Ret.retData(i).stimRT{j}(1)];
        
        if iscell(rawDat.Ret.retData(i).stimresp{j})
            ret.resp = [ret.resp; rawDat.Ret.retData(i).stimresp{j}{1}];
        else
            ret.resp = [ret.resp; rawDat.Ret.retData(i).stimresp{j}];
        end
    end
    ret.sess = [ret.sess; i*ones(size(rawDat.Ret.retData(i).stimRT))'];
end


%When making conditions here, refer to the condition used in the encoding
%list. 
for c = 1:length(ret.item)
    EncIdx = strmatch(ret.item{c}, enc.Item);
    ret.cond(c) = enc.Cond(EncIdx);
    ret.conf(c) = enc.Resp(EncIdx);
end

%actual person/scene condition
idx.ret.person =  (ret.cond==1)';
idx.ret.scene = (ret.cond==2)';

%response was "person" or "scene"
idx.ret.respPerson = strcmp(ret.resp, thisHand.ret);
idx.ret.respScene = strcmp(ret.resp, thatHand.ret);

%was the subject confident of their encoding of this item
idx.ret.conf = strcmp(ret.conf, thisHand.enc)';

%sensical responses, i.e. only "scene" or "person," no missed responses
%or other button presses
idx.ret.cleanResp = idx.ret.respPerson + idx.ret.respScene;

% indices of confident/nonconfident, correct/incorrect, person/scene
% responses
idx.ret.respPersonCorConf = (idx.ret.person.*idx.ret.respPerson.*idx.ret.conf);
idx.ret.respSceneCorConf = (idx.ret.scene.*idx.ret.respScene.*idx.ret.conf);
idx.ret.respPersonCorNonConf = (idx.ret.person.*idx.ret.respPerson.*~idx.ret.conf);
idx.ret.respSceneCorNonConf = (idx.ret.scene.*idx.ret.respScene.*~idx.ret.conf);
idx.ret.respPersonIncConf = (idx.ret.person.*idx.ret.respScene.*idx.ret.conf);
idx.ret.respSceneIncConf = (idx.ret.scene.*idx.ret.respPerson.*idx.ret.conf);

idx.ret.cor = (idx.ret.person.*idx.ret.respPerson + idx.ret.scene.*idx.ret.respScene);
idx.ret.inc = (idx.ret.person.*idx.ret.respScene + idx.ret.scene.*idx.ret.respPerson);

idx.ret.confCor = idx.ret.respSceneCorConf+idx.ret.respPersonCorConf;
idx.ret.confInc = idx.ret.respSceneIncConf+idx.ret.respPersonIncConf;

idx.ret.sess = ret.sess;

%% Subsequent memory indices
for c = 1:length(enc.Item)
    RetIdx = strmatch(enc.Item{c}, ret.item);
    if isempty(RetIdx)
        subsMem(c) = -1;
    elseif idx.ret.cor(RetIdx)==1;
        subsMem(c) = 1;
    elseif idx.ret.inc(RetIdx)==1;
        subsMem(c) = 0;
    else
        subsMem(c) = -1;
    end
end

idx.enc.subsMem = (subsMem==1)';
idx.enc.subsForgotten = (subsMem==0)';

%% Retrieval results

%percent correct for all person and scene trials
dat.sub.Ret.pctCor.person = sum(idx.ret.person.*idx.ret.respPerson)/sum(idx.ret.person.*idx.ret.cleanResp);
dat.sub.Ret.pctCor.scene = sum(idx.ret.scene.*idx.ret.respScene)/sum(idx.ret.scene.*idx.ret.cleanResp);

%percent correct for all person and scene trials that were rated
%confidently in the encoding phase
dat.sub.Ret.pctCor.personConf = sum(idx.ret.respPersonCorConf)/sum(idx.ret.person.*idx.ret.cleanResp.*idx.ret.conf);
dat.sub.Ret.pctCor.sceneConf = sum(idx.ret.respSceneCorConf)/sum(idx.ret.scene.*idx.ret.cleanResp.*idx.ret.conf);
dat.sub.Ret.pctCor.personNonConf = sum(idx.ret.respPersonCorNonConf)/sum(idx.ret.person.*idx.ret.cleanResp.*~idx.ret.conf);
dat.sub.Ret.pctCor.sceneNonConf = sum(idx.ret.respSceneCorNonConf)/sum(idx.ret.scene.*idx.ret.cleanResp.*~idx.ret.conf);
dat.sub.Ret.pctCor.allConf = sum(idx.ret.respPersonCorConf+idx.ret.respSceneCorConf)/sum(idx.ret.cleanResp.*idx.ret.conf);

% breakdown of accuracy by session
for tl = 1:size(rawDat.Ret.retData,2)
    iS = [(tl*28-27):(tl*28)];
    dat.sub.Ret.pctCor.bySessPersonConf(tl) = sum(idx.ret.respPersonCorConf(iS))/sum(idx.ret.person(iS).*idx.ret.cleanResp(iS).*idx.ret.conf(iS));
    dat.sub.Ret.pctCor.bySessSceneConf (tl) =  sum(idx.ret.respSceneCorConf(iS))/sum(idx.ret.scene(iS).*idx.ret.cleanResp(iS).*idx.ret.conf(iS));
end

%allRTs
dat.sub.Ret.RT.allRTs = ret.respRT;

%RTs for person correct, person all, and person confident correct
dat.sub.Ret.RT.personCor = median(ret.respRT(find(idx.ret.person.*idx.ret.respPerson)));
dat.sub.Ret.RT.personAll = median(ret.respRT(find(idx.ret.person)));
dat.sub.Ret.RT.personCorConf = median(ret.respRT(find(idx.ret.person.*idx.ret.respPerson.*idx.ret.conf)));

%RTs for scene correct, scene all, and scene confident correct
dat.sub.Ret.RT.sceneCor = median(ret.respRT(find(idx.ret.scene.*idx.ret.respScene)));
dat.sub.Ret.RT.sceneAll = median(ret.respRT(find(idx.ret.scene)));
dat.sub.Ret.RT.sceneCorConf = median(ret.respRT(find(idx.ret.scene.*idx.ret.respScene.*idx.ret.conf)));

dat.sub.Ret.RT.allConfByTrial = ret.respRT(find(idx.ret.confCor));
dat.sub.Ret.RT.CorConf = median(ret.respRT(find(idx.ret.person.*idx.ret.respPerson.*idx.ret.conf + idx.ret.scene.*idx.ret.respScene.*idx.ret.conf)));

dat.sub.Ret.RT.personIncConf = median(ret.respRT(find(idx.ret.person.*idx.ret.respScene.*idx.ret.conf)));
dat.sub.Ret.RT.sceneIncConf = median(ret.respRT(find(idx.ret.scene.*idx.ret.respPerson.*idx.ret.conf)));
dat.sub.Ret.RT.IncConf =  median(ret.respRT(find(idx.ret.person.*idx.ret.respScene.*idx.ret.conf + idx.ret.scene.*idx.ret.respPerson.*idx.ret.conf)));


