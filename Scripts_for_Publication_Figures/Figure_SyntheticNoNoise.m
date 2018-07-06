% script to build Figure 3: changes between Sleep epochs
% general format:
% run plotting function
% tidy plot (ticks, ranges)
% export
% exportfig(h,[figpath 'Fig_' figID],'Color',color,'Format',format,'Resolution',dpi)


clear all; close all

% where are the data?
if ispc
    filepath = 'C:\Users\lpzmdh\Dropbox\Analyses\PfC Sampling hypothesis\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/PfC Sampling hypothesis/';    
end

% style sheet
run figure_properties

% which data for examples
type = 'Learn';
N = 35;
iBin = 4; % 5 ms
iSession = 7;

lW = 0.12;  % spacing of strip plots around the binsize tick mark

smallstrp = [10 15 4.5 3.5];

%% panel: words in common between sleep epochs


load([filepath 'DataWords_And_Counts_N' num2str(N) '_' type],'binsizes')
load([filepath 'UniqueWord_Data_N' num2str(N) '_Learn'])
LearnUWord = UWord; NLearn = numel(LearnUWord);
load([filepath 'UniqueWord_Data_N' num2str(N) '_Stable85'])
StableUWord = UWord; NStable = numel(StableUWord);


for iB = 1:numel(binsizes)
    for iL = 1:NLearn
        Learn.Pre.PSleep(iL,iB) = LearnUWord(iL).Bins(iB).Pre.PSleepEpochs;
        Learn.Post.PSleep(iL,iB) = LearnUWord(iL).Bins(iB).Post.PSleepEpochs;
    end
    for iL = 1:NStable
        Stable.Pre.PSleep(iL,iB) = StableUWord(iL).Bins(iB).Pre.PSleepEpochs;
        Stable.Post.PSleep(iL,iB) = StableUWord(iL).Bins(iB).Post.PSleepEpochs;
    end
end


figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',smallstrp);
for iB = 1:numel(binsizes)
    hold on
    plot(zeros(NStable,1)+binsizes(iB) + binsizes(iB)*lW,100*Stable.Pre.PSleep(:,iB),'o',...
            'MarkerFaceColor',colours.stable.marker,'MarkerEdgeColor',colours.stable.edge,'Markersize',M);
    plot(zeros(NLearn,1)+binsizes(iB) - binsizes(iB)*lW,100*Learn.Pre.PSleep(:,iB),'o',...
            'MarkerFaceColor',colours.learning.marker,'MarkerEdgeColor',colours.learning.edge,'Markersize',M);
    
end

set(gca,'XScale','log','XTick',xtick,'XTickLabel',strXlabel,'XLim',[xmin xmax],'XMinorTick','off');
set(gca,'YLim',[0 100])
xlabel('Bin size (ms)'); 
ylabel('Common words (%)'); 
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath 'Fig3_PreWordsInPost'],'-dsvg'); 

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',smallstrp);
for iB = 1:numel(binsizes)
    hold on
    plot(zeros(NStable,1)+binsizes(iB) + binsizes(iB)*lW,100*Stable.Post.PSleep(:,iB),'o',...
            'MarkerFaceColor',colours.stable.marker,'MarkerEdgeColor',colours.stable.edge,'Markersize',M);
    plot(zeros(NLearn,1)+binsizes(iB) - binsizes(iB)*lW,100*Learn.Post.PSleep(:,iB),'o',...
            'MarkerFaceColor',colours.learning.marker,'MarkerEdgeColor',colours.learning.edge,'Markersize',M);
    
end

set(gca,'XScale','log','XTick',xtick,'XTickLabel',strXlabel,'XLim',[xmin xmax],'XMinorTick','off');
set(gca,'YLim',[0 100])
xlabel('Bin size (ms)'); 
ylabel('Common words (%)'); 
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath 'Fig3_PostWordsInPre'],'-dsvg'); 

%% panel: proportion of time on common words between sleep
load([filepath 'DataWords_And_Counts_N' num2str(N) '_' type],'binsizes')
load([filepath 'PropTime_Data_N' num2str(N) '_Learn'])
LearnPTime = PTime; NLearn = numel(LearnPTime);
load([filepath 'PropTime_Data_N' num2str(N) '_Stable85'])
StablePTime = PTime; NStable = numel(StablePTime);


for iB = 1:numel(binsizes)
    for iL = 1:NLearn
        Learn.Pre.Common(iL,iB) = LearnPTime(iL).Bins(iB).Pre.SleepPtime;
        Learn.Post.Common(iL,iB) = LearnPTime(iL).Bins(iB).Post.SleepPtime;
    end
     for iS = 1:NStable
        Stable.Pre.Common(iS,iB) = StablePTime(iS).Bins(iB).Pre.SleepPtime;
        Stable.Post.Common(iS,iB) = StablePTime(iS).Bins(iB).Post.SleepPtime;
    end   
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',smallstrp);
for iB = 1:numel(binsizes)
    hold on
    plot(zeros(NStable,1)+binsizes(iB) + binsizes(iB)*lW,100*Stable.Pre.Common(:,iB),'o',...
            'MarkerFaceColor',colours.stable.marker,'MarkerEdgeColor',colours.stable.edge,'Markersize',M);

    plot(zeros(NLearn,1)+binsizes(iB) - binsizes(iB)*lW,100*Learn.Pre.Common(:,iB),'o',...
            'MarkerFaceColor',colours.learning.marker,'MarkerEdgeColor',colours.learning.edge,'Markersize',M);
    
end

set(gca,'XScale','log','XTick',xtick,'XTickLabel',strXlabel,'XLim',[xmin xmax],'XMinorTick','off');
set(gca,'YLim',[0 100])
xlabel('Bin size (ms)'); 
ylabel('Time by common words(%)'); 
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath 'Fig3_PreTimeCommon'],'-dsvg'); 

%% panel: example P(Pre) vs P(Post) scatter
load([filepath 'Pword_Data_N' num2str(N) '_' type], 'PWord');  % data Pword, and the binary IDs of each word

pPre = PWord(iSession).Bins(iBin).Pre.P; IDpre = PWord(iSession).Bins(iBin).Pre.binaryIDs;
pPost = PWord(iSession).Bins(iBin).Post.P; IDpost = PWord(iSession).Bins(iBin).Post.binaryIDs;
[~,cmnPre,cmnPost] = intersect(IDpre,IDpost);
plotScatter(pPre(cmnPre),pPost(cmnPost),[],figsize,colours.learning,widths,fontsize,fontname,'P(Pre)','P(Post)',1,3)
set(gca,'XScale','log','YScale','log','XLim',[1e-6 1],'yLim',[1e-6 1])
set(gca,'XTick',[1e-6 1e-4 1e-2 1],'YTick',[1e-6 1e-4 1e-2 1])
print([exportpath 'Fig3_PreVsPost_Scatter'],'-dsvg'); 

%% panel: scatter of D(Pre|Post) vs D(Pre*|Post*) [with CIs] for Learning - 5 ms

load([filepath 'DeltaSleep_Data_N' num2str(N) '_Learn'])
Nsessions = numel(DeltaSleep);

% Distance for the chosen example:
ExampleDistance = DeltaSleep(iSession).Bins(iBin).D_Pre_Post

% collate D(data), mean[D(permute)], and CI[D(permute)] over the sessions
Ddata = zeros(Nsessions,1); MDperm = zeros(Nsessions,1); CIperm = zeros(2,Nsessions);
for iS = 1:Nsessions
    Ddata(iS) = DeltaSleep(iS).Bins(iBin).D_Pre_Post;
    MDperm(iS) = DeltaSleep(iS).Bins(iBin).Perm.M_Pre_Post;
    ci99 = DeltaSleep(iS).Bins(iBin).Perm.CI_Pre_Post(2);
    CIperm(:,iS) = [MDperm(iS)-ci99; MDperm(iS)+ci99]';
end

plotScatter(Ddata,MDperm,CIperm,figsize,colours.learning,widths,fontsize,fontname,'Data: D(Pre|Post)','Null model: D(Pre*|Post*)',1,M)
print([exportpath 'Fig3_Learn_SleepChangeDiff_Null_Model'],'-dsvg');     


%% panel: scatter of D(Pre|Post) vs D(Pre*|Post*) [with CIs] for Stable -  5 ms

load([filepath 'DeltaSleep_Data_N' num2str(N) '_Stable85'])
Nsessions = numel(DeltaSleep);

% collate D(data), mean[D(permute)], and CI[D(permute)] over the sessions
Ddata = zeros(Nsessions,1); MDperm = zeros(Nsessions,1); CIperm = zeros(2,Nsessions);
for iS = 1:Nsessions
    Ddata(iS) = DeltaSleep(iS).Bins(iBin).D_Pre_Post;
    MDperm(iS) = DeltaSleep(iS).Bins(iBin).Perm.M_Pre_Post;
    ci99 = DeltaSleep(iS).Bins(iBin).Perm.CI_Pre_Post(2);
    CIperm(:,iS) = [MDperm(iS)-ci99; MDperm(iS)+ci99]';
end

plotScatter(Ddata,MDperm,CIperm,figsize,colours.stable,widths,fontsize,fontname,'Data: D(Pre|Post)','Null model: D(Pre*|Post*)',1,M)
print([exportpath 'Fig3_Stable_SleepChangeDiff_Null_Model'],'-dsvg');     


%% panel: summary over binsizes
lW = 0.075;  % spacing of strip plots around the binsize tick mark

load([filepath 'DataWords_And_Counts_N' num2str(N) '_' type],'binsizes')
load([filepath 'DeltaSleep_Data_N' num2str(N) '_Learn'])
LearnDelta = DeltaSleep; NLearn = numel(LearnDelta);
load([filepath 'DeltaSleep_Data_N' num2str(N) '_Stable85'])
StableDelta = DeltaSleep; NStable = numel(StableDelta);
load([filepath 'DataWords_And_Counts_N' num2str(N) '_Learn'],'binsizes')

for iB = 1:numel(binsizes)
    Diff(iB).Learn = zeros(NLearn,1);
    for iS = 1:NLearn
        Diff(iB).Learn(iS) = LearnDelta(iS).Bins(iB).D_Pre_Post - LearnDelta(iS).Bins(iB).Perm.M_Pre_Post;    
    end
    
    Diff(iB).Stable = zeros(NStable,1);
    for iS = 1:NStable
        Diff(iB).Stable(iS) = StableDelta(iS).Bins(iB).D_Pre_Post - StableDelta(iS).Bins(iB).Perm.M_Pre_Post;            
    end
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 7 4]);
line([eps xmax],[0 0],'Color',colours.shuf.line,'Linewidth',widths.axis); hold on
for iB = 1:numel(binsizes)
    hold on
    plot(zeros(NStable,1)+binsizes(iB) + binsizes(iB)*lW,Diff(iB).Stable,'o',...
            'MarkerFaceColor',colours.stable.marker,'MarkerEdgeColor',colours.stable.edge,'Markersize',M);
    plot(zeros(NLearn,1)+binsizes(iB) - binsizes(iB)*lW,Diff(iB).Learn,'o',...
            'MarkerFaceColor',colours.learning.marker,'MarkerEdgeColor',colours.learning.edge,'Markersize',M);
    
end
set(gca,'XScale','log','XTick',xtick,'XTickLabel',strXlabel,'XLim',[xmin xmax],'XMinorTick','off');
set(gca,'YLim',[-0.1 0.5])
xlabel('Bin size (ms)'); 
ylabel('Difference: D(data) - D(null)'); 
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath 'Fig3_AllBinsizes_SleepChangeDiff_Null_Model'],'-dsvg');     



