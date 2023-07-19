%% code taken from Han, Bonin 2022 paper

clear all
% zscored values takes average zscored per session
%load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\pupil_data.mat','pupil_data')

% this load processed data for checking statistical tests
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat','new_pupil_data')
pupil_data=new_pupil_data

session_list=unique(pupil_data.sessions_id)
animals_id=unique(pupil_data.animal_id)


%%
nconditions=12

data_all={}; % cells containing the parameters of inidividual areas

%data_in=pupil_average_data % averaged over sessions
%add_val=2;
%data_in=pupil_data.zscored_data % all trials
data_in=pupil_data.median_norm_diam_stims_offset
add_val=3;    

for icond=1:nconditions
    data_all{icond}=data_in(:,icond+add_val);
end


%% Fig 3. example hierarchical KS test
threshold=0.05; % default significance level 95%

options.nruns = 1000; % # random sampling
options.nanimal = 5;  % # mice per sampling
options.npopsize = 150; % # neurons per mouse
options.isplot = 0; % plot intermediate steps or not

nruns = 1000; % # random sampling
nanimal = 5;  % # mice per sampling
npopsize = 150; % # neurons per mouse
isplot = 0; % plot intermediate steps or not

%%
ngroup=numel(data_all);
mat_sig_now=zeros(ngroup);

pairs=nchoosek(1:ngroup,2);

for ipair=1:size(pairs,1)
    clear data1
    clear data2
    
    pair_now=pairs(ipair,:);
    
    ttest2stat=nan(nruns,1);
    data_rand=[];
    
    data1=data_all{pair_now(1)};
    data2=data_all{pair_now(2)};
    
    for i =1:nruns
        
    num_lev1 = size(data1,1);
    temp = NaN(nanimal,npopsize);
    rand_lev1 = randi(num_lev1,nanimal,1);
    for j = 1:length(rand_lev1);
        num_lev2 = find(~isnan(data1(rand_lev1(j),:)),1,'last'); %We need to calculate this again here because there is a different number of trials for each neuron
        rand_lev2 = randi(num_lev2,1,npopsize); %Resample only from trials with data but same number of sample trials for all
        temp(j,:) = data1(rand_lev1(j),rand_lev2);
    end
    
    data_rand(:,1)=vector(temp);
    %

    num_lev1 = size(data2,1);
    temp = NaN(nanimal,npopsize);
    rand_lev1 = randi(num_lev1,nanimal,1);
    for j = 1:length(rand_lev1)
        num_lev2 = find(~isnan(data2(rand_lev1(j),:)),1,'last'); %We need to calculate this again here because there is a different number of trials for each neuron
        rand_lev2 = randi(num_lev2,1,npopsize); %Resample only from trials with data but same number of sample trials for all
        temp(j,:) = data2(rand_lev1(j),rand_lev2);
    end
    data_rand(:,2)=vector(temp);

    end
    
    % run ttest2 test
    %[~,~,ttest2stat(i)]=ttest2(data_rand(:,1),data_rand(:,2));
    
    [H(ipair),P(ipair),CI{ipair}]=ttest2(data_rand(:,1),data_rand(:,2));
    %[P(ipair),H(ipair),STATS(ipair)] = ranksum(data_rand(:,1),data_rand(:,2));

    if P(ipair)>0.05
        sig_level=1
    elseif P(ipair)<=0.05 & P(ipair)>0.01
        sig_level=2
    elseif P(ipair)<=0.01 & P(ipair)>0.001
        sig_level=3
    elseif P(ipair)<=0.001 & P(ipair)>0.0001
        sig_level=4
    else sig_level=5;
        
    end
    
        pval_sig_now(pair_now(1),pair_now(2))=P(ipair);
        mat_sig_now(pair_now(1),pair_now(2))=sig_level*2;

  
end

%%
%mat_sig_now(mat_sig_now==2)=1;
mat_sig_now=squareform(squareform(mat_sig_now'));
 
mat_sig_now(find(eye(ngroup)==1))=1;
mat_sig=mat_sig_now;

%%
figure
[XX,YY]=meshgrid([1:ngroup+1]-0.5,[1:ngroup+1]-0.5);
hp=pcolor(XX,YY,padarray(mat_sig,[1 1],0,'post'));
set(hp,'edgecolor','w');
colorspace=brewermap(8,'Blues');
colormap(colorspace)
set(gca,'clim',[0 8]);
axis square ij
box on
set(gca,'xtick',[],'ytick',[])

diams = 5; % coefficient of marker size
%Obtain the axes size (in axpos) in Points
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);

markersize = diams/diff(xlim)*min(axpos(3),axpos(4))*1.5; % Calculate Marker width in points
list_conditions={'0','30','60','90','120','150','180','210','240','270','300','330','control'}

set(gca,'xtick',[1:length(list_conditions)],'xticklabel',list_conditions);
set(gca,'ytick',[1:length(list_conditions)],'yticklabel',list_conditions);
xtickangle(90);
