%% without bootstrapping only ttest
% average zscored pupil size for each session (in total 40)

%% code taken from Han, Bonin 2022 paper

clear all
% zscored values takes average zscored per session
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\pupil_data.mat','pupil_data')

session_list=unique(pupil_data.sessions_id)
animals_id=unique(pupil_data.animal_id)

pupil_average_data=nan(length(pupil_data.sessions_id),15);
j=0
for iAn=1:length(animals_id)

    indecies_animal=find(pupil_data.zscored_data(:,1)==animals_id(iAn));
    indecies_session=unique(pupil_data.zscored_data(indecies_animal,2));
    
for isess=1:length(indecies_session)
    j=j+1
    indecies=intersect(find(pupil_data.zscored_data(:,2)==indecies_session(isess)),indecies_animal);
    tmp=pupil_data.zscored_data(indecies,:);
    pupil_average_data(j,:)=[animals_id(iAn),session_list(isess), nanmean(tmp(:,4:16),1)];
end

end

%%
nconditions=13

data_all={}; % cells containing the parameters of inidividual areas

%data_in=pupil_average_data % averaged over sessions
%add_val=2;
data_in=pupil_data.zscored_data % all trials
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


ngroup=numel(data_all);
mat_sig_now=zeros(ngroup);
pval_sig_now=zeros(ngroup);

pairs=nchoosek(1:ngroup,2);
ks2stat=nan(length(pairs),1);
    
for ipair=1:size(pairs,1)
    clear data1
    clear data2
    pair_now=pairs(ipair,:);    
    data1=data_all{pair_now(1)};
    data2=data_all{pair_now(2)};
    % run KS2 test
    %[~,~,ks2stat(ipair)]=kstest2(data1(:),data2(:)); 
    %[H(ipair),P(ipair),CI{ipair}]=ttest2(data1(:),data2(:)); 
    %[P(ipair),H(ipair),STATS(ipair)] = ranksum(data1(:),data2(:));
    
    [H(ipair),P(ipair),KSSTAT{ipair}]=kstest2(data1(:),data2(:));
   
    
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

%
%mat_sig_now(mat_sig_now==2)=1;
pval_sig_now2=squareform(squareform(pval_sig_now'));
 
pval_sig_now2(find(eye(ngroup)==1))=1;
pval_sig=pval_sig_now2;

%mat_sig_now(mat_sig_now==2)=1;
mat_sig_now2=squareform(squareform(mat_sig_now'));
 
mat_sig_now2(find(eye(ngroup)==1))=1;
mat_sig=mat_sig_now2;

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
%%
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

%     hold on
%     colors_area=getColorAreas;
%     for iarea=1:ngroup
%         scatter(iarea,iarea,markersize,'markerfacecolor',colors_area(list_areas{iarea),:),...
%             'markeredgecolor','w');
%     end

%%
% data1=[data_all{[1,2,3,4,11,12]}];
% data2=[data_all{[5,6,7,8,9,10]}];
% 
% [H,P,CI]=ttest2(data1(:),data2(:));
% pval=ranksum(data1(:),data2(:))
% 
% [H_ks,P_ks, KSstats] = kstest2(data1(:),data2(:))