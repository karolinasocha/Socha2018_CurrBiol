
% select same number of animals
% select same number of sessions

clear all
% zscored values takes average zscored per session
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\pupil_data.mat','pupil_data')

session_list=unique(pupil_data.sessions_id)
animals_id=unique(pupil_data.animal_id)

animal_id_list=pupil_data.animal_id
sessions_id_list=pupil_data.sessions_id
%% 
nAnimal_sel=5
nruns = 100
numSelections=3

%permuted_pupil_data_cell=cell(nAnimal_sel, nruns);

i=0;

for iruns=1:nruns;
    i=i+1;
    
    permuted_pupil_data_cell={};
    
    %step 1: select 7 animal's id randomly
    permuted_in = animals_id(randperm(numel(animals_id)));
    permuted_animal_id=permuted_in(1:nAnimal_sel);
    
    % step2: get one session per animal and randomly selected session if
    % there are multiple sessions
    
    j=0;
    for iAn=1:length(permuted_animal_id);
        % find sessions related to animal id; some animals have only 1 session
        j=j+1;
        clear animal_tmp
        clear index
        clear sessions_selected
        clear selected_trials
        clear permuted_selected_trials
        clear selected_indexes
        
        animal_tmp=find(animal_id_list==permuted_animal_id(iAn));
        
        % randomly select one experimental session from that animal
        sessions_selected=sessions_id_list(animal_tmp);
        index = randperm(length(animal_tmp));
        %session_chosen= animal_tmp(index(1))
        
        selected_animal=permuted_animal_id(iAn);
        selected_session=index(1);
        % select randomly selected 5 trials with a widthrown for randomly select animal and session
        selected_trials=intersect(find(pupil_data.zscored_data(:,2)==index(1)), find(pupil_data.zscored_data(:,1)==selected_animal));
        
        for i = 1:numSelections
            index_trial = randi(length(selected_trials));  % Generate a random index
            selected_indexes(i) = index_trial;  % Add the selected index to the selected_indexes array
        end
        
        permuted_selected_trials=selected_trials(selected_indexes);
        % create new cell with randomly chosen data
        permuted_pupil_data_cell{j}=pupil_data.zscored_data(permuted_selected_trials,:);
        %permuted_pupil_data_cell{iAn, iruns}=pupil_data.zscored_data(permuted_selected_trials,:);

    end

    permuted_pupil_data=cell2mat(permuted_pupil_data_cell');
    
    nconditions=12
    data_all={}; % cells containing the parameters of inidividual areas
    %data_in=pupil_average_data % averaged over sessions
    %add_val=2;

    data_in=permuted_pupil_data; % all trials
    add_val=3;    

    for icond=1:nconditions
        data_all{icond}=data_in(:,icond+add_val);
    end
    
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
        [P(ipair),H(ipair),STATS(ipair)] = ranksum(data1(:),data2(:));

        %[H(ipair),P(ipair),KSSTAT{ipair}]=kstest2(data1(:),data2(:));


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
    
    pvals_test{iruns}=pval_sig_now;
    mat_sign{iruns}=mat_sig_now;



end
%%
tmp=zeros(12,12)
for i=1:100
tmp=pvals_test{i}+tmp;
end

tmp2=tmp/100
%%
ngroup=numel(data_all);
mat_sig=zeros(ngroup);
    
for ipair=1:size(pairs,1)

    pair_now=pairs(ipair,:);  

    P(ipair)=tmp2(pair_now(1),pair_now(2))
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

    mat_sig(pair_now(1),pair_now(2))=sig_level*2;
end

%mat_sig_now(mat_sig_now==2)=1;
mat_sig_now2=squareform(squareform(mat_sig'));
 
mat_sig_now2(find(eye(ngroup)==1))=1;
mat_sig_final=mat_sig_now2;

%
figure
[XX,YY]=meshgrid([1:ngroup+1]-0.5,[1:ngroup+1]-0.5);
hp=pcolor(XX,YY,padarray(mat_sig_final,[1 1],0,'post'));
set(hp,'edgecolor','w');
colorspace=brewermap(8,'Blues');
colormap(colorspace)
set(gca,'clim',[0 8]);
axis square ij
box on
set(gca,'xtick',[],'ytick',[])
%
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

%%