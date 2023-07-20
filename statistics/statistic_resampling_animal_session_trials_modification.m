
% select same number of animals
% select same number of sessions

clear all
% zscored values takes average zscored per session
%load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\pupil_data.mat','pupil_data')
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')
pupil_data= new_pupil_data

session_list=unique(pupil_data.sessions_id)
animals_id=unique(pupil_data.animal_id)

animal_id_list=pupil_data.animal_id
sessions_id_list=pupil_data.sessions_id

%data_selected_in=new_pupil_data.median_norm_diam_stims_offset;

%data_selected_in=new_pupil_data.raw_diam_stims_offset

data_selected_in=new_pupil_data.raw_relative_diam_delta;
%% select values for iterations
nAnimal_sel=5
nruns = 1000
numSelections=5

%permuted_pupil_data_cell=cell(nAnimal_sel, nruns);
permuted_pupil_data_cell={};

j=0;
for iruns=1:nruns;
    %step 1: select 7 animal's id randomly
    permuted_in = animals_id(randperm(numel(animals_id)));
    permuted_animal_id=permuted_in(1:nAnimal_sel);
    
    % step2: get one session per animal and randomly selected session if
    % there are multiple sessions
    for iAn=1:length(permuted_animal_id);
        % find sessions related to animal id; some animals have only 1 session
        j=j+1;
        clear animal_tmp
        clear index
        clear sessions_selected
        clear selected_trials
        clear permuted_selected_trials
        clear selected_indexes
        
        % select sessions from that animal
        animal_tmp=find(animal_id_list==permuted_animal_id(iAn));
        
        % randomly select one experimental session from that animal
        sessions_selected=sessions_id_list(animal_tmp);
        % chose sessions
        clear index
        index = randperm(length(animal_tmp));
        %session_chosen= animal_tmp(index(1))
        
        % SELECT ANIMAL ID
        selected_animal=permuted_animal_id(iAn);
        % select SESSION ID ONE SESSION
        selected_session=index(1);
        % select randomly selected 5 trials with a widthrown for randomly select animal and session
        selected_trials=intersect(find(data_selected_in(:,2)==index(1)), find(data_selected_in(:,1)==selected_animal));
        
        for i = 1:numSelections
            index_trial = randi(length(selected_trials));  % Generate a random index
            selected_indexes(i) = index_trial;  % Add the selected index to the selected_indexes array
        end
        
        permuted_selected_trials=selected_trials(selected_indexes);
        % create new cell with randomly chosen data
        permuted_pupil_data_cell{j}=data_selected_in(permuted_selected_trials,:);
        %permuted_pupil_data_cell{iAn, iruns}=pupil_data.zscored_data(permuted_selected_trials,:);

    end
    
end

permuted_pupil_data=cell2mat(permuted_pupil_data_cell');

'DONE'

% TEST DATA

%
nconditions=12

data_all={}; % cells containing the parameters of inidividual areas

%data_in=pupil_average_data % averaged over sessions
%add_val=2;

data_in=permuted_pupil_data; % all trials
add_val=3;    

for icond=1:nconditions
    data_all{icond}=data_in(:,icond+add_val);
end

% PLOT RESULTS

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
    %[H(ipair),P(ipair),CI{ipair}]=ttest2(data1(:),data2(:), 'Tail','right'); 
    [P(ipair), H(ipair), stat{ipair}]=signrank(data1(:),data2(:), 'Tail','right'); 
    %[P(ipair),H(ipair),STATS(ipair)] = ranksum(data1(:),data2(:));

    %[H(ipair),P(ipair),KSSTAT{ipair}]=kstest2(data1(:),data2(:));

    if P(ipair)>0.01
        sig_level=1
    elseif P(ipair)<=0.01 & P(ipair)>0.001
        sig_level=2
    elseif P(ipair)<=0.001 & P(ipair)>0.0001
        sig_level=3
    elseif P(ipair)<=0.0001 & P(ipair)>0.00001
        sig_level=4
    else sig_level=5;
        
    end
        pval_sig_now(pair_now(1),pair_now(2))=P(ipair);
        mat_sig_now(pair_now(1),pair_now(2))=sig_level*2;
    
    if P(ipair) <= alpha_Bonferroni
        Bonferroni_correction(ipair)=1;
        disp(['Iteration ', num2str(i), ': Reject null hypothesis']);
    else
        disp(['Iteration ', num2str(i), ': Fail to reject null hypothesis']);
        Bonferroni_correction(ipair)=0;
    end
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

%
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


set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];
print(gcf,'-dpdf',[filepathanalysis, 'Figure1D_statistical_test_12directions_btstrp.pdf']);


%% Bonferroni correction

    alpha = 0.05;
    alpha_Bonferroni = alpha / nruns;

    % Compare each p-value to the Bonferroni-corrected significance level
    for i = 1:nruns
        if P(ipair) <= alpha_Bonferroni
            Bonferroni_correction(i)=1;
            disp(['Iteration ', num2str(i), ': Reject null hypothesis']);
        else
            disp(['Iteration ', num2str(i), ': Fail to reject null hypothesis']);
            Bonferroni_correction(i)=0;
        end
    end





