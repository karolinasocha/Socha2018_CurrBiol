
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

data_selected_in=new_pupil_data.raw_relative_diam_delta;

% general test to find siginifcance in nasal vs temporal directions
%% select values for iterations
nAnimal_sel=10 % number of trials to be selected randomly
numSession=1 % number of sessions to be selected randomly
numSelections=7 % number of trials to be selected randomly
nruns = 1000 % number of repeats 
ngroup=2 % 2 directions : nasal vs temporal

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

permuted_pupil_data_cell={};

for iruns=1:nruns;

    data_rand=[];

    %step 1: select animal's id randomly
    permuted_in = animals_id(randperm(numel(animals_id)));
    permuted_animal_id=permuted_in(1:nAnimal_sel); % selected 5 animals randomly

    % step2: get one session per animal and randomly selected session if
    % there are multiple sessions

    for iAn=1:length(permuted_animal_id);
        % find sessions related to animal id; some animals have only 1 session

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
        index = randperm(length(sessions_selected));

        % SELECT ANIMAL ID
        selected_animal=permuted_animal_id(iAn);
        % select SESSION ID ONE SESSION
        selected_session=index(numSession);
        % select trials from randomly chosen animal and session
        selected_trials=intersect(find(data_selected_in(:,2)==selected_session), find(data_selected_in(:,1)==selected_animal));

        % select randomly 5 (numSelections) trials with a widthrown from randomly selected animal and session
        rand_selected_trials=selected_trials(randi(length(selected_trials),1,numSelections));

        % create new cell with randomly chosen data from one animal
        rand_pupil_data_cell{iAn}=data_selected_in(rand_selected_trials,:);
    end

    % randomly selected data
    rand_pupil_data=cell2mat(rand_pupil_data_cell');
    rand_pupil_data_all{iruns}=rand_pupil_data;

    for_test_nasal_data=rand_pupil_data(:,3+order_nasal);
    for_test_temporal_data=rand_pupil_data(:,3+order_temporal);

    [pval(iruns), h0(iruns), stats(iruns)]=signrank(100*for_test_nasal_data(:), 100*for_test_temporal_data(:),'Tail','right');

end

% Original significance level
alpha = 0.05;
alpha_Bonferroni = alpha / nruns;

% Compare each p-value to the Bonferroni-corrected significance level
for i = 1:nruns
    if pval(i) <= alpha_Bonferroni
        Bonferroni_correction(i)=1;
        disp(['Iteration ', num2str(i), ': Reject null hypothesis']);
    else
        disp(['Iteration ', num2str(i), ': Fail to reject null hypothesis']);
        Bonferroni_correction(i)=0;
    end
end

sum(Bonferroni_correction)

%% TRY TO COMPARE EFFECT ACROSS DIFFERENT ANIMALS

%% compare multiple animals
data_multiple_animals=cell2mat(rand_pupil_data_all');
ngroup=length(animals_id);
mat_sig_now=zeros(ngroup);

pairs=nchoosek(1:ngroup,2);

for ipair=1:length(pairs)
    pair_now=pairs(ipair,:);
    animal_compare=animals_id(pair_now);

    data_animal1=data_multiple_animals(find(data_multiple_animals(:,1)==animal_compare(1)),3+order_nasal);
    data_animal2=data_multiple_animals(find(data_multiple_animals(:,1)==animal_compare(2)),3+order_nasal);
    
    [pval_animal(ipair), h0_animal(ipair), stats_animal(ipair)]=ranksum(100*data_animal1(:), 100*data_animal2(:));

    if pval_animal(ipair)>0.05
        sig_level=1
    elseif pval_animal(ipair)<=0.05 & pval_animal(ipair)>0.01
        sig_level=2
    elseif pval_animal(ipair)<=0.01 pval_animal(ipair)>0.001
        sig_level=3
    elseif pval_animal(ipair)<=0.001 & pval_animal(ipair)>0.0001
        sig_level=4
    else sig_level=5;
        
    end
    
        pval_sig_now(pair_now(1),pair_now(2))=pval_animal(ipair);
        mat_sig_now(pair_now(1),pair_now(2))=sig_level*2;
    
end

%%
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
% 
% order_nasal=[11,12,1,2,3,4];
% order_temporal=[5,6,7,8,9,10];
% 
% try_different_test_data=cell2mat(rand_pupil_data_all');
% 
% for_test_nasal_data=try_different_test_data(:,3+order_nasal);
% for_test_temporal_data=try_different_test_data(:,3+order_temporal);
% 
% [pval_all, h0_all, stats_all]=signrank(100*for_test_nasal_data(:), 100*for_test_temporal_data(:),'Tail','right');
% 
%     

