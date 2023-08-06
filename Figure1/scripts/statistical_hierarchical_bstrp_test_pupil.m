


%%

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
%% select values for iterations
nAnimal_sel=10 % number of trials to be selected randomly
numSession=1 % number of sessions to be selected randomly
numSelections=5 % number of trials to be selected randomly
nruns = 100 % number of repeats 
ngroup=12 % 12 directions
%permuted_pupil_data_cell=cell(nAnimal_sel, nruns);

mat_sig_now=zeros(ngroup);
pval_sig_now=zeros(ngroup);

pairs=nchoosek(1:ngroup,2);
ks2stat_mean=nan(length(pairs),1);

for ipair=1:length(pairs)
    disp(['STARTS RUN=', sprintf('%d', ipair)]);
    
    pair_sel=pairs(ipair,:);
    
    permuted_pupil_data_cell={};
    ks2stat=nan(nruns,1);       
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

        nconditions=12;
        data_all={}; % cells containing the parameters of inidividual areas
        add_val=3;    

        for icond=1:nconditions
            data_all{icond}=rand_pupil_data(:,icond+add_val);
        end
        
        data_rand(:,1)=vector(data_all{pair_sel(1)});
        data_rand(:,2)=vector(data_all{pair_sel(2)});
        % run KS2 test
        [~,~,ks2stat(iruns)]=kstest2(data_rand(:,1),data_rand(:,2));
    
    end

    % ref: http://fcaglp.unlp.edu.ar/~observacional/papers/PDFs/statistics/Critical_KS.pdf
    % critical value for 2-sampled KS test
    mat_alpha_coefficient=[ 0.1 0.05 .025 .01 .005 .001;
        1.22 1.36 1.48 1.63 1.73 1.95];
    
    alpha_sel=2;
    coefficent_sel=mat_alpha_coefficient(2,alpha_sel);
    % critical value Da for large sample sizes.
    n1=numel(data_rand(:,1));
    n2=numel(data_rand(:,2));

    Da=coefficent_sel*sqrt((n1+n2)/(n1*n2));
    vec_Da=mat_alpha_coefficient(2,:).*sqrt((n1+n2)/(n1*n2));

    mean_now=mean(ks2stat);
    sem_now=std(ks2stat);
    lowCV=prctile(ks2stat,2.5);
    highCV=prctile(ks2stat,97.5);

    if lowCV>vec_Da(6)
        sig_level=4;
    elseif lowCV>vec_Da(4)
        sig_level=3;
    elseif lowCV>vec_Da(2)
        sig_level=2;
    else sig_level=1;
    end
    
    mat_sig_now(pair_sel(1),pair_sel(2))=sig_level*2;

end

%%
mat_sig_now(mat_sig_now==2)=1;
mat_sig_now=squareform(squareform(mat_sig_now'));
 
mat_sig_now(find(eye(ngroup)==1))=1;
mat_sig=mat_sig_now;

% PLOT RESULTS
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
order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

try_different_test_data=cell2mat(rand_pupil_data_all');

for_test_nasal_data=try_different_test_data(:,3+order_nasal);
for_test_temporal_data=try_different_test_data(:,3+order_temporal);

[pval, h0, stats]=signrank(for_test_nasal_data(:), for_test_temporal_data(:),'Tail','right');

% %% PAIRS AND BOOTSTRAPPING
% for ipair=1:size(pairs,1)
%     clear data1
%     clear data2
%     
%     pair_now=pairs(ipair,:);
%     
%     ks2stat=nan(nruns,1);
%     data_rand=[];
%     
%     data1=data_all{pair_now(1)};
%     data2=data_all{pair_now(2)};
% 
%     % FIRST PAIR    
%     num_lev1 = size(data1,1);
%     temp = NaN(nanimal,npopsize);
%     rand_lev1 = randi(num_lev1,nanimal,1);
%     num_lev2 = find(~isnan(data1(rand_lev1,:))); % remove NaNs
%     temp(iAn,1:length(num_lev2))=data1(rand_lev1(num_lev2),:);
% 
%     
%     data_rand(:,1)=vector(temp);
%     % SECOND PAIR  
%     num_lev1 = size(data2,1);
%     temp = NaN(nanimal,npopsize);
%     rand_lev1 = randi(num_lev1,nanimal,1);
%     num_lev2 = find(~isnan(data2(rand_lev1,:))); % remove NaNs
%     temp(iAn,1:length(num_lev2))=data2(rand_lev1(num_lev2),:);
% 
%     data_rand(:,2)=vector(temp);
% 
%     % KS TEST
%     [~,~,ks2stat(i)]=kstest2(data_rand(:,1),data_rand(:,2));
% 
% end
% 
% mat_alpha_coefficient=[ 0.1 0.05 .025 .01 .005 .001;
%     1.22 1.36 1.48 1.63 1.73 1.95];
% alpha_sel=2;
% coefficent_sel=mat_alpha_coefficient(2,alpha_sel);
% % critical value Da for large sample sizes.
% n1=numel(data_rand(:,1));
% n2=numel(data_rand(:,2));
% 
% Da=coefficent_sel*sqrt((n1+n2)/(n1*n2));
% vec_Da=mat_alpha_coefficient(2,:).*sqrt((n1+n2)/(n1*n2));
% 
% mean_now=mean(ks2stat);
% sem_now=std(ks2stat);
% lowCV=prctile(ks2stat,2.5);
% highCV=prctile(ks2stat,97.5);
% 
% %%