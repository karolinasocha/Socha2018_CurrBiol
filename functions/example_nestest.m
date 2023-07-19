%% transform the dataset to a hierarchially organization
vec_para=rand(1000,1);
vec_resp=rand(1000,1)>0.5;
vec_area=repmat(1:5,100,1); vec_area=vec_area(:);vec_area=cat(1,vec_area,vec_area);
vec_animal=repmat(1:2,500,1);vec_animal=vec_animal(:);
data_all=prepareHierData(vec_para,vec_resp,vec_area,vec_animal);

%% generate hierarchial bootstrap dataset
nrun=1000;  % number of bootstrap
nanimal=5;  % number of animal per bootstrap
npopsize=50; % number of cells per animal
nrandsample=nanimal*npopsize; % number of cells per bootstrap
data_hier_all={};
tic
for igroup=1:numel(data_all)
    
    data_now=data_all{igroup};
    mat_now=[];
    for irun =1:nrun
        seed_selection_now=irun;
        rng(seed_selection_now);
        
        num_lev1 = size(data_now,1);
        temp = NaN(nanimal,npopsize); % new bootstrap dataset
        rand_lev1 = randi(num_lev1,nanimal,1);
        for j = 1:length(rand_lev1)
            num_lev2 = find(~isnan(data_now(rand_lev1(j),:)),1,'last'); %We need to calculate this again here because there is a different number of trials for each neuron
            rand_lev2 = randi(num_lev2,1,npopsize); %Resample only from trials with data but same number of sample trials for all
            temp(j,:) = data_now(rand_lev1(j),rand_lev2);
        end
        
        idx_rand=temp(:);
        mat_now(:,irun)=idx_rand;
    end
    
    data_hier_all{igroup}=mat_now;
end

toc
%% you can compare the distribution of bootstrap datasets as you wish 


%% nested function
function data_all=prepareHierData(vec_para,vec_resp,vec_area,vec_animal)
% Input variables:
% vec_para: a vector of parameter values
% vec_resp: a vector of bollean values showing whether the cell is responsive (1,0,0,1,1,...)
% vec_area: a vector of cells' areal id (area 1, 2, 3,...)
% vec_animal: a vector of cell's animal id (mouse 1,2,3,...)
% output: a cell variable where each cell is a MxN matrix describing each area's neurons activity, where M is the number of animal, N is the maximal number of neurons across the M animals. The rest of the matrix is padded with NaN.
%%
nanimal=sum(unique(vec_animal)>0);
list_area=unique(vec_area(find(vec_resp)));
list_animal=unique(vec_animal(find(vec_resp)));
narea=numel(list_area);

ncell_all=[];
for iarea=1:narea
    for ianimal=1:nanimal
        idx_sel=mintersect(find(vec_resp),find(vec_area==list_area(iarea)),find(vec_animal==list_animal(ianimal)));
        if ~isempty(idx_sel)
            ncell_all(ianimal,iarea)=numel(idx_sel);
        end
    end
end
data_all={};
for iarea=1:narea
    nmax=max(ncell_all(:,iarea));
    mat_now=[];
    for ianimal=1:nanimal
        idx_sel=mintersect(find(vec_resp),find(vec_area==list_area(iarea)),find(vec_animal==list_animal(ianimal)));
        if ~isempty(idx_sel) && numel(idx_sel)>10 % discard empty and low cell count cases
            vec_now=nan(1,nmax);
            vec_now(1:numel(idx_sel))=vec_para(idx_sel);
            mat_now=cat(1,mat_now,vec_now);
        end
    end
    data_all{iarea}=mat_now;
end
end