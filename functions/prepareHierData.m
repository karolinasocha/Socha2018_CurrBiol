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