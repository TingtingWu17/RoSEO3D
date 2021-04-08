function [lc_center,indx_out] = define_pathes(lc,thred)
indx_out = zeros(size(lc,2),1);

if size(lc,2)>1
for ii = 1:size(lc,2)
    lc_cur = lc(:,ii);

    distance = sqrt(sum((repmat(lc_cur,1,size(lc,2))-lc).^2,1));
    indx = distance<thred;
    lc_similar_box(:,ii) = mean(lc(:,indx),2);
    std_similar_box(:,ii) = std(lc(:,indx),1,2);
    indx_save(:,ii) = indx.';
        
end

ii = 0;
while sum(indx_save(:))~=0
    ii = ii+1;
    [~,cdt_idx] = max(sum(indx_save,1));
    lc_center(:,ii) = lc_similar_box(:,cdt_idx);
    indx_out(indx_save(cdt_idx,:))=ii;
    indx_save(:,indx_save(cdt_idx,:))=0;
end
lc_center = round(lc_center);

else
  lc_center = lc;
  indx_out = 1;
end