function ISI=compute_ISI(units,startt,endt)
%%ISI is a N element vector (N=number of units) containing the mean ISI 
% of each unit
matrix=spikest2vector(units,startt,endt);
ISI=NaN(size(matrix,1),1);
ISItmp=[];
for i=1:size(matrix,1)
idx=find(matrix(i,:)>0);
% if numel(idx)>1
% ISI(i)=mean(diff(idx));
% end

ISItmp=[ISItmp;diff(idx)'];
end
ISI=ISItmp;

end