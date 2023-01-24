function matrix=spikest2vector(units,startt,endt)
%matrix rows contain units and columns represent t within trial
matrix=zeros(size(units,2),round((endt-startt)*1000));
for unit=1:size(units,2)
    idxsp=find(((units(unit).ts>=startt) & (units(unit).ts<endt)));
    tsp=round((units(unit).ts(idxsp)-startt)*1000)+1;
    matrix(unit,tsp)=1;
    
    if numel(unique(tsp))<numel(tsp)
        disp('double spike')
        N = histcounts(tsp,0:max(tsp));
        double_spike=find(N>1);
    end
end

end