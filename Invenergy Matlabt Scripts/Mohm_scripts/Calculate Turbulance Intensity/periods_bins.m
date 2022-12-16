function [indices,lim] = periods_bins(data,nMonths,nHours)
%jbej 18032020
%Divide angles in sectors ( nsector+1 = all )

%data = 1 + (360-1).*rand(10000,1); 
%north=300;
%nSector=1

for m = 1:nMonths+1

    if m == nMonths+1
        indices_m{nMonths+1}=transpose([1:1:length(data)]);
        limM{m}="yearly";
    else
        indices_m{m}=find(month(data)==m);
        limM{m}=strcat('month:','{}',num2str(m));
    end
end

for h = 1:nHours+1

    if h == nHours+1
        indices_h{nHours+1}=transpose([1:1:length(data)]);
        limH{h}="daily";
    else
        indices_h{h} = find( hour(data)==(h-1));
        limH{h}=strcat('hour:','{}',num2str(h-1));
    end
end

for m = 1:nMonths+1
    for h = 1:nHours+1
        lim{m,h}=strcat(limM{m},'{}',limH{h});
        indices{m,h}=intersect(indices_m{m},indices_h{h});
    end
end

%disp('ending period')

