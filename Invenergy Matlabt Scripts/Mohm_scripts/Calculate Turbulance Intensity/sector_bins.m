function [indices,lim] = sector_bins(data,nSector,north)
%jbej 18032020
%Divide angles in sectors ( nsector+1 = all )

%data = 1 + (360-1).*rand(10000,1); 
%north=300;
%nSector=12

data=mod(data,360);

if nSector==0
    dSector=-999
    shift=-999
else
    dSector=360/nSector;
    shift=dSector/2-north;
end

for i = 1:nSector+1
    if i == nSector+1
        indices{nSector+1}=transpose([1:1:length(data)]);
        lim{nSector+1}='omnidirectional';
    else
        indices{i} = find( (i-1)*dSector <= mod(data+shift,360) & mod(data+shift,360) < i*dSector );
        ul=mod(i*dSector-shift,360);
        if ul==0 
            ul=360;
        end
        lim{i}=strcat(num2str(mod((i-1)*dSector-shift,360)),'\circ -',num2str(ul),'\circ');
    end
%    histogram(mod(data(indices{i}),360),[0:5:360])
%    xticks([0:5:360])
%    title(i)
end


end
