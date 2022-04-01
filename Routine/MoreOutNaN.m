function UX = MoreOutNaN(Ux,X1)
GB_L = [0.95 0.75 0.5 0.25 0.1 0];         iO=1;           rEduced = 0.26;
while rEduced > 0.25 && iO <= length(GB_L)   
if length(X1)~=1
    GB_Out = ceil(GB_L(iO)/(abs(X1(2)-X1(1))));  %a GB is at least 0.75um in length
else
    GB_Out = ceil(GB_L(iO)/X1);
end
Ux(Ux==0)=NaN;          UX = Ux;
if GB_L(iO)~=0
for iG = 1:size(Ux,2)
    f = ind2sub(size(Ux(:,iG)),find(isnan(Ux(:,iG)))); % find NaNs
    if ~isempty(f) && length(f)<size(Ux,1)
        if min(f) == 1 && max(f) == size(Ux,1)
            iV = 1; while f(end-iV)==f(end)-iV && iV<length(f)-1; iV=iV+1; end
            iV = f(end-iV+1);
            To_Remove(iG,1) = iV-GB_Out; % from up
            iV = 1; while f(iV+1)==f(iV)+1 && iV<length(f)-1;     iV=iV+1; end
            To_Remove(iG,2) = iV+GB_Out; % from down
            if To_Remove(iG,1)<1;           To_Remove(iG,1)=1;          end
            if To_Remove(iG,2)>size(Ux,2);  To_Remove(iG,2)=size(Ux,2); end
            F = [1:To_Remove(iG,2), To_Remove(iG,1):size(Ux,1)]';
        elseif max(f) == size(Ux,1)
            if length(f) == 1;  iV = f; 
            else
                iV = 1; while f(end-iV)==f(end)-iV && iV<length(f)-1; iV=iV+1; end
                iV = f(iV);
            end
            To_Remove(iG,1) = iV-GB_Out; % from up
            if To_Remove(iG,1)<1;           To_Remove(iG,1)=1;          end
            F = [To_Remove(iG,1):size(Ux,1)]';
        elseif min(f) == 1
            if length(f) == 1;  iV = f; 
            else
                iV = 1; while f(iV+1)==f(iV)+1 && iV<length(f)-1;     iV=iV+1; end
            end
            iV = f(iV);
            To_Remove(iG,2) = iV+GB_Out; % from down
            if To_Remove(iG,2)>size(Ux,1);  To_Remove(iG,2)=size(Ux,1); end
            F = [1:To_Remove(iG,2)]';
        end
        try
        UX(F,ones(size(F)).*iG)=NaN;
        end
%         contourf(UX)
    end
end

for iG = 1:size(Ux,1)
    f = ind2sub(size(Ux(iG,:)),find(isnan(Ux(iG,:)))); % find NaNs
    if ~isempty(f) && length(f)<size(Ux,2)
        if min(f) == 1 && max(f) == size(Ux,2)
            iV = 1; while f(end-iV)==f(end)-iV && iV<length(f)-1; iV=iV+1; end
            iV = f(end-iV+1);
            To_Remove(iG,1) = iV-GB_Out; % from up
            iV = 1; while f(iV+1)==f(iV)+1 && iV<length(f)-1;     iV=iV+1; end
            To_Remove(iG,2) = iV+GB_Out; % from down
            if To_Remove(iG,1)<1;           To_Remove(iG,1)=1;          end
            if To_Remove(iG,2)>size(Ux,2);  To_Remove(iG,2)=size(Ux,2); end
            F = [1:To_Remove(iG,2), To_Remove(iG,1):size(Ux,2)]';
        elseif max(f) == size(Ux,2)
            if length(f) == 1;  iV = f; 
            else
                iV = 1; while f(end-iV)==f(end)-iV && iV<length(f)-1; iV=iV+1; end
                iV = f(iV);
            end
            
            To_Remove(iG,1) = iV-GB_Out; % from up
            if To_Remove(iG,1)<1;           To_Remove(iG,1)=1;          end
            F = [To_Remove(iG,1):size(Ux,2)]';
        elseif min(f) == 1
            if length(f) == 1;  iV = f; 
            else
                iV = 1; while f(iV+1)==f(iV)+1 && iV<length(f)-1;     iV=iV+1; end
            end
            iV = f(iV);
            To_Remove(iG,2) = iV+GB_Out; % from down
            if To_Remove(iG,2)>size(Ux,2);  To_Remove(iG,2)=size(Ux,2); end
            F = [1:To_Remove(iG,2)]';
        end
        UX(ones(size(F)).*iG,F)=NaN;
    end
end
end
    C_Ux = Ux;      C_Ux(isnan(C_Ux))=0;    C_el = nnz(C_Ux);
    C_UX = UX;      C_UX(isnan(C_UX))=0;    C_El = nnz(C_UX);
    rEduced = (C_el-C_El)/C_el;             iO=iO+1;
end
end
     