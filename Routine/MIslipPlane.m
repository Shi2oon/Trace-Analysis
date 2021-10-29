function [MI]=MIslipPlane(CS,ebsd)
%% Identify the slip plane!!
    % These are systems of hkl being ing
    % 100, 110, 111, 200, 211
for i=1:length(CS)
    h = [...
    Miller(CS{i}.aAxisRec.hkl(1),CS{i}.aAxisRec.hkl(2),CS{i}.aAxisRec.hkl(3),...
        ebsd(ebsd.mineralList{ebsd.indexedPhasesId(i)}).CS), ...
    Miller(CS{i}.bAxisRec.hkl(1),CS{i}.bAxisRec.hkl(2),CS{i}.bAxisRec.hkl(3),...
        ebsd(ebsd.mineralList{ebsd.indexedPhasesId(i)}).CS), ....
    Miller(CS{i}.cAxisRec.hkl(1),CS{i}.cAxisRec.hkl(2),CS{i}.cAxisRec.hkl(3),...
        ebsd(ebsd.mineralList{ebsd.indexedPhasesId(i)}).CS)...
                                ];
     eval(sprintf('MI.h%d=h;',i));
end