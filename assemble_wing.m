function  assembled_wing = assemble_wing(wing1,wing2)
% assembla due tronchi di ala se il numero di punti in corda dei due
% tronchi è uguale (tensori con lo stesso numero di righe);
    if size(wing1,1) ~= size(wing2,1)
        disp('ERRORE! le matrici devono avere lo stesso numero di rige')
        return
    elseif wing1(:,end,:)~=wing2(:,1,:)
            disp('ERRORE! il primo elemento della seconda ala deve essere')
            disp('uguale all''ultimo della prima')
            return
    end
    assembled_wing = [wing1 wing2(:,2:end,:)];
    %prova
    
end