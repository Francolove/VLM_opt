function winglet = build_winglet(root_chord, tip_chord, height,...
                                 R_raccordo, n_chord, n_height, cant,...
                                 sweep, toe_out, NACAs, wing, twist,up)
% costruisce un tensore N_corda x N_aperture x 3 che contiene le coordinate
% dei nodi dello winglet

% lo costruisco come estensione di ala e poi lo ruoto

    if length(NACAs) == 1
            NACAs = repmat(NACAs,1,n_height);
    end
    
    chord_vect = linspace(root_chord,tip_chord,n_height);
    % quante sezioni stanno nel tratto circolare?
    last_y = R_raccordo*(pi/2-cant);
    last_n_height = round(n_height*last_y/height);% suddivido in percentuale
                                                  % le sezioni tra raccordo
                                                  % e winglet
    yy = [linspace(0,last_y,last_n_height)...
          linspace(last_y,height,n_height-last_n_height+1)];
    %elimino il punto doppio
    yy = [yy(1:last_n_height) yy(last_n_height+2:end)];
    for i = 1:n_height
        ml = naca_mean_line(NACAs(i), n_chord, chord_vect(i));
        winglet(:,i,1) = ml(1,:);
        winglet(:,i,3) = ml(2,:);
    end
    % aggiungo toe_out
    winglet(:,:,2) = repmat(yy,n_chord,1);
    toe_out_vect = linspace(twist,twist+toe_out,n_height);   
    
    winglet(:,:,1) = winglet(:,:,1).*repmat(cos(toe_out_vect),n_chord,1) +...
                  winglet(:,:,3).*repmat(sin(toe_out_vect),n_chord,1);
    winglet(:,:,3) = -winglet(:,:,1).*repmat(sin(toe_out_vect),n_chord,1) +...
                  winglet(:,:,3).*repmat(cos(toe_out_vect),n_chord,1);
    % allineo il punto medio del profilo 
    winglet(:,:,1) = winglet(:,:,1) + repmat((root_chord -chord_vect)/2,n_chord,1);

    % aggiungo lo sweep 
    winglet(:,:,1) = winglet(:,:,1) + yy*sin(sweep);
    
    % ruoto lo winglet creando la circonferenza
    theta = -up*[linspace(0,pi/2-cant,last_n_height)...
             (pi/2-cant)*ones(1,n_height-last_n_height)];
    for i = 1:n_height-1
        winglet(:,i+1:end,2) = winglet(:,i,2) +...
                              (winglet(:,i+1:end,2)-winglet(:,i,2))*cos(theta(i+1)-theta(i)) + ...
                              (winglet(:,i+1:end,3)-winglet(:,i,3))*sin(theta(i+1)-theta(i));  
        winglet(:,i+1:end,3) = winglet(:,i,3) -...
                              (winglet(:,i+1:end,2)-winglet(:,i,2))*sin(theta(i+1)-theta(i)) + ...
                              (winglet(:,i+1:end,3)-winglet(:,i,3))*cos(theta(i+1)-theta(i));  
    end
    
    
    
    winglet(:,:,1) = winglet(:,:,1)+wing(1,end,1);
    winglet(:,:,2) = winglet(:,:,2)+wing(1,end,2);
end