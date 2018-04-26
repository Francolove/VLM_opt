function ml = naca_mean_line(NACA, n_nodes, chord)
% genera la linea media di un NACA 4 cifre
%
   if floor(NACA/10000)==0
      M = floor(NACA/1000);
      P = floor((NACA - 1000*M)/100);
      XX = floor(NACA - 1000*M - 100*P);

      p = P/10;
      m = M/100;

      x = linspace(0,chord,n_nodes);
      % x=chord*0.5*(1-cos(pi*linspace(0,1,n_nodes)));
      x_lead = x(x<p*chord);
      x_trail = x(x>=p*chord);
      ml(1,:) = x;
      ml(2,:) = [m/p^2 * ( 2*p*x_lead/chord-(x_lead/chord).^2),...
               m/(1-p)^2 *(1-2*p + 2*p*x_trail/chord-(x_trail/chord).^2)]*chord;
   else % 5cifre
      D1 = floor(NACA/10000);
      D2 = floor((NACA-1e4*D1)/1000);
      D3 = floor((NACA-1e4*D1-1e3*D2)/100);
      D4 = floor((NACA-1e4*D1-1e3*D2-1e2*D3)/10);
      D5 = NACA-1e4*D1-1e3*D2-1e2*D3-1e1*D4;
      DIGITS = D2*10+D3;

      x = linspace(0,chord,n_nodes);
      if D3==0; % no reflex
         DIGIT = [10 20 30 40 50];
         R = [0.0580 0.1260 0.2025 0.2900 0.3910];
         K1 = [361.400 51.640 15.957 6.643 3.230];

         r = R(DIGIT==DIGITS);
         k1 = K1(DIGIT==DIGITS);
         x_lead = x(x<r*chord);
         x_trail = x(x>=r*chord);
         ml(1,:) = x;
         ml(2,:) = [k1/6*((x_lead/chord).^3-3*r*(x_lead/chord).^2+r^2*(3-r)*(x_lead/chord)),...
                    k1*r^3/6*(1-(x_trail/chord))]*chord;

      else % reflex

         DIGIT = [21 31 41 51];
         R = [0.1300 0.2170 0.3180 0.4410];
         K1 = [51.990 15.793 6.520 3.191];
         K2_K1 = [0.000764 0.00677 0.0303 0.1355];

         r = R(DIGIT==DIGITS);
         k1 = K1(DIGIT==DIGITS);
         k2_k1 = K2_K1(DIGIT==DIGITS);

         x_lead = x(x<r*chord);
         x_trail = x(x>=r*chord);
         ml(1,:) = x;
         ml(2,:) = [k1/6*(((x_lead/chord)-r).^3-k2_k1*(1-r)^3*(x_lead/chord)-r^3*(x_lead/chord)+r^3),...
                    k1/6*(3*k2_k1*((x_trail/chord)-r).^3-k2_k1*(1-r)^3*(x_trail/chord)-r^3*(x_trail/chord)+r^3)];



      end
   end
end
