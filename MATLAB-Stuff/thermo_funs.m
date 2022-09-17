% per suggerimenti scrivetemi pure a fcapelizza@gmail.com
% per usarlo scrivete
% (nome scelto da voi = thermo_funs; 
% (risultato ottenuto) = th.(nomefunzione in thermo_funs)(input)
% EG
% th = thermo_funs;
% h_340K = th.delh(methane,delHmeth,340,'NIST')



classdef thermo_funs
    methods (Static)
        function out=delg(species,delGo,delHof,T,fcorr,npoints)
            % per dubbi/suggerimenti/bug mandatemi una mail a fcapelizza@gmail.com
            % fcorr e npoints sono argomenti opzionali
            % delg: energia libera di gibbs di una specie in fase GAS IDEALE a T in j/kmol
            % species: vettore con le costanti per le correlazioni dei cp, sono
            % supportati tre range di lavoro diversi per ciascuna correlazione, la
            % prima correlazione può essere o come quella del perry (esponenziale) o come quella del
            % NIST (polinomiale), le costanti da C1 a C5 sono quelle relative alla
            % prima correlazione, i dati da inserire sono quelli dati dal testo o dal
            % NIST, non vanno manipolate, la funzione fa già tutto.
            % La seconda e la terza correlazione sono sempre del tipo NIST, sono
            % richiesti i coefficienti dalla A alla E.
            % gli ultimi due elementi del vettore species (da 17 elementi) sono le
            % temperature massime di lavoro della prima e della seconda correlazione,
            % se non ci sono 3 correlazioni sul NIST (mi sembra ci siano solo per i gas
            % biatomici) mettete al posto dei 5 coefficienti della terza correlazione
            % degli 0, in realtà potete mettere quello che volete, non cambia nulla.
            % delGo: energia libera di gibbs della specie a 298.15 K
            % delHof: l'entalpia della specie a 298.15 K in j/kmol
            % T: temperatura finale di integrazione in K
            % fcorr: STRINGA che può essere o 'PERRY' o 'NIST' che specifica il tipo di correlazione
            % usata nella prima integrazione, se non mettete niente usa quella del
            % perry
            % npoints: numero di elementi usati per l'integrazione numerica della
            % funzione delh, se non impostato usa 200 punti di integrazione

            % BASI TEORICHE: integro l'equazione di Van't Hoff su un range di T facendo
            % variare l'intervallo per l'integrazione fatta da delh.
            
            %VALIDAZIONE: tramite diagrammi di francis per metano, etano, propano,
            %etilene e acetilene.

            %ATTENZIONE: se provate a vedere come varia delg in funzione della
            %temperatura vedrete che per i composti con delHof<0 sale e poi riscende,
            %con un andamento parabolico, mentre per composti elementari (delHof=0)
            %questo scende dopo i 298.15 K, è normale, è dato dal fatto che
            %d(delg/(R*T))/dT =-delH/(R*T^2) , se delH(T) è positivo delg/(R*T)
            %diminuisce col variare della temperatura.

            if nargin<6
                npoints = 200;
            elseif nargin<5
                npoints = 200;
                fcorr = 'PERRY';
            elseif nargin<4
                error('numero di input non sufficienti')
            else
                %do nothing
            end
            Tref = 298.15;
            R = 8314.5;
            Xs = linspace(Tref,T,npoints);
            Ys = zeros(1,npoints);
            th = thermo_funs;
            for i=1:length(Ys)
                Ys(i) = th.delh(species,delHof,Xs(i),fcorr)./(R*Xs(i).^2);
            end
            val_int = trapz(Xs,Ys);
            out = (delGo/(R*Tref) - val_int)*R*T;
        end

        function output=delh(species,delHof,Tend,fcorr)
            % per dubbi/suggerimenti/bug mandatemi una mail a fcapelizza@gmail.com
            % fcorr è opzionale come input
            % delh: entalpia di GAS IDEALE in j/kmol
            % species: vettore con le costanti per le correlazioni dei cp, sono
            % supportati tre range di lavoro diversi per ciascuna correlazione, la
            % prima correlazione può essere o come quella del perry (esponenziale) o come quella del
            % NIST (polinomiale), le costanti da C1 a C5 sono quelle relative alla
            % prima correlazione, i dati da inserire sono quelli dati dal testo o dal
            % NIST, non vanno manipolate, la funzione fa già tutto.
            % La seconda e la terza correlazione sono sempre del tipo NIST, sono
            % richiesti i coefficienti dalla A alla E.
            % gli ultimi due elementi del vettore species (da 17 elementi) sono le
            % temperature massime di lavoro della prima e della seconda correlazione,
            % se non ci sono 3 correlazioni sul NIST (mi sembra ci siano solo per i gas
            % biatomici) mettete al posto dei 5 coefficienti della terza correlazione
            % degli 0, in realtà potete mettere quello che volete, non cambia nulla.
            % delHof: l'entalpia della specie a 298.15 K in j/kmol
            % Tend: temperatura finale di integrazione in K
            % fcorr: STRINGA che può essere o 'PERRY' o 'NIST' che specifica il tipo di correlazione
            % usata nella prima integrazione, se non mettete niente usa quella del
            % perry
            if nargin<4
                fcorr = 'PERRY';
            else
                % do nothing
            end
            C1 = species(1);
            C2 = species(2);
            C3 = species(3);
            C4 = species(4);
            C5 = species(5);
            A =species(6);
            B =species(7);
            C =species(8);
            D = species(9);
            E = species(10);
            Tmax1 = species(16);
            Tmax2 = species(17);
            A2 =species(11);
            B2 =species(12);
            C_2 =species(13);
            D2 = species(14);
            E2 = species(15);
            Tref = 298.15;
            sp = strcmp('PERRY',fcorr);
            sn = strcmp('NIST',fcorr);
            if sn==0 && sp==0
                error('problema in correlazione, riprova con "NIST" o "PERRY".')
            else
                %do nothing
            end
            if sp==1
                if Tend > Tmax1 && Tend<Tmax2
                    intcp1 = C1*Tmax1 - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tmax1) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tmax1) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
                    intcp2 = ((B*Tend.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tend.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tend.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)./Tend + (1000000*E)./Tmax1 + A*Tend - A*Tmax1)*10^3;
                    intcp = intcp1+intcp2;
                elseif Tend > Tmax1 && Tend>Tmax2
                    intcp1 = C1*Tmax1 - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tmax1) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tmax1) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
                    intcp2 = ((B*Tmax2.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tmax2.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tmax2.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)./Tmax2 + (1000000*E)./Tmax1 + A*Tmax2 - A*Tmax1)*10^3;
                    intcp3 =((B2*Tend.^2)/2000 - (B2*Tmax2.^2)/2000+ (C_2*Tend.^3)/3000000 - (C_2*Tmax2.^3)/3000000 + (D2*Tend.^4)/4000000000 - (D2*Tmax2.^4)/4000000000 - (1000000*E2)./Tend + (1000000*E2)/Tmax2 + A2*Tend - A2*Tmax2)*10^3;
                    intcp = intcp1+intcp2+intcp3;
                else
                    intcp1 = C1*Tend - C1*Tref - (2*C2*C3)/(exp(-(2*C3)/Tend) - 1) - (2*C4*C5)/(exp(-(2*C5)/Tend) + 1) + (2*C2*C3)/(exp(-(2*C3)/Tref) - 1) + (2*C4*C5)/(exp(-(2*C5)/Tref) + 1);
                    intcp = intcp1;
                end
            elseif sn==1
                if Tend > Tmax1 && Tend<Tmax2
                    intcp1 = ((C2*Tmax1.^2)/2000 + (C3*Tmax1.^3)/3000000 - (1000000*C5)/Tmax1 + (C4*Tmax1.^4)/4000000000 - (C2*Tref.^2)/2000 - (C3*Tref.^3)/3000000 + (1000000*C5)/Tref - (C4*Tref.^4)/4000000000 + C1*Tmax1 - C1*Tref)*10^3;
                    intcp2 = ((B*Tend.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tend.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tend.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)/Tend + (1000000*E)/Tmax1 + A*Tend - A*Tmax1)*10^3;
                    intcp = intcp1+intcp2;
                elseif Tend > Tmax1 && Tend>Tmax2
                    intcp1 = ((C2*Tmax1.^2)/2000 + (C3*Tmax1.^3)/3000000 - (1000000*C5)/Tmax1 + (C4*Tmax1.^4)/4000000000 - (C2*Tref.^2)/2000 - (C3*Tref.^3)/3000000 + (1000000*C5)/Tref - (C4*Tref^4)/4000000000 + C1*Tmax1 - C1*Tref)*10^3;
                    intcp2 = ((B*Tmax2.^2)/2000 - (B*Tmax1.^2)/2000+ (C*Tmax2.^3)/3000000 - (C*Tmax1.^3)/3000000 + (D*Tmax2.^4)/4000000000 - (D*Tmax1.^4)/4000000000 - (1000000*E)/Tmax2 + (1000000*E)/Tmax1 + A*Tmax2 - A*Tmax1)*10^3;
                    intcp3 =((B2*Tend.^2)/2000 - (B2*Tmax2.^2)/2000+ (C_2*Tend.^3)/3000000 - (C_2*Tmax2.^3)/3000000 + (D2*Tend.^4)/4000000000 - (D2*Tmax2.^4)/4000000000 - (1000000*E2)/Tend + (1000000*E2)/Tmax2 + A2*Tend - A2*Tmax2)*10^3;
                    intcp = intcp1+intcp2+intcp3;
                else
                    intcp1 = ((C2*Tend.^2)/2000 + (C3*Tend.^3)/3000000 - (1000000*C5)/Tend + (C4*Tend.^4)/4000000000 - (C2*Tref.^2)/2000 - (C3*Tref.^3)/3000000 + (1000000*C5)/Tref - (C4*Tref.^4)/4000000000 + C1*Tend - C1*Tref)*10^3;
                    intcp = intcp1;
                end
            end
            output = delHof + intcp;
        end
        function z = ZVdW(index,T,P,ph,Pc,Tc)
            % P and Pc must have the same unit
            R = 8.3144621;
            a = 0.421875*R^2*Tc(index)^2/Pc(index);
            b = 0.125*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1-B A -A*B];
            z = roots(polinomio); %cerco le radici
            % z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
        end

        function z = ZRK(index,T,P,ph,Pc,Tc)
            R = 8.3144621;
            k = (Tc(index)/T)^0.5;
            a = 0.42748*R^2*Tc(index)^2*k/Pc(index);
            b = 0.08664*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
        end
        function z = ZRKS(index,T,P,ph,Pc,Tc,w)
            R = 8.3144621;
            S = 0.48+1.574*w(index)-0.176*w(index)^2;
            k = (1+S*(1-(T./Tc(index)).^0.5)).^2;
            a = 0.42748*R^2*Tc(index)^2*k/Pc(index);
            b = 0.08664*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapor
            end
        end
        function z = ZPR(index,T,P,ph,Pc,Tc,w)
            R = 8.3144621;
            S = 0.37464+1.54226*w(index)-0.26992*w(index)^2;
            k = (1+S*(1-(T/Tc(index))^0.5))^2;
            a = 0.45724*R^2*Tc(index)^2*k/Pc(index);
            b = 0.0778*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
        end
        function z = ZVIR(index,T,P,Pc,Tc,w)
            Tr = T/Tc(index);
            pr = P/Pc(index);
            B0 = 0.083 - 0.422*Tr^(-1.6);
            B1 = 0.139 - 0.172*Tr^(-4.2);
            beta = B0 + B1*w(index);
            z = 1 + beta*pr/Tr;
        end

        % functions for hr
        
        function hr = hrVdW(index,T,P,ph,Pc,Tc)
            R = 8.3144621;
            a = 0.421875*R^2*Tc(index)^2/Pc(index);
            b = 0.125*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1-B A -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            hr = z-1-A/z;
            hr = hr*R*T;
        end
        function hr = hrRK(index,T,P,ph,Pc,Tc)
            R = 8.3144621;
            k = (Tc(index)/T)^0.5;
            a = 0.42748*R^2*Tc(index)^2*k/Pc(index);
            b = 0.08664*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            % z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            hr = z-1-3*A*log((z+B)/z)/(2*B);
            hr = hr*R*T;
        end

        function hr = hrRKS(index,T,P,ph,Pc,Tc,w)
            R = 8.3144621;
            Tr = T/Tc(index);
            S = 0.48+1.574*w(index)-0.176*w(index)^2;
            k = (1+S*(1-(T/Tc(index))^0.5))^2;
            a = 0.42748*R^2*Tc(index)^2*k/Pc(index);
            b = 0.08664*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            hr = z-1-A*(1+S*Tr^0.5/(1+S*(1-Tr^0.5)))*log((z+B)/z)/B;
            hr = hr*R*T;
        end
        function hr = hrPR(index,T,P,ph,Pc,Tc,w)
            R = 8.3144621;
            Tr = T/Tc(index);
            S = 0.37464+1.54226*w(index)-0.26992*w(index)^2;
            k = (1+S*(1-(T/Tc(index))^0.5))^2;
            a = 0.45724*R^2*Tc(index)^2*k/Pc(index);
            b = 0.0778*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            hr = z-1-A*(1+S*Tr^0.5/(1+S*(1-Tr^0.5)))*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B);
            hr = hr*R*T;
        end
        function hr = hrVIR(index,T,P,Pc,Tc,w)
            R = 8.3144621;
            Tr = T/Tc(index);
            pr = P/Pc(index);
            B0 = 0.083 - 0.422*Tr^(-1.6);
            B1 = 0.139 - 0.172*Tr^(-4.2);
            beta = B0 + B1*w(index);
            B = beta*R*Tc(index)/Pc(index);
            z = 1 + beta*pr/Tr;
            dBdT = R*(0.675/Tr^2.6 + w(index)*0.722/Tr^5.2)/Pc(index);
            hr = P*(B-T*dBdT);
        end
        % functions for phi
        function phi = phiVdW(index,T,P,ph,Pc,Tc)
            R = 8.3144621;
            a = 0.421875*R^2*Tc(index)^2/Pc(index);
            b = 0.125*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1-B A -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            phi = z-1-A/z-log(z-B);
            phi = exp(phi);
        end
        
        function phi = phiRK(index,T,P,ph,Pc,Tc)
            R = 8.3144621;
            k = (Tc(index)/T)^0.5;
            a = 0.42748*R^2*Tc(index)^2*k/Pc(index);
            b = 0.08664*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            phi = z-1-A*log((z+B)/z)/B-log(z-B);
            phi = exp(phi);
        end
        function phi = phiRKS(index,T,P,ph,Pc,Tc,w)
            R = 8.3144621;
            Tr = T/Tc(index);
            S = 0.48+1.574*w(index)-0.176*w(index)^2;
            k = (1+S*(1-(T/Tc(index))^0.5))^2;
            a = 0.42748*R^2*Tc(index)^2*k/Pc(index);
            b = 0.08664*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            phi = z-1-A*log((z+B)/z)/B-log(z-B);
            phi = exp(phi);
        end
        function phi = phiPR(index,T,P,ph,Pc,Tc,w)
            R = 8.3144621;
            Tr = T/Tc(index);
            S = 0.37464+1.54226*w(index)-0.26992*w(index)^2;
            k = (1+S*(1-(T/Tc(index))^0.5))^2;
            a = 0.45724*R^2*Tc(index)^2*k/Pc(index);
            b = 0.0778*R*Tc(index)/Pc(index);
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if strcmp(ph,'L')==1
                z = min(z); %z liquid
            else
                z = max(z); %z vapore
            end
            phi = z-1-A*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B)-log(z-B);
            phi = exp(phi);
        end

        function phi = phiVIR(index,T,P,Pc,Tc,w)
            R = 8.3144621;
            Tr = T/Tc(index);
            pr = P/Pc(index);
            B0 = 0.083 - 0.422*Tr^(-1.6);
            B1 = 0.139 - 0.172*Tr^(-4.2);
            beta = B0 + B1*w(index);
            B = beta*R*Tc(index)/Pc(index);
            phi = P*B/(R*T);
            phi = exp(phi);
        end
        % functions entalphy for non ideal mixes
        function hr = hrVdWmix(T,P,x,state,Pc,Tc)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                a(i) = 0.421875*R^2*Tc(i)^2/Pc(i);
                b(i) = 0.125*R*Tc(i)/Pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*P/(R*T)^2;
            B = bmix*P/(R*T);
            polinomio = [1 -1-B A -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            hr = z-1-A/z;
            hr = hr*R*T;
        end

        function hr = hrRKmix(T,p,x,state,pc,Tc)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                Tr(i) = T/Tc(i);
                k(i) = 1/Tr(i)^0.5;
                a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i);
                b(i) = 0.08664*R*Tc(i)/pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*p/(R*T)^2;
            B = bmix*p/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            hr = z-1-3*A*log((z+B)/z)/(2*B);
            hr = hr*R*T;
        end

        function hr = hrRKSmix(T,p,x,state,pc,Tc,w)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                S(i) = 0.48+1.574*w(i)-0.176*w(i)^2;
                Tr(i) = T/Tc(i);
                k(i) = (1+S(i)*(1-Tr(i)^0.5))^2;
                a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i);
                b(i) = 0.08664*R*Tc(i)/pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*p/(R*T)^2;
            B = bmix*p/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            coeff = 0;
            for i=1:length(x)
                coeff = coeff + x(i)*S(i)*(a(i)*Tr(i)/k(i))^0.5;
            end
            hr = z-1-A*(1+coeff/amix^0.5)*log((z+B)/z)/B;
            hr = hr*R*T;
        end

        function hr = hrPRmix(T,p,x,state,pc,Tc,w)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                S(i) = 0.37464+1.54226*w(i)-0.26992*w(i)^2;
                Tr(i) = T/Tc(i);
                k(i) = (1+S(i)*(1-Tr(i)^0.5))^2;
                a(i) = 0.45724*R^2*Tc(i)^2*k(i)/pc(i);
                b(i) = 0.0778*R*Tc(i)/pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*p/(R*T)^2;
            B = bmix*p/(R*T);
            polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            coeff = 0;
            for i=1:length(x)
                coeff = coeff + x(i)*S(i)*(a(i)*Tr(i)/k(i))^0.5;
            end
            hr = z-1-A*(1+coeff/amix^0.5)*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B);
            hr = hr*R*T;
        end

        function hr = hrVIRmix(Temp,press,x,pc,Tc,w,zc)
            R = 8.3144621;
            vc = R*Tc.*zc./pc;
            %Z
            for i=1:length(Tc)
                for j=1:length(Tc)
                    Z(i,j) = (zc(i)+zc(j))/2;
                end
            end
            %V
            for i=1:length(Tc)
                for j=1:length(Tc)
                    V(i,j) = ((vc(i)^(1/3)+vc(j)^(1/3))/2)^3;
                end
            end
            %k - OKKIO KE NON SEMPRE C'E'
            for i=1:length(Tc)
                for j=1:length(Tc)
                    k(i,j) = 1-(V(i,i)*V(j,j))^(1/2)/V(i,j);
                end
            end
            %T
            for i=1:length(Tc)
                for j=1:length(Tc)
                    T(i,j) = (Tc(i)*Tc(j))^(1/2)*(1-k(i,j));
                end
            end
            %P
            for i=1:length(Tc)
                for j=1:length(Tc)
                    P(i,j) = Z(i,j)*R*T(i,j)/V(i,j);
                end
            end
            %W
            for i=1:length(Tc)
                for j=1:length(Tc)
                    W(i,j) = (w(i)+w(j))/2;
                end
            end
            TR = Temp./T;
            pR = press./P;
            B0 = 0.083 - 0.422*TR.^(-1.6);
            B1 = 0.139 - 0.172*TR.^(-4.2);
            beta1 = B0 + B1.*W;
            beta1 = beta1.*pR./TR;
            beta = 0;
            for i=1:length(Tc)
                for j=1:length(Tc)
                    beta = beta + x(i)*x(j)*beta1(i,j);
                end
            end
            z = 1 + beta;
            B2 = 0.675./TR.^2.6;
            B3 = 0.722./TR.^5.2;
            dBdtpR = pR.*(B2 + W.*B3);
            dBdt = 0;
            for i=1:length(Tc)
                for j=1:length(Tc)
                    dBdt = dBdt + x(i)*x(j)*dBdtpR(i,j);
                end
            end
            hr = (beta - dBdt)*R*Temp;
        end

        % functions for phi in non ideal mixes

        function phi = phiVdWmix(T,p,x,index,state,pc,Tc,w)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                a(i) = 0.421875*R^2*Tc(i)^2/pc(i);
                b(i) = 0.125*R*Tc(i)/pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*p/(R*T)^2;
            B = bmix*p/(R*T);
            Ai = a(index)*p/(R*T)^2;
            Bi = b(index)*p/(R*T);
            polinomio = [1 -1-B A -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            phi = Bi/(z-B)-2*(Ai*A)^0.5/z-log(z-B);
            phi = exp(phi);
        end

        function phi = phiRKmix(T,p,x,index,state,pc,Tc,w)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                Tr(i) = T/Tc(i);
                k(i) = 1/Tr(i)^0.5;
                a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i);
                b(i) = 0.08664*R*Tc(i)/pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*p/(R*T)^2;
            B = bmix*p/(R*T);
            Ai = a(index)*p/(R*T)^2;
            Bi = b(index)*p/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            phi = Bi*(z-1)/B+A*(Bi/B-2*(Ai/A)^0.5)*log((z+B)/z)/B-log(z-B);
            phi = exp(phi);
        end

        function phi = phiRKSmix(T,p,x,index,state,pc,Tc,w)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                S(i) = 0.48+1.574*w(i)-0.176*w(i)^2;
                Tr(i) = T/Tc(i);
                k(i) = (1+S(i)*(1-Tr(i)^0.5))^2;
                a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i);
                b(i) = 0.08664*R*Tc(i)/pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*p/(R*T)^2;
            B = bmix*p/(R*T);
            Ai = a(index)*p/(R*T)^2;
            Bi = b(index)*p/(R*T);
            polinomio = [1 -1 A-B-B^2 -A*B];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            coeff = 0;
            for i=1:length(x)
                coeff = coeff + x(i)*S(i)*(a(i)*Tr(i)/k(i))^0.5;
            end
            phi = Bi*(z-1)/B+A*(Bi/B-2*(Ai/A)^0.5)*log((z+B)/z)/B-log(z-B);
            phi = exp(phi);
        end

        function phi = phiPRmix(T,p,x,index,state,pc,Tc,w)
            amix = 0;bmix = 0;
            R = 8.3144621;
            for i=1:length(x)
                S(i) = 0.37464+1.54226*w(i)-0.26992*w(i)^2;
                Tr(i) = T/Tc(i);
                k(i) = (1+S(i)*(1-Tr(i)^0.5))^2;
                a(i) = 0.45724*R^2*Tc(i)^2*k(i)/pc(i);
                b(i) = 0.0778*R*Tc(i)/pc(i);
            end
            for i=1:length(x)
                for j=1:length(x)
                    amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5;
                    bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2;
                end
            end
            A = amix*p/(R*T)^2;
            B = bmix*p/(R*T);
            Ai = a(index)*p/(R*T)^2;
            Bi = b(index)*p/(R*T);
            polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3];
            z = roots(polinomio); %cerco le radici
            z = z(imag(z)==0); %elimino le radici immaginarie
            if state=='L'
                z = min(z); %z liquido
            else
                z = max(z); %z vapore
            end
            phi = Bi*(z-1)/B+A*(Bi/B-2*(Ai/A)^0.5)*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B)-log(z-B);
            phi = exp(phi);
        end

        function phi = phiVIRmix(Temp,press,x,index,pc,Tc,w)
            R = 8.3144621;
            zc = zeros(1,length(Tc));
            for i=1:length(Tc)
                zc(i) = 1 - 0.339 - 0.033*w(i);
            end
            vc = R*Tc.*zc./pc;
            %Z
            for i=1:length(Tc)
                for j=1:length(Tc)
                    Z(i,j) = (zc(i)+zc(j))/2;
                end
            end
            %V
            for i=1:length(Tc)
                for j=1:length(Tc)
                    V(i,j) = ((vc(i)^(1/3)+vc(j)^(1/3))/2)^3;
                end
            end
            %k - OKKIO KE NON SEMPRE C'E'
            for i=1:length(Tc)
                for j=1:length(Tc)
                    k(i,j) = 1-(V(i,i)*V(j,j))^(1/2)/V(i,j);
                end
            end
            %T
            for i=1:length(Tc)
                for j=1:length(Tc)
                    T(i,j) = (Tc(i)*Tc(j))^(1/2)*(1-k(i,j));
                end
            end
            %P
            for i=1:length(Tc)
                for j=1:length(Tc)
                    P(i,j) = Z(i,j)*R*T(i,j)/V(i,j);
                end
            end
            %W
            for i=1:length(Tc)
                for j=1:length(Tc)
                    W(i,j) = (w(i)+w(j))/2;
                end
            end
            TR = Temp./T;
            pR = press./P;
            B0 = 0.083 - 0.422*TR.^(-1.6);
            B1 = 0.139 - 0.172*TR.^(-4.2);
            beta1 = B0 + B1.*W;
            beta1 = beta1.*pR./TR;
            beta = 0;
            for i=1:length(Tc)
                for j=1:length(Tc)
                    beta = beta + x(i)*x(j)*beta1(i,j);
                end
            end
            beta2 = 0;
            for i=1:length(Tc)
                beta2 = beta2 + x(i)*beta1(index,i);
            end
            phi = -beta + 2*beta2;
            phi = exp(phi);
        end
        %WILSON
        function gam=gamma(x,Aij)
            gam=exp(1-log(x*Aij')-sum((x'.*Aij)./(x*Aij')'));
        end
    end
end