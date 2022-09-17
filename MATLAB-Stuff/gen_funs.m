classdef gen_funs
    methods (Static)
        function T=Tk(Tm,um)
            if nargin<2
                um='C';
            else
                %do nothing
            end
            s1='C';
            s2='F';
            sc=strcmp(s1,um);
            sf=strcmp(s2,um);
            if sc==1
                T=Tm+273.15;
            elseif sf==1
                T=(Tm+459.67)*5/9;
            else
                T=Tm+273.15;;
            end
        end
        function Pv=TensVapW(T,um)
            if nargin<2
                um='Pa';
            else
                %do nothing
            end
            s1='Pa';
            s2='mmHg';
            s3='atm';
            s4='bar';
            spa=strcmp(um,s1);
            shg=strcmp(um,s2);
            satm=strcmp(um,s3);
            sbar=strcmp(um,s4);
            Pv=exp(18.3036-3816.44/(T-46.13));
            if spa==1
                Pv=Pv/760*101325;
            elseif shg==1
                Pv=Pv;
            elseif satm==1
                Pv=Pv/760;
            elseif sbar==1
                Pv=Pv/760*1.01325;
            else
                %do nothing
            end
        end
        function s=Sc(mu,rho,diff)
            s=mu/(rho*diff);
        end
        function R=Rvalue(um)
            s1='SI';
            s2='L*atm';
            ssi=strcmp(um,s1);
            slatm=strcmp(um,s2);
            if ssi==1
                R=8.31446261815324;
            elseif slatm==1
                R=0.082057366080960;
            else
                R=8.31446261815324;
            end
        end
        function rh=rho(P,MM,T,perc)
            %INSTRUCTIONS

            %MM is in g/mol
            %P is in Pascals
            %T is in Kelvin
            %perc is expressed as 0. not as %
            if nargin<4
                perc=1;
            else
                %do nothing
            end
            MM=MM/1000;
            rh=perc*P*MM/(8.31446261815324*T);
        end
        function RePrSc=ReyPrSc(v,L,dv,rho,cp,k,diffusivity,reminder)
            if nargin<8
                reminder='no';
            else
                %do nothing
            end
            Re=v*L*rho/dv;
            Pr=dv*cp/k;
            Sc=dv/(diffusivity*rho);
            RePrSc=[Re;Pr;Sc];
            s1='yes';
            s2='no';
            syes=strcmp(reminder,s1);
            sno=strcmp(reminder,s2);
            if syes==1
                disp('first Reynolds, Second Prandtl,third Schmidt');
            elseif sno==1
                %do nothing
            else
                %do nothing
            end
        end
        function R=Re(rho,v,L,mu)
            R=v*rho*L/mu;
        end
        function q=RadiantFlux(T,eps)
            if nargin<2
                eps=1;
            else
                %do nothing
            end
            sig=5.670367e-8;
            q=eps*sig*T^4;
        end
        function P=Pa(Pm,um)
            if nargin<2
                um='atm';
            else
                %do nothing
            end
            s1='atm';
            s2='bar';
            s3='torr';
            s4='psi';
            s5='mmHg';
            satm=strcmp(s1,um);
            sbar=strcmp(s2,um);
            storr=strcmp(s3,um);
            spsi=strcmp(s4,um);
            smm=strcmp(s5,um);
            if satm==1
                P=Pm*101325;
            elseif sbar==1
                P=Pm*10^5;
            elseif storr==1
                P=Pm*101325/760;
            elseif smm==1
                P=Pm*101325/760;
            elseif spsi==1
                P=Pm*6894.76;
            else
                P=Pm*101325;
            end
        end
        function n = sigdigits(x)
            if( ~isnumeric(x) || ~isfinite(x) || isa(x,'uint64') )
                error('Need any finite numeric type except uint64');
            end
            if( x == 0 )
                n = 0;
                return;
            end
            x = abs(x);
            y = num2str(x,'%25.20e'); % Print out enough digits for any double
            z = [' ' y]; % Pad beginning to allow rounding spillover
            n = find(z=='e') - 1; % Find the exponent start
            e = n;
            while( str2double(y) == str2double(z) ) % While our number is still equal to our rounded number
                zlast = z;
                c = z(e); % Least significant printed digit
                if( c == '.' )
                    e = e - 1;
                    c = z(e);
                end
                z(e) = '0'; % 0 the least significant printed digit
                e = e - 1;
                if( c >= '5' ) % Round up if necessary
                    c = z(e);
                    if( c == '.' )
                        e = e - 1;
                        c = z(e);
                    end
                    while( true ) % The actual rounding loop
                        if( c == ' ' )
                            z(e) = '1';
                            break;
                        elseif( c < '9' )
                            z(e) = z(e) + 1;
                            break;
                        else
                            z(e) = '0';
                            e = e - 1;
                            c = z(e);
                            if( c == '.' )
                                e = e - 1;
                                c = z(e);
                            end
                        end
                    end
                end
            end
            n = n - 1;
            z = zlast(1:n); % Get rid of exponent
            while( z(n) == '0' ) % Don't count trailing 0's
                n = n - 1;
            end
            n = n - 2; % Don't count initial blank and the decimal point.
        end
        function M=MM(s,um)
            if nargin<2
                um='SI';
            else
                % do nothing
            end
            s4='SI';
            s5='g/mol';
            s1='air';
            s2='CO2';
            s3='oxygen';
            s6='water';
            s7='He';
            s8='H2';
            sair=strcmp(s,s1);
            sco2=strcmp(s,s2);
            sox=strcmp(s,s3);
            ssi=strcmp(um,s4);
            sg=strcmp(um,s5);
            sw=strcmp(s,s6);
            she=strcmp(s,s7);
            sh=strcmp(s,s8);
            if ssi==1
                if sair==1
                    M=0.0289647;
                elseif sox==1
                    M=0.032000;
                elseif sco2==1
                    M=0.044;
                elseif sw==1
                    M=0.018;
                elseif she==1
                    M=0.004;
                elseif sh==1
                    M=0.002;
                end
            elseif sg==1
                if sair==1
                    M=28.9647;
                elseif sox==1
                    M=32.000;
                elseif sco2==1
                    M=44;
                elseif sw==1
                    M=18;
                elseif she==1
                    M=4;
                elseif sh==1
                    M=2;
                end
            else
                if sair==1
                    M=0.0289647;
                elseif sox==1
                    M=0.032000;
                elseif sco2==1
                    M=0.044;
                elseif sw==1
                    M=0.018;
                elseif she==1
                    M=0.004;
                elseif sh==1
                    M=0.002;
                end
            end
        end
        function out=logb(x,base)
            if nargin<2
                base=exp(1);
            else
                %do nothing
            end
            out=log(x)/log(base);
        end
        function fvpa=fanning(v,L,dv,rho,eps,shape)
            % if the shape of the duct section isn't circular then you have to
            %input the hydraulic diameter, which is defined as the ratio of 4 times the
            %area of the section and the perimeter of the section (4*A/P)
            %the viscosity is the kinematic viscosity, written as nu.
            %the velocity has to be written as meters per seconds.
            %L is either the diameter or the hydraulic diameter, it depends on the
            %string input for shape, if you write 4 inputs it's the diameter else if
            %you input 'duct-like' or 'circular' it's considered to be the diameter of
            %the duct. If you wanna use the hydraulic diameter you have to input as
            %shape either 'non-duct-like' or 'non'.
            %eps is the roughness of the ducts.

            if nargin<5
                fvpa=0;
                disp('more inputs')
                return
            elseif nargin<6
                shape='duct-like';
            end
            s1='duct-like';
            s2='non-duct-like';
            s3='non';
            s4='circular';
            sd1=strcmp(shape,s1);
            sd2=strcmp(shape,s2);
            sd3=strcmp(shape,s3);
            sd4=strcmp(shape,s4);
            if sd1==1
                Re=v*L*rho/dv;
            elseif sd4==1
                Re=v*L*rho/dv;
            elseif sd3==1
                Re=v*L*rho/dv;
            elseif sd2==1
                Re=v*L*rho/dv;
            else
                fvpa=0;
                return
            end
            f0=(-4*log10((eps/L)/3.7)).^(-2);
            syms f
            eqn=1/sqrt(f)+4*log10((eps/L)/3.7+1.255/(Re*sqrt(f)));
            fvpa=vpasolve(eqn==0,f,f0);
        end
        function deltaT=deltaTln(deltain,deltaout)
            ln=log(deltain/deltaout);
            deltaT=(deltain-deltaout)/ln;
        end
        function R=CondRes(L1,L2,K,shape,measure,L3)
            if nargin<4
                R=0;
                disp('need more input')
                return
            elseif nargin<4
                shape='sheet';
                L3=L2;
                measure='d';
            elseif nargin<5
                measure='d';
                L3=L2;
            elseif nargin<6
                L3=L2;
            else
                %do nothing
            end
            s1='hollow spheres';
            s2='sheet';
            s3='duct';
            s4='r';
            s5='d';
            sd=strcmp(s5,measure);
            sr=strcmp(s4,measure);
            shs=strcmp(shape,s1);
            ss=strcmp(shape,s2);
            sduct=strcmp(shape,s3);
            if sr==1
                %do nothing
            elseif sd==1 && ss==0
                L1=L1/2;
                L2=L2/2;
            end
            if shs==1
                R=(L2-L1)/(4*pi*K*L1*L2);
            elseif sduct==1
                R=log(L2/L1)/(2*pi*L3*K);
            elseif ss==1
                R=L1/(K*L3);
            end
        end
    end
end

