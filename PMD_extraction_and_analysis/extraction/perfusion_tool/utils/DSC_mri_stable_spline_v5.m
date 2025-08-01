function [res_ss]=DSC_mri_stable_spline_v5(conc,aif,mask,options)
% ultima modifica: Marco Castellaro 28/10/2010

%Funzione del pacchetto DSC_mri - DSC_mri_stable_spline_v5
%Autore: Castellaro Marco Gianluigi Pilonetto Denis Peruzzo - Universit� di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Flow (CBF) per un soggetto
%e per la deconvoluzione utilizza il metodo STABLE SPLINE
%
%Parametri in ingresso - conc (Matrice 4D), contiene gli andamenti delle
%                        concentrazioni DSC di tutti i voxel.
%                      - aif, andamento delle concentrazioni nel sito
%                        scelto come arteriale
%                      - mask (Matrice 3D), contiene la matrice utilizzata
%                        per mascherare il volume cerebrale da analizzare
%Options � la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%

% FASE 1: Preparazione dei dati
if size(aif,1)==1
    aif=aif';
end

nT=length(aif);
G = toeplitz(aif,[aif(1) zeros(1,nT-1)]); % Matrice della risposta impulsiva.

res_ss.map = zeros(options.nR,options.nC,options.nS);
res_ss.parametri = zeros(options.nR,options.nC,options.nS,3);
res_ss.ritardo = zeros(options.nR,options.nC,options.nS);
if options.deconv.SS.residual
    res_SS.residual=zeros(size(conc));
end

% Parametri algoritmo
bias_space          = 0; % 1: utilizza il bias space, 0: non lo utilizza
ritardo_max         = 0; % non stima il ritardo. 10; % Ritardo massimo accettabile dall'algoritmo
dati_x_stima_rumore = 2/3; %Frazione dei dati da usare per la stima del rumore


if options.waitbar
    hw_ss=waitbar(0,'Computing CBF by SS');
end



for r=1:options.nR
    if options.waitbar
        waitbar((r-1)/(options.nR),hw_ss);
    end
    for c=1:options.nC
        for s=1:options.nS
            if mask(r,c,s)
                
                % Specifiche dei parametri in ingresso/uscita
                % datiIn: variabile strutturata contenente
                %  - campioni: i dati da fittare
                %  - aif:      aif
                %
                % datiOut: variabile strutturata contente
                %  - residuo:       il residuo stimato
                %  - ritardo:       il ritardo stimato
                %  - varianza_dati: la varianza stimata sui dati
                
                campioni=reshape(conc(r,c,s,:),options.nT,1);
                
                if size(campioni,1)==1
                    campioni=campioni';
                end
                
                % FASE 2: stima del rumore nei dati
                % Gridotta=G(:,round(dati_x_stima_rumore*nT)); % Matrice della risposta impulsiva ridotta
                % [x,s,var_dati] = lscov(Gridotta,campioni); % Raw deconvolution: var_dati contiene una stima della varianza nei dati
                var_stimata=var_dati(campioni);
                sigma_dati=eye(nT,nT)*var_stimata;
                
                
                % FASE 3: stima del ritardo
                griglia_rit = (0:1:ritardo_max)';
                tab_residui   = zeros(nT,ritardo_max+1);
                tab_obiettivo = zeros(1,ritardo_max+1);
                tab_iper      = zeros(3,ritardo_max+1);
                for cont=1:(ritardo_max+1)
                    rit=griglia_rit(cont);
                    [res_stimato,obj_stimato,iper]=DNP(campioni',G,sigma_dati,rit,bias_space);
                    tab_residui(:,cont) = res_stimato;
                    tab_obiettivo(cont) = obj_stimato;
                    tab_iper(:,cont)    = iper;
                end
                [~,minPos]=min(tab_obiettivo);
                
                ritardo_stimato=griglia_rit(minPos);
                
                % FASE 4: stima della funzione residuo
                residuo = tab_residui(:,minPos);
                iper_par= tab_iper(:,minPos);
                
                % FASE 5: preparazione dell'output
                res_ss.parametri(r,c,s,:) = iper_par;
                res_ss.map(r,c,s)=max(residuo);
                res_ss.ritardo(r,c,s)=ritardo_stimato;
                if options.deconv.SS.residual
                    res_SS.residual(r,c,s,:)=residuo;
                end
                
            end
        end
    end
end

if options.waitbar
    delete(hw_ss)
end
end

function [var_stimata]=var_dati(campioni)
% Stima il livello di rumore (varianza) dei dati di DSC forniti in
% ingresso.
% Si basa sull'ipotesi che prima dell'arrivo del picco la concentrazione di
% tracciante dovrebbe essere nulla e usa quei campioni per stimare la
% varianza dei dati.
% L'arrivo del tracciante � stimato come l'ultimo campione prima del picco
% che sia inferiore al 15% dell'altezza del picco.

soglia=0.15; % l'arrivo del bolo � stimato quando la concentrazione di
% tracciante � inferiore al 15% dell'altezza del picco

[picco,posPicco]=max(campioni);
posArrivo=find(campioni(1:posPicco)<=soglia.*picco,1,'last');

var_stimata=var(campioni(1:posArrivo));
end

function [eta,obj,iper]=DNP(y,G,C_v,alfa,sel)

N=size(G,1);
Nl=size(G,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Hyperparameter vector estimation%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ipstart=[20  -2.4   -3.2];
[iper,obj]=fminsearch('minloglikDNP',ipstart,foptions,y,G,C_v,alfa,sel);
l2=exp(iper(1));beta=exp(iper(2));sigma2=1;%exp(iper(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Function estimation%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bias=exp([zeros(1,alfa) -beta*[1:(Nl-alfa)]])';
C=G*bias;
[C_z]=creaC_z(Nl,alfa,beta);
E=l2*G*C_z*G'+sigma2*C_v;
Ei=inv(E);
der=(sel==1)*inv(C'*Ei*C)*C'*Ei*y';
eta=der*bias+l2*C_z*G'*Ei*(y'-G*der*bias);
end

