% % Paso 0
% extrae los datos para guardar en un esquema de datos por anclaje
% GD, Marzo 2023
% Modificado para extraer los archivos pero ya en las fechas
% correspondientes para que todos sean horarios y estén sincronizadas las
% observaciones
% close all; clc; clear all;
% global seccion
% agrega al path de matlab la carpeta donde se encuentran los programas del
% procesado:
addpath(genpath('C:\Users\nucle\Documents\GitHub\Canek_VelocidadPorAnclaje\extrn\'));

% Ruta base donde estan los datos
base_adata = 'I:\TesisDoc\datos\CNKs_raw\';

% Ruta donde se van a guardar los datos generados por el programa
destino = 'I:\TesisDoc\datos\CNK_extractions\P0\';

seccion = 'Yucatan';
% seccion = 'Florida';

switch seccion
    case 'Yucatan'
        ancla = {'YUC1', 'YUC2', 'YUC3', 'YUC4', 'YUC5', 'YUCI5', 'YUC6', 'YUCI6',...
            'YUC7', 'YUCI7', 'YUC8', 'YUCI8', 'YUC9', 'YUCI9', 'YUC10'};
        %
    case 'Florida'
        ancla = {'EFL1', 'EFL2', 'EFL3', 'EFL4', ...
            'EFL5', 'EFLI5', 'EFL6', 'EFL7'};
end

% fechas de los caneks
cnkd = {['05-Jul-2012';'10-Jun-2013'];
    ['20-Jun-2013';'01-Jul-2014'];
    ['05-Jul-2014';'01-Jul-2015'];
    ['08-Jul-2015';'01-Aug-2016'];
    ['14-Aug-2016';'16-Jul-2018'];
    ['05-Aug-2018';'10-Nov-2020'];
    ['10-Aug-2022';'10-Jul-2023']};

CNK = [29  34   37  39  42  48  57]; % nombres de las campanias CNK
%===============================================================================
%%
for an = 4:4%length(ancla) % corre para cada anclaje
    
    for kk = 1 : length(CNK) % corre para cada canek
        
        anchor = struct(); q = 1; % estructura con datos de cada anclaje por canek
        cnk = num2str(CNK(kk)); %camapania canek
        
        % A).- Base del arbol. Aqui implemento una separacion en la
        % estructura del arbol ,quiza no es necesaria para fines practicos
        ruta_b = [base_adata, cnk,'\RECUPERACIONES\', seccion, '\'];
        fdi = dir(ruta_b);
        nam = {fdi.name}';
        
        %         encuentra el indice de la carpeta del anclaje, anade un gion para
        %         evitar duplicados como YUC1 y YUC10
        anid = find(~cellfun('isempty',strfind(nam,  [char(ancla(an)), '-'])));
        % =====================================================================
        if ~isempty(anid) % no todos los anclajes estan en todos los caneks
            
            nam = char(nam(anid)); % nombre de la carpeta del an-esimo anclaje
            % Primer nivel de la carpeta del anclaje
            ruta_1 = [ruta_b, nam, '\'];
            fdi = dir(ruta_1);
            [~, Lista]= eliminaDots( fdi ); % elimina directorios de punto y regresa solo la pura lista de nombres de directorios sin puntos
            % en este punto tenemos los
            
            % Aqui pasaremos al segundo nivel del directorio de anclaje
            % como es extraccion solo de velocidades, descartamos carpetas
            % de microcats o termistores si los hubiera.
            % Descarta otros sensores que no sean de velocidad (proximas 3 lineas):
            %                         instid = [find(~cellfun('isempty',strfind(Lista, 'ADCPS'))), ...
            %                         find(~cellfun('isempty',strfind(Lista, 'CORRPTS')))]; % encuentra indices de las carpetas de correntometria
            %                         Lista = Lista(instid); % descarta lo que no seacorrentometria
            
            % nos movemos en esta lista para ir bajando los niveles hasta llegar a los datos .mat
            for jj = 1 : length(Lista)
                
                inst = [ruta_1, char(Lista(jj)), '\'];
                fdi = dir([inst, '*.mat']); % en este caso no usamos los datos horarios
                inL = {fdi.name}'; % Lista de instrumentos
                
                if ~isempty(inL) % Por si en realidad no hubiese instrumentos o sus .mat
                    for qq = 1 :length(inL)
                        % aqui vamos a Zins
                        arMF = matfile([inst, char(inL(qq))]);
                        tiempo = arMF.jd;
                        
                        [tiempo, noRep] = Revisa_tiempo(tiempo); % a veces el tiempo tiene repetidos O.o
                        fecs = datenum(char(cnkd(kk, :)));
                        
                        tt = (fecs(1):1/24:fecs(end)); tt= tt(:);
                        Lon = arMF.Lon; Lat = arMF.Lat; % coordenadas del sensor
                        coors = nanmean([Lon, Lat]);
                        if isnan(sum(coors))
                            str1 = ['Ninguno de los instrumentos tiene coordenadas vaya a la \n'];
                            str2 = ['jodida bitacora y levante las coordenadas del anclaje: \n', char(ancla(an)), ' de CNK' num2str(cnk), '\n'];
                            fprintf(str1);
                            fprintf(str2);
                            Lat = input('latitud: ', 's');
                            Lon = input('longitud: ', 's');
                            coors = [str2num(lono), str2num(lato)];
                        end
                        
                        try
                            profe = abs(arMF.ProfEstimada);
                        catch
                            profe = abs(arMF.ProfDiseno);
                        end
                        
                        profd = abs(arMF.ProfDiseno);
                        Zins = abs(profe);
                        
                        banderasIVel % extrae variablers y hace todo el jale
                        
                        anchor(q).cnk = ['CNK-', cnk]; % campagna canek
                        anchor(q).sec = seccion; % seccion (Yucatan o Florida)
                        anchor(q).name = arMF.nam; % nombre del sensor
                        anchor(q).anc = char(ancla(an)); % nombre del anclaje

                        anchor(q).znom = abs(Zins); % z nominal del sensor
                        anchor(q).lola = [Lon, Lat]; % coordenadas del anclaje
                        
                        anchor(q).Pflg = P_f; % tiene sensor de preson? P_f = 1, otherwise P_f = 0
                        anchor(q).DPflg = 0; % De-drifted Pressure Flag
                        anchor(q).RPflg = 0; % Discarded Pressure Flag
                        anchor(q).pdis = abs(profd - profe); % diferencias entre prof de diseno y estimada
                        anchor(q).time = tt; % tiempo de mediciones
                        % banderas para rechazar Vel o Temp despues
                        anchor(q).RegVel = 0;
                        anchor(q).RegTemp = 0;
                        
                        anchor(q).p = Pr; % presion
                        anchor(q).z = Zr; % z
                        %  coordenadas para los bins
                        anchor(q).x = Xr; % lon
                        anchor(q).y = Yr; % lat
                        
                        anchor(q).u = Ur; % u
                        anchor(q).v = Vr; % v
                        anchor(q).w = Wr; % w
                        
                        anchor(q).T = Tr; % temperatura
                        anchor(q).S = Sr; % sal
                        anchor(q).O = Or; % oxigeno
                        
                        anchor(q).varis = misva; % lista de variables que en efecto existen (lo demas son nan)
                        anchor(q).duration = [min(tt), max(tt)]; % duracion de las mediciones
                        
                        q = q + 1 ;
                    end
                end
            end
            
            if isfield(anchor, 'lola')
                coors = [];
                for jjo = 1 : length(anchor)
                    coors(jjo, :) = [anchor(jjo).lola];
                end
                coors = nanmean(coors, 1);
                if isnan(sum(coors))
                    str1 = ['Ninguno de los instrumentos tiene coordenadas vaya a la \n'];
                    str2 = ['jodida bitacora y levante las coordenadas del anclaje: \n', char(ancla(an)), ' de CNK' num2str(cnk), '\n'];
                    fprintf(str1);
                    fprintf(str2);
                    lato = input('latitud: ', 's');
                    lono = input('longitud: ', 's');
                    coors = [str2num(lono), str2num(lato)];
                end
                for ww = 1 : length(anchor)
                    anchor(ww).lola = coors;
                end
            end
            
            tosa  = [destino, 'CNK',cnk, '_', char(ancla(an))];
            save(tosa, 'anchor', '-v7.3'); % guarda un archivo por anclaje por canek
            [char(ancla(an)), '. CNK ', cnk]
            
        end
    end
end
%anchor(1)