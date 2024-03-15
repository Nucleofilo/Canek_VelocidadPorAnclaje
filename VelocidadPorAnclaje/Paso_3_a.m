
% programas de procesado para base de datos por anclaje
% Paso 3a
% extrae los datos para guardar en un esquema de datos por anclaje
% GD, Marzo 2023. Modificado y actualizado en Enero del 2024
% Pone series de presion a sensores que no tienen interpolando el perfil de migracion vertical a cada tiempo

close all; clc; clear all;
% global seccion
% agrega al path de matlab la carpeta donde se encuentran los programas del
% procesado:
addpath(genpath('C:\Users\nucle\Documents\GitHub\Canek_VelocidadPorAnclaje\extrn\'));

cnkd = {['05-Jul-2012';'10-Jun-2013'];
    ['20-Jun-2013';'01-Jul-2014'];
    ['05-Jul-2014';'01-Jul-2015'];
    ['08-Jul-2015';'01-Aug-2016'];
    ['14-Aug-2016';'16-Jul-2018'];
    ['05-Aug-2018';'10-Nov-2020'];
    ['10-Aug-2022';'10-Jul-2023']};

seccion = 'Yucatan';
% seccion = 'Florida';
switch seccion
    case 'Yucatan'
        ancla = {'YUC1', 'YUC2', 'YUC3', 'YUC4', 'YUC5', ...
            'YUCI5', 'YUC6', 'YUCI6', 'YUC7', 'YUCI7', 'YUC8', ...
            'YUCI8', 'YUC9', 'YUCI9', 'YUC10'};
        % en yucatan comenzar en el 5
        ani = 1:length(ancla);
        
    case 'Florida'
        ancla = {'EFL1', 'EFL2', 'EFL3', 'EFL4', 'EFL5', 'EFLI5', 'EFL6', 'EFL7'};
        % en florida desde 2 hasta end-1
%         ani = 1:length(ancla)-1;
end

clc
CNK = [29  34, 37  39	42  48  57];
cols = gray(9); 
sflg = lower(seccion);
sflg = sflg(1:3);

raiz = 'I:\TesisDoc\datos\CNK_extractions\P0\';
destino = 'I:\TesisDoc\datos\CNK_extractions\P3\A\';
%%
warning off
for an = 13 : 13%length(ancla)
    
    figure('pos', [27 323 1579 581], 'color', 'w');
    cla;

    for kk = 4 : 4%length(CNK)
        
        cnk = num2str(CNK(kk));
        fecs = datenum(char(cnkd(kk, :)));
        
        hold on;
        ptc = patch( [ fecs(1), fecs(1), fecs(2), fecs(2), fecs(1) ], [200, -2200, -2200, 200, 200], cols(kk+1, :) );
        ptc.FaceAlpha = 0.2; ptc.EdgeColor = 'none';
        text( mean(fecs), 60, ['CNK-', cnk], 'HorizontalAlignment', 'center');
        drawnow
        colorbar

        try
            clear anchor
            load([raiz, 'CNK', cnk, '_', char(ancla(an))])
        catch
        end
        
        if exist('anchor') && ~isempty(fieldnames(anchor))
             name = anchor(1).name;
             tira = find(int16(name) == 45);
             tiran = str2num(name(tira(1)+2:tira(2)-1));
            
            zpo = [anchor(:).znom];
            flg = [anchor(:).Pflg];
            [zpo, ai] = sort(zpo, 'ascend');
            anchor = anchor(ai);
            flg = flg(ai);
            td =[]; z = []; v = [];
            t = anchor(1).time;

            flgp = logical([anchor(:).Pflg]);
            
            if any(~flgp) & tiran > 500
                anchor_si = anchor(flgp); % instrumentos que SI tienen sensor de presion
                zposi = [anchor_si(:).znom]';
                [zposi, ai] = sort(zposi, 'ascend');
                anchor_si = anchor_si(ai);
                pdis_no = nan(length(find(~flgp)), length(t));
                
                for tt = 1:length(t)
                    tiempo = fix(anchor(1).time(tt));
                    pnom = nan(size(anchor_si')); % profundìdad nominal de los isntrumentos
                    pinst = nan(size(anchor_si')); % desplazamiento de los instrumentos (prof. instantanea)
                    velis = nan(size(anchor_si'));
                    
                    for jj = 1 : length(anchor_si)
                        pnom(jj) = min(anchor_si(jj).p);
                        pinst(jj) =  anchor_si(jj).p(tt);
                        velis(jj) = sqrt(anchor_si(jj).v(tt, end)^2 + anchor_si(jj).u(tt, end)^2);
                    end
                    % si los instrumentos se encuentran a la misma profundidad se promedian
                    % sus desplazamientos:
                    %======================================================================
                    idx = find( abs( pnom - [pnom(2:end); 0] )  <= 5);
                    
                    for hh = 1 : length(idx)
                        pnom(idx(hh)+1) = mean(pnom([idx(hh), idx(hh) + 1]));
                        pnom(idx(hh)) = [];
                        pinst(idx(hh)+1) = mean(pinst([idx(hh), idx(hh) + 1]));
                        pinst(idx(hh)) = [];
                        velis(idx(hh)+1) = mean(velis([idx(hh), idx(hh) + 1]));
                        velis(idx(hh)) = [];
                    end
                    %======================================================================
                    % desplazamiento instantaneo de los instrumentos en el tiempo tt
                    pdis = (pnom(:)-pinst(:));
                    % instrumentos que no tienen sensor de presion:
                    anchor_no = anchor(~flgp);  % instrumentos que NO tienen sensor de presion
                    zpo = [anchor_no(:).znom]';
                    try
                        pdis_no(:, tt) = interp1([pnom; tiran+50], [pdis; 0], zpo, 'pchip', 'extrap');
                    catch
                        pdis_no(:, tt) = nan(size(zpo));
                    end
                    [tt, length(t)]
                end
                p_rec =  zpo + abs(pdis_no);
            
            %============ALGORITMO PARA PASAR PRESION A SENSORES QUE NO TIENEN============
            
            nos = find(~flgp);
            
            for kno = 1:length(nos)
                anchor(nos(kno)).p = p_rec(kno, :)';
                v = anchor(nos(kno)).v;
                
                nbin = size(v, 2);
                
                if nbin == 1
                    z = p_rec(kno, :)';
                    
                else
                    z = anchor(nos(kno)).z;
                    
                    for ll = 1 : nbin
                        z(:, ll) =  abs(z(:, ll)) +  abs(pdis_no(kno, :))';
                    end
                end
                anchor(nos(kno)).z = -abs(z);
            end
            end
            %==========================================================================
            
            for lll = 1 : length(anchor)
                v = anchor(lll).v;
                z = anchor(lll).z;
                
                z(abs(z)<70) = NaN;
                td = repmat( t(:), [1, size(v, 2)] );
                hold on;
                scatter(td(1:24:end), -abs(z(1:24:end)), [], v(1:24:end), 'filled');
                drawnow
            end
        tosav = [destino, 'CNK', cnk, '_', char(ancla(an)), '.mat'];
        save(tosav, 'anchor', '-v7.3')
        clear anchor
            
        end

    end
% end
% 
title([char(ancla(an)), '. Magnitud de la velocidad horizontal'])
xlim(datenum(['05-Jul-2012';'10-Jul-2023']))
datetick('x', 'mmm-yyyy');
xlim(datenum(['05-Jul-2012';'20-Jul-2023']))
ylim([-(tiran+100), 200])

box on;
grid on
cb = colorbar;
ylabel(cb, 'Rapidez (m/s)', 'Interpreter', 'latex')
% 
ylabel('Profundidad (m)')
xlabel('Tiempo (dias)')
caxis([-0.2, 1.3])
colormap jet;
% %     cmocean('balance', 'pivot', 0)
set(findall(gcf,'-property','Interpreter'),'Interpreter', 'latex');
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex');
Fs = findall(gcf,'-property','FontSize');
% 
for k = 1 : length(Fs)
    if strcmp(get(Fs(k), 'type'), 'axes')
        if Fs(k).FontSize < 12
            Fs(k).FontSize = 12;
        end
    end
end


drawnow
pause(0.1)
export_fig(['I:\TesisDoc\datos\CNK_extractions\P3\A\figs\', char(ancla(an))], '-m3')
pause(0.5)
% close all;
end
%% plot((zpo*0)+736556-20, -zpo, 'ok')
