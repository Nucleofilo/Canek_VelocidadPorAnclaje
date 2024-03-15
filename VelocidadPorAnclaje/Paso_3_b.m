% PRIMERA APROXIMACION A INTERPOLAR a mismas profundidades
% Basado en los trabajos de Kansow y Uwe Send
% viernes 5 de mayo 2023
%
% como parte de mis esfuerzos para calcular las auto-covarianzas espaciales
% Usando diferentes filtros temporales
%==========================================================================
% primera parte, hallar perfiles superficiales usando datos validos, para
% extrapolar
clear all; close all; clc
% break

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

cols = gray(9);
CNK = [29, 34, 37, 39, 42, 48, 57];
origen = 'I:\TesisDoc\datos\CNK_extractions\P3\A\';
%%
blanc = 1;
thresh = 100;

for an =  5 : 5%length(ancla)
    zpos_obs = {}; tpos_obs = {}; zpos_ins = {};tpos_ins = {};
    cou = 1; Vmat = []; tax = []; Umat = []; Vmask = [];
    
    figure('pos', [-1919 265 1536 739], 'color', 'w');
    axes('pos', [0.0606 0.3721 0.665311458333333 0.5694]);
    for kk = 1 : length(CNK)
        
        yys = []; xxs = []; ll = [];
        cnk = num2str(CNK(kk));
        fecs = datenum(char(cnkd(kk, :)));
        hold on;
        
        try
            clear anchor
            load([origen, 'CNK', cnk, '_', char(ancla(an))])
        catch
            
        end
        text( mean(fecs), 60, ['CNK-', cnk], 'HorizontalAlignment', 'center');
        if exist('anchor') && ~isempty(fieldnames(anchor))
            
            regected = logical([anchor(:).RegVel]');
            anchor(regected) = [];
            
            zpo = [anchor(:).znom];
            [zpo, ai] = sort(zpo, 'ascend');
            
            anchor = anchor(ai);
            td =[]; z = []; v = []; u = []; % por seguridad resetea variables principales
            
            name = anchor(1).name;
            tira = find(int16(name) == 45);
            tiran = str2num(name(tira(1)+2:tira(2)-1));
            
            if sum(seccion == 'Yucatan')  % si es Yucatan hay que tratar distinto a YUC1, 2 y 3
                
                if tiran >= 500
                    start = 70;
                    dz = 10;
                elseif an == 1
                    start = 5;
                    dz = 1;
                elseif an == 2
                    start = 10;
                    dz = 2;
                elseif an == 3
                    start = 20;
                    dz = 2;
                end
                
            else % si es florida
                start = 30;
                dz = 10;
            end
            
            dumx= [anchor(:).x];
            lon = nanmean(dumx(:));
            
            
            Zn = (-start:-dz:-2100)';
            tri = anchor(1).time;
            %
            for tti = 1:length(tri)
                yy = []; xx = []; gg = [];
                for i = 1:length(anchor)
                    t = anchor(i).time(tti);
                    x = anchor(i).x(tti,:);
                    z = -abs(anchor(i).z(tti,:));
                    p = anchor(i).p(tti,:);
                    v = anchor(i).v(tti,:);
                    u = anchor(i).u(tti,:);
                    
                    znom = anchor(i).znom;
                    dv = v(:); dz = -z(:); du = u(:);
                    
                    tax(cou) = t;
                    %                     dz(dz> -50) = NaN;
                    inba = isnan(dv+dz); dv(inba) =[]; dz(inba) = [];
                    if ~isempty(dv)
                        yy = [yy;v(:)]; gg = [gg; u(:)]; xx = [xx; z(:)];
                    end
                end
                if ~isempty(yy)
                    
                    [Zi,ia,idx] = unique(xx,'rows');
                    Zi = -abs(Zi);
                    
                    if tti == 200
                        zpos_obs{kk} = Zi(:);
                        tpos_obs{kk} = ones(size(Zi(:))).*fecs(1);
                    end
                    Vi = accumarray(idx,yy,[],@mean);
                    Ui = accumarray(idx,gg,[],@mean);
                    
                    for rr = 1 : length(Zn)
                        dip(rr) = min( abs(abs(Zn(rr)) - abs(Zi))  );
                    end
                    
                    if length(Vi) > 3
                        
                        Vi = interp1([Zi; -(tiran+100)], [Vi;0], Zn, 'pchip');
                        Ui = interp1([Zi; -(tiran+100)], [Ui;0], Zn, 'pchip');
                        
                        inba = find( Zn > roundn( max(Zi) , 1 ) |  Zn < roundn( min(Zi) , 1 ));
                        Vi(inba) = NaN;
                        Ui(inba) = NaN;
                        mask = Zn./Zn; mask(inba) = NaN;
                        if blanc
                            Vi(dip >= thresh) = NaN;
                            Ui(dip >= thresh) = NaN;
                        end
                    else
                        Vi = Zn*NaN;
                        Ui = Zn*NaN;
                    end
                else
                    Vi = Zn*NaN;
                    Ui = Zn*NaN;
                    mask = Zn*NaN;
                end
                
                Vmat(:, cou) = Vi;
                Umat(:, cou) = Ui;
                Vmask(:, cou) = mask;
                cou = cou + 1;
            end
        end
        [kk, length(CNK)]
    end
    %
    ti = tax(1):1/24:tax(end);
    undi = find(ismember(ti, tax));
    mask = nan([length(Zn), length(ti)]);
    mask(:, undi) = Vmask;
    Vmu = mask;
    Mag = sqrt(Vmat.^2 + Umat.^2);
    Vmu(:, undi) = smoothdata(Mag, 2, 'gaussian', 48);
    
    hold on;
    pcolor(ti(1:24:end), Zn, Vmu(:, 1:24:end));
    shading interp
    for kk = 1 : length(zpos_obs)
        zzzi = -abs(zpos_obs{kk});
        zzzi(zzzi > -start) = NaN;
        plot(tpos_obs{kk},zzzi, '.k', 'markersize', 9)
    end
    
    title(['(A).- ', char(ancla(an)), '. Magnitud de la velocidad horizontal'])
    xlim(datenum(['29-Jun-2012';'20-Jul-2023']))
    datetick('x', 'mmm-yyyy');
    xlim(datenum(['29-Jun-2012';'20-Jul-2023']))
    ylim([-(tiran+50), 100])
    
    box on;
    grid on
    cb = colorbar;
    ylabel(cb, 'Rapidez (m/s)', 'Interpreter', 'latex')
    xlabel('Tiempo (dias)', 'Interpreter', 'latex')
    ylabel('Profundidad (m)', 'Interpreter', 'latex')
    
    caxis([0, 1])
    set(gca,'TickLength', [1e-5, 1e-5],'pos', [0.0606 0.3721 0.675077083333333 0.5694]);
    
    axes('pos', [0.796895833333333 0.3721 0.1445 0.5694]);
    meprof = nanmean(Mag, 2);
    stprof = nanstd(Mag, [], 2);
    
    coy = nan(size(Zn));
    for liu = 1 : length(Zn)
        coy(liu) = length(find(~isnan(Vmask(liu, :))));
    end
    coy = coy < 365*24;
    ppr = nanmean(Mag, 2);
    
    ppr(coy) = NaN;
    meprof(coy) = NaN;
    
    yst = [meprof(:)+stprof(:); meprof(end:-1:1)-stprof(end:-1:1); meprof(1)+stprof(1)];
    zst = [Zn(:); Zn(end:-1:1); Zn(1)];
    
    inba = isnan(yst+ zst);
    yst(inba) = []; zst(inba) = [];
    
    hold on;
    pol = polyshape(yst, zst);
    pu = plot(pol); pu.EdgeColor = 'none';
    
    
    plot(ppr, -abs(Zn), 'k','linewidth', 2);
    
    ylim([-(tiran+50), 100])
    
    box on; set(gca, 'YAxisLocation', 'right')
    
    xlabel('Rapidez (m/s)')
    ylabel('Profundidad (m)')
    title([char(ancla(an)), '. Promedio en tiempo'])
    grid on;
    %=====================================================
    set(findall(gcf,'-property','Interpreter'),'Interpreter', 'latex');
    set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex');
    set(findall(gcf,'-property','FontSize'), 'FontSize',13);
    clear anchor
    drawnow
    
    set(findall(gcf,'-property','Interpreter'),'Interpreter', 'latex');
    set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex');
    colormap(parula(50));
%     if blanc
%         %         tosav = ['T:\Papers\Paper_I\Figures\Programas\ScalesEvaluation\', char(ancla(an)), '_bl', '.mat'];
%     else
%         %         tosav = ['T:\Papers\Paper_I\Figures\Programas\ScalesEvaluation\', char(ancla(an)), '.mat'];
%     end
    %         save( tosav, 'Umat', 'Vmat', 'Zn', 'tax', 'lon', 'Vmask', '-v7.3')
    %         export_fig(['T:\Papers\Paper_I\Figures\Programas\ScalesEvaluation\', char(ancla(an))])
    %     pause(0.5)
    %         close all;
    [an, length(ancla)]
end
%%
figure
plot(gradient(ppr), Zn)
hold on;
plot(ppr, Zn)
grid on
% figure;






%% plot((zpo*0)+736556-20, -zpo, 'ok')


%
% close all
% subplot(2, 1, 1)
% plot(tax, Vmat(20, :), '.')
% xlim(datenum(['05-Jul-2012';'10-Nov-2020']))
% datetick('x', 'mmm-yyyy');
% xlim(datenum(['05-Jul-2012';'10-Nov-2020']))
% % ylim([-(tiran+50), 50])
%
% subplot(2, 1, 2)
% plot(tax, (tax*0)+1, '.')
% xlim(datenum(['05-Jul-2012';'10-Nov-2020']))
% datetick('x', 'mmm-yyyy');
% xlim(datenum(['05-Jul-2012';'10-Nov-2020']))