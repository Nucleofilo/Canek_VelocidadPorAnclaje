% clear all; close all; clc
global seccion
% fechas de los caneks
cnkd = {['05-Jul-2012';'10-Jun-2013'];
    ['20-Jun-2013';'01-Jul-2014'];
    ['05-Jul-2014';'01-Jul-2015'];
    ['08-Jul-2015';'01-Aug-2016'];
    ['14-Aug-2016';'16-Jul-2018'];
    ['05-Aug-2018';'10-Nov-2020'];
    ['10-Aug-2022';'10-Jul-2023']};

% seccion = 'Yucatan';
% seccion = 'Florida';

switch seccion
    
    case 'Yucatan'
        ancla = {'YUC1', 'YUC2', 'YUC3', 'YUC4', 'YUC5', 'YUCI5', 'YUC6', 'YUCI6',...
            'YUC7', 'YUCI7', 'YUC8', 'YUCI8', 'YUC9', 'YUCI9', 'YUC10'};
    case 'Florida'
        ancla = {'EFL1', 'EFL2', 'EFL3', 'EFL4', 'EFL5', 'EFLI5', 'EFL6', 'EFL7'};
end
cols = gray(8);
CNK = [29  34	37  39	42  48  57];

origen = 'I:\TesisDoc\datos\CNK_extractions\P0\';

destino = 'I:\TesisDoc\datos\CNK_extractions\P1\FigPress\'; % donde guardara las figuras con
% las series de tiempo de presion

%%
for an = 1 : length(ancla) % corre para cada anclaje
    for kk = 1 : length(CNK) % corre para cada canek
        close all
        cnk = num2str(CNK(kk));
        try
            clear anchor
            load([origen, 'CNK', cnk, '_', char(ancla(an))])
            figure('pos', [30 3 1781 975], 'color', 'w');
        catch
            
        end
        
        if exist('anchor') && ~isempty(fieldnames(anchor))
            fecs = datenum(char(cnkd(kk, :)));
            
            zpo = [anchor(:).znom];
            flg = [anchor(:).Pflg];
            [zpo, ai] = sort(zpo, 'ascend');
            anchor = anchor(ai);
            
            flg = flg(ai);
            
            anchor = anchor(logical(flg));
            
            
            for i = 1 : length(anchor)
                
                t = anchor(i).time;
                p = anchor(i).p;
                znom = anchor(i).znom;
                subplot(length(anchor), 1, i);
                plot(t, p);
                datetick('x', 'mmm/yyyy');
                xlim([fecs])
                title([anchor(i).name, '. Pnom = ', num2str(znom)])
                drawnow
            end
        export_fig( [destino, char(ancla(an)),  '_','CNK',cnk] )
        pause(0.3)
        end
    end
end



