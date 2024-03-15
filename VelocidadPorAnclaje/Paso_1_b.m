% clear all; close all; clc
global seccion
CNK = [29  34	37  39	42  48  57];

cnkd = {['05-Jul-2012';'10-Jun-2013'];
    ['20-Jun-2013';'01-Jul-2014'];
    ['05-Jul-2014';'01-Jul-2015'];
    ['08-Jul-2015';'01-Aug-2016'];
    ['14-Aug-2016';'16-Jul-2018'];
    ['05-Aug-2018';'10-Nov-2020'];
    ['10-Aug-2022';'10-Jul-2023']};

origen = 'I:\TesisDoc\datos\CNK_extractions\P0\';
destinof = 'I:\TesisDoc\datos\CNK_extractions\P1\FigPress\';
destinoa = 'I:\TesisDoc\datos\CNK_extractions\P1\F1\';

Ta = readtable(['I:\TesisDoc\datos\CNK_extractions\P1\FaPos.xlsx']);

secc = Ta.Seccion;
anc = Ta.Anclaje;
cnk = (Ta.CNK); 
NS = Ta.NS;
%
for k = 1:length(NS)  
    arch = [origen, 'CNK',num2str(cnk(k)), '_', char(  anc(k) )];
    
    load(arch);
    
    C = {anchor.name}';
    inst = find(contains(C, NS(k)));
    anchor(inst).Pflg = 0;
    anchor(inst).RPflg = 1; % Discarded Pressure Flag
    save(arch, 'anchor', '-v7.3')
    anc(k)
end

%%
for k = 1:length(NS)  
    close all;
    figure('pos', [30 3 1781 975], 'color', 'w');
    arch = [origen, 'CNK',num2str(cnk(k)), '_', char(  anc(k) )];

    load(arch)
    
%  fecs = datenum(char(cnkd(k, :)));
            
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
                xlim([min(t), max(t)])
                title([anchor(i).name, '. Pnom = ', num2str(znom)])
 
            end
        export_fig( [destinof, char(  anc(k) ),  '_','CNK',num2str(cnk(k)), 'rp'] )
        pause(0.3)

end

