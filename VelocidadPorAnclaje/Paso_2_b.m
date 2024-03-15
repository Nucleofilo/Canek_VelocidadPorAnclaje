% clear all; close all; clc

addpath(genpath('C:\Users\nucle\Documents\GitHub\Canek_VelocidadPorAnclaje\extrn\'));

CNK = [29  34	37  39	42  48  57];
cnkd = {['05-Jul-2012';'10-Jun-2013'];
    ['20-Jun-2013';'01-Jul-2014'];
    ['05-Jul-2014';'01-Jul-2015'];
    ['08-Jul-2015';'01-Aug-2016'];
    ['14-Aug-2016';'16-Jul-2018'];
    ['05-Aug-2018';'10-Nov-2020']
    ['10-Aug-2022';'10-Jul-2023']};


origen = 'I:\TesisDoc\datos\CNK_extractions\P0\';
destino = 'I:\TesisDoc\datos\CNK_extractions\P2\FigDDPress\';
Ta = readtable(['I:\TesisDoc\datos\CNK_extractions\/P2/Drift.xlsx']);

anc = Ta.Anclaje;
cnk = (Ta.CNK);
NS = Ta.NS;
Ty = Ta.Tipo;

% exponencial lineal
f3 = fittype('a_1*(1-exp(a_2.*x)) + a_3.*x + a_4');
% exponencial
f2_a = fittype('b_1*exp(b_2.*x)');
% exponencial con shifting
f2_b = fittype('b_1*exp(b_2.*x) + b_3');

%%
%%
for k = 1:length(NS)
    
    arch = [origen, 'CNK',num2str(cnk(k)), '_', char(  anc(k) )];
    
    load(arch);
    
    C = {anchor.name}';
    inst = find(contains(C, {['-', char(NS(k)), '-']}));
    
    %============================================
    t = anchor(inst).time;
    p = anchor(inst).p;
    znom = anchor(inst).znom;
    %             subplot(n, 1, a);
    %=====================
    y_ = p-nanmean(p);
    x_ = t;
    x_ = x_ - min(x_);
    inba = isnan(y_ + x_);
    y_(inba) = []; x_(inba) = [];
    %============================================
    
    tipo = Ty(k);
    figure('pos', [30 3 1781 975], 'color', 'w');
    subplot(2, 1, 1)
    plot(x_, y_, 'k')
    xlim([min(x_), max(x_)]);
    title(anchor(inst).name)
    switch tipo
        
        case 1
            
            Bo = polyfit(x_, y_, 1);
            hold on;
            rf = refline(Bo); rf.LineWidth = 2;
            rf.Color = 'r';
            deriva = polyval(Bo, t(:)-t(1));

        case 2
            % solo exponencial pero para ver cual se ajusta mejor se hace
            % el caso con shifting tambien
            [c1,~] = fit(x_, smoothdata(y_, 'movmean', 20*24), f2_a, 'StartPoint',[10, -1e-4]);
            A1 = coeffvalues(c1);
            
            trende = A1(1)*(exp(A1(2).*x_));
            
            hold on;
            plot(x_, trende, 'r', 'linewidth', 2);
            
            [c2,~] = fit(x_, smoothdata(y_, 'movmean', 20*24), f2_b, 'StartPoint',[A1, 0]);
            A2 = coeffvalues(c2);
            
            trende = A2(1)*(exp(A2(2).*x_)) + A2(3);
            plot(x_, trende, 'r--', 'linewidth', 2);
            
%             fab = ...
%                 input...
%                 ('Indique cual de los dos ajustes es más adecualdo: \n 1: linea continua o 2: linea quebrada '...
%                 );
fab = 2;
            
            switch fab 
                case 1
                    deriva = A1(1)*(exp(A1(2).*(t(:)-t(1)))  );
                case 2
                    deriva = A2(1)*(exp(A2(2).*(t(:)-t(1)))  ) + A2(3);
            end   
        case 3
            % exponencial lineal, es medio un reto encontrar un buen set de
            % parametros iniciales
            % para dar parametro inicial de la parte exponencial:
                [c1,~] = fit(x_, smoothdata(y_, 'movmean', 20*24),f2_a,'StartPoint',[0.1, -1e-4]);
                fge = coeffvalues(c1);
                
                Bo = polyfit(x_, y_, 1);
                
                try
                [c2,gof] = fit(x_, smoothdata(y_, 'movmean', 20*24), f3,'StartPoint',[fge(1:2) 0 mean(y_)]);
                A = coeffvalues(c2);
                trende = A(1)*(1-exp(A(2).*x_)) + A(3).*x_ + A(4);
                hold on;
                plot(x_, trende, 'r', 'linewidth', 2)
                catch
                [c2,gof] = fit(x_, smoothdata(y_, 'movmean', 20*24), f3,'StartPoint',[1, -0.0001 Bo(1) mean(y_)]);
                A = coeffvalues(c2);
                trende = A(1)*(1-exp(A(2).*x_)) + A(3).*x_ + A(4);
                hold on;
                plot(x_, trende, 'r', 'linewidth', 2)
                end
                
                deriva =  A(1)*(1-exp(A(2).*(t(:)-t(1)))) + A(3).*(t(:)-t(1)) + A(4);
                        
    end
    
    detre = p(:) - deriva(:);
    
    
    subplot(2, 1, 2)
    plot(t, detre);
    xlim([min(t), max(t)]);
    datetick('x', 'mmm/yyyy')
    xlim([min(t), max(t)]);
    title('De-drifted pressure series')
    drawnow
    pause(0.1);
%     export_fig( [destino char(  anc(k) ),  '_','CNK',num2str(cnk(k)), 'ddp'] )
%     pause()
    pause(0.1);
    close all;
    anchor(inst).p = detre;
    anchor(inst).DPflg = 1;
    
    save(arch, 'anchor', '-v7.3')
end
