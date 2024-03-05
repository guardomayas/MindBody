%%% Análisis tasa de descuento estudio mente y cuerpo %%%
% silvia 07.10.2020
% santiago 08.29.2023

% Specifying Data Pathways with MATLAB Drive. I will modify it using only
% relative pathways. 
data_folder = 'data'; %Specify with your own folder
data_file_folders = dir(fullfile(data_folder,'datos*'));
sprintf('Tenemos %d carpetas de datos', length(data_file_folders))

% extraer datos de análisis de cada sesión que son csv 
for f = 1:length(data_file_folders) %Uncomment length when you want to analyze pavlovia data
    % esta_sesion = fullfile(data_folder,data_file_folders(f).name,'convert');
    esta_sesion = fullfile(data_folder,data_file_folders(f).name);
    addpath(esta_sesion)
    prep = dir(fullfile(esta_sesion,'*.csv')); %Check for the file that crashes here 
    %with prep(d).name
    for d = 1:length(prep)
        getBDataEstudiantes(prep(d).name)
    end
end 
% extraer datos de análisis de cada sesión que son xlsx 
for f = 1:length(data_file_folders) %Uncomment length when you want to analyze pavlovia data
    % esta_sesion = fullfile(data_folder,data_file_folders(f).name,'convert');
    esta_sesion = fullfile(data_folder,data_file_folders(f).name);
    addpath(esta_sesion)
    prep = dir(fullfile(esta_sesion,'*.xlsx'));
    for d = 1:length(prep)
        getBDataEstudiantes(prep(d).name)
    end
end 
% correr algoritmo de optimización para parámetros tasa de descuento en
% todos los archivos

exData_folder = fullfile(pwd,'bData');
addpath(exData_folder)
exFiles = dir(fullfile(exData_folder,'BE*'));

for s = 1:length(exFiles)
    esteParticipante = exFiles(s).name;
    ITClh(esteParticipante)
end

% organizar los datos relevantes en una tabla

parametros = nan(length(exFiles),7);
param_folder = fullfile(pwd,'estLH');
addpath(param_folder)
pFiles = dir(fullfile(param_folder,'LH*'));

for s = 1:length(exFiles)
    load(pFiles(s).name)
    parametros(s,1) = subject;
    parametros(s,2) = session;
    parametros(s,3) = b(2);
    parametros(s,4) = log(b(2));
    parametros(s,5) = b(1);
    parametros(s,6) = r2;
    parametros(s,7) = percentimp;
end

% organizar por orden de número de participante
parametros = sortrows(parametros,1);

% crear una tabla y guardarla
tablaParam = array2table(parametros);
tablaParam.Properties.VariableNames = {'Num' 'Sesion' 'kappa' 'ln_kappa' 'beta' 'r2' 'percentImp'};

save('parametrosMC.mat','tablaParam')

% importar tabla general de datos mente y cuerpo
% addpath('/Users/silvia/Desktop')
T = readtable('BDmente.xlsx');
T = T(find(~isnan(T.num)), :);

T.kappa = nan(height(T),1);
T.ln_kappa = nan(height(T),1);
T.beta = nan(height(T),1);
T.r2 = nan(height(T),1);
T.percentImp = nan(height(T),1);

% adicionar datos parametros
%Tab = join(tablaParam,T,'keys','Num');
for p = 1:height(T)
    i = find(T.num(p)==tablaParam.Num);
    if isempty(i) == 0 
        T.kappa(p) = tablaParam.kappa(i);
        T.ln_kappa(p) = tablaParam.ln_kappa(i);
        T.beta(p) = tablaParam.beta(i);
        T.r2(p) = tablaParam.r2(i);
        T.percentImp(p) = tablaParam.percentImp(i);
    else
        T.kappa(p) = nan;
        T.ln_kappa(p) = nan;
        T.beta(p) = nan;
        T.r2(p) = nan;
        T.percentImp(p) = nan;
    end
end

writetable(T,'datosMCmasDescuento_20220224.xlsx','FileType','spreadsheet')
%% Resultados
% distribucion de tasa de descuento y media
figure
%h = histogram(T.ln_kappa,20);
h = histogram(T.kappa,20);
h.Normalization = 'probability';
rgb1 = [82 117 181]./255; 
h.FaceColor = rgb1;
h.EdgeColor = 'w';
h.FaceAlpha = 0.8;
% xlabel('tasa de descuento (ln(\kappa))')
xlabel('tasa de descuento (\kappa)')
ylabel('proporción')
set(gca,'TickDir','Out','FontSize',16,'FontWeight','Normal')
box off
%set(gca,'YLim',[0 0.15])
% hold on
% p1 = plot(mean(T.ln_kappa,'omitnan'),0.13,'v','Color',rgb1);
% p1.MarkerSize = 14;
% p1.MarkerFaceColor = rgb1;
% t = text(-2,0.13,['promedio = ' sprintf('%0.2f',mean(T.ln_kappa,'omitnan'))]);
% t.FontSize = 20;


% Análisis multivariado (modelo lineal generalizado)
model_Dem = fitglm(T, 'ln_kappa ~ 1 + edad + sexo');
display(model_Dem)

model_Adict = fitglm(T, 'ln_kappa ~ 1 + tabaco + alcohol'); %queda pendiente recodificar la variable del alcohol
display(model_Adict)

model_Metab = fitglm(T, 'ln_kappa ~ 1 + imc + colesterol + glucemia');
display(model_Metab)

model_Psicol = fitglm(T, 'ln_kappa ~ 1 + estrs + iritotal + atencinemocional + claridaddesentimientos + reparacinemocional');
display(model_Psicol)

model_Cardio = fitglm(T, 'ln_kappa ~ 1 + scvscore');
display(model_Cardio)

model_final = fitglm(T, 'ln_kappa ~ 1 + inteligencia_emocional_tmms24_co * estrs + edad + imc + colesterol + glucemia + iritotal');
display(model_final)

[r,p] = corr(T.estrs,T.inteligencia_emocional_tmms24_co,'rows','complete')
%%
s = scatter(T.fc,T.ln_kappa)
