%% CALCOLO DELLA MASSA RILASCIATA FIGURA 3 per il DNR
%Massa rilasciata

steps=51;
[MrelRph7, tRph7] = mass_released("datiR_2701", steps);

steps=76;
[MrelRph5, tRph5] = mass_released("datiR_2701ph5", steps);

steps=101;
[MrelRph4, tRph4] = mass_released("datiR_2701ph4", steps);

% Riscalamento del tempo per ottenere le curve della figura 
tRph7 = tRph7 * 100/0.66;
tRph5 = tRph5 * 22/1.08;
tRph4 = tRph4 * 12/1.88;

close("all")
plot(tRph7, MrelRph7, color = "g")
hold on
plot([tRph5; 75], [MrelRph5, MrelRph5(end)], color = "r")
plot([tRph4; 75], [MrelRph4, MrelRph4(end)], color = "magenta")

ylim([0, 0.6])
xlim([0, 75])
grid on
legend("ph7.4", "ph5.5", "ph4")

%% CALCOLO DELLA MASSA RILASCIATA FIGURA 4 per il CTP
steps=101; % numero di istanti di tempo 

[MrelRph7, tRph7] = mass_released("datiR2_2701ph7", steps);
[MrelRph4, tRph4] = mass_released("datiR2_2701ph4", steps);

tRph7 = tRph7 * 10 / 37.5;
tRph4 = tRph4 * 5 / 9.2;
plot([tRph7; 25], [MrelRph7, MrelRph7(end)], color = "g")
hold on
plot([tRph4; 25], [MrelRph4, MrelRph4(end)], color = "magenta")
xlim([0, 25])
ylim([0, 0.6])
grid on

%% FITTING a PH variabile
% Inserimento dei DATI del paper

% DNR
ph = [7.4, 5.5, 4.0];
D = [113.05, 61.50, 35.60];
Da = [1.10, 0.96, 0.63];
S = [5.60, 3.35, 2.64];
k = [10.89, 4.65, 0.47];

% CTP
ph2 = [7.4, 4.0];
D2 = [86.91, 48.59];
Da2 = [0.67, 0.14];
S2 = [0.11, 2.20];
k2 = [6.31, 0.09];

% Fitting per il  DNR
coefD = polyfit(ph, D, 2);
coefDa =polyfit(ph, Da, 2);
coefS = polyfit(ph, S, 2);
coefk = polyfit(ph, k, 2);


%Fitting per il CTP

coefD2 = polyfit(ph2, D2, 1);
coefDa2 =polyfit(ph2, Da2, 1);
coefS2 = polyfit(ph2, S2, 1);
coefk2 = polyfit(ph2, k2, 1);

% Estrazione dei valori per le simulazioni a diversi pH
cases = 18;
phval = linspace(4, 7.4, cases);
Dval = polyval(coefD, phval);
Daval = polyval(coefDa, phval);
Sval = polyval(coefS, phval);
kval = polyval(coefk, phval);

D2val = polyval(coefD2, phval);
Da2val = polyval(coefDa2, phval);
S2val = polyval(coefS2, phval);
k2val = polyval(coefk2, phval);

%% MASSA RILASCIATA A DIVERSI PH
% DNR
Mrel = mass_pH("D_DNR", cases);
plot(phval, Mrel)
hold on
% CTP
Mrel = mass_pH("D_CTP", cases);
plot(phval, Mrel)
legend("DNR", "CTP")
%%
function [Mrel, t] = mass_released(name, steps)

    % Impostando le opzione per l'inserimento dei valori di COMSOL 

    opts_1 = delimitedTextImportOptions("NumVariables", 5);
    opts_1.DataLines = [2, Inf];
    opts_1.Delimiter = " ";
    opts_1.VariableNames = ["R", "t", "c", "b_0", "b_1"];
    opts_1.VariableTypes = ["double", "double", "double", "double", "double"];
    opts_1.ExtraColumnsRule = "ignore";
    opts_1.EmptyLineRule = "read";
    opts_1.ConsecutiveDelimitersRule = "join";
    opts_1.LeadingDelimitersRule = "ignore";
    
    % Importazione dei dati

    tbl = readtable("/Users/andrearosini/Politecnico Di Torino Studenti Dropbox/Andrea Rosini/modelli fisio/progetto/" + name + "clean.csv", opts_1);
    R = tbl.R;
    t = tbl.t;
    c = tbl.c;
    b_0 = tbl.b_0;
    b_1 = tbl.b_1;
    
    % Creazione vettori di tempo e spazio per i domini

    T = steps;
    t = t(1:T);
    R = R(1:T:end);
    L = length(R);
    R_0 = R(1:L/2);
    R_1 = R(L/2 + 1: end);
    
    % Formattazione dei valori delle concentrazioni

    c = reshape(c, [T,L]);
    b_0 = reshape(b_0, [T,L]);
    b_1 = reshape(b_1, [T, L]);
    b_0 = b_0(:, 1:L/2);
    b_1 = b_1(:, L/2+1:end);
    c_0 = c(:, 1:L/2);
    c_1 = c(:, L/2+1:end);
    
    % Calcolo dell'integranda per la massa

    b_0_Integrand = R_0 .^ 2 .* b_0';
    b_1_Integrand = R_1 .^ 2 .* b_1';
    c_0_Integrand = R_0 .^ 2 .* c_0';
    c_1_Integrand = R_1 .^ 2 .* c_1';
    
    % Integrazione con il metodo del trapezio
   
    B_0 = trapz(R_0, b_0_Integrand)* 4 * pi;
    B_1 = trapz(R_1, b_1_Integrand)* 4 * pi;
    C_0 = trapz(R_0, c_0_Integrand) * 4 * pi;
    C_1 = trapz(R_1, c_1_Integrand) * 4 * pi;
    
    % Calcolo della massa rilasciata
    
    Mtot = B_0 + B_1 + C_0 + C_1;
    M0 = (B_0(1) + B_1(1)) * ones(size(Mtot));
    Mrel = (M0 - Mtot) ./ M0;
    
end

%% Set up the Import Options and import the data
function [Mrel] = mass_pH(name, ~)

opts = delimitedTextImportOptions("NumVariables", 179);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ["\t", " "];

% Specify column names and types
opts.VariableNames = ["R", "D", "Da", "S", "k", "t", "c", "b_0", "b_1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";


% Import the data
    tbl = readtable("/Users/andrearosini/Politecnico Di Torino Studenti Dropbox/Andrea Rosini/modelli fisio/progetto/" + name + "clean.csv", opts);
    R = tbl.R;
    
    c = tbl.c;
    b_0 = tbl.b_0;
    b_1 = tbl.b_1;
    cases = 18;
    L = 102;
    R = R(1:L);
    R_0 = R(1:L/2);
    R_1 = R(L/2 + 1: end);
    
    c = reshape(c, [L, cases]);
    b_0 = reshape(b_0, [L,cases]);
    b_1 = reshape(b_1, [L, cases]);
    b_0 = b_0(1:L/2, :);
    b_1 = b_1(L/2+1:end, :);
    c_0 = c(1:L/2, :);
    c_1 = c(L/2+1:end, :);
    
    b_0_Integrand = R_0 .^ 2 .* b_0;
    b_1_Integrand = R_1 .^ 2 .* b_1;
    c_0_Integrand = R_0 .^ 2 .* c_0;
    c_1_Integrand = R_1 .^ 2 .* c_1;
    
    B_0 = trapz(R_0, b_0_Integrand)* 4 * pi;
    B_1 = trapz(R_1, b_1_Integrand)* 4 * pi;
    C_0 = trapz(R_0, c_0_Integrand) * 4 * pi;
    C_1 = trapz(R_1, c_1_Integrand) * 4 * pi;
    
    
        

    Mtot = B_0 + B_1 + C_0 + C_1;
    
    if name == "D_CTP"
        M0 = 0.958;
        Mtot = Mtot(end:-1:1); %I dati sono stati presi in ordine storto per il CTP
    else
        M0 = 1.1834;
    end

    M0 = M0 * ones(size(Mtot));
    Mrel = (M0 - Mtot) ./ M0;


end
