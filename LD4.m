%% 4 Laboratorinis darbas
% Vaizdo atpazinimas naudojant dirbtiniu neuronu tinkla
% Atliko KTFM21 grupes studentas - Vsevolod Kapustin
clear all
clc

pavadinimas = 'pav.jpg';
pozymiai_tinklo_mokymui = pozymiai_raidems_atpazinti(pavadinimas, 1);
%% Atpaþintuvo kûrimas
%% Development of character recognizer
% poþymiai ið celiø masyvo perkeliami á matricà
% take the features from cell-type variable and save into a matrix-type variable
x = cell2mat(pozymiai_tinklo_mokymui);

% Gaunamas masyvas 9 elementu su 35 pozymiais, reiskia turime 35 iejimus


%% RBF sluoksnio skaiciavimas
% 1 neurono isejimas (iejimas skaiciuojamas be svoriu)
% neuronu parametrai
% Spindulio parametras ir centro parametras(sudaromas masyvas is tiek spinduliu kiek yra
% pozymiu), Sio atveju 35

for n = 1 : length(x(:, 1))
    k = 0;
    c(n) = randn(1);
    r(n) = randn(1);
    %% Skaiciojama spindulio bazes funkcija kiekvienam iejimui ir
% kiekvienam pozymiui
    c = c(n);
    r = r(n);
    for k = 1 : length(x(n, :))
        f(n, k) = rbf_aktyvavimas(x(n, k), c, r);
         %% Svoriu generavimas sinapsems iki isejimo neurono
         w(n, k) = randn(1);
         b(k) = randn(1);
        
    end    
end 


%% Sinapiniu svoriu skaiciavimas kiekvienam pozymiui
%% Sinaptiniu svoriu turi buti tiek pat kiek raidziu
for k = 1 : length(f(1, :))
    % Vienai raidei parametrai
          funk = f(:, k);
          wunk = w(:, k);
          bunk = b;
    for n = 1 : length(f(1, :))
       if n == 1
          v(k) = funk(1)*wunk(1);
          continue
       end
       v(k) = v(k) + funk(n)*wunk(n);
    end
   % paskutinio pozymio atveju prie svoriu pridedamas "bias"
   v(k) = v(k) + bunk(k);
end

%% isejimo pozymiu vekoriu formavimas
% Isejimu verciu logika
% 0,ASCII(DEC)
%ASCII    V  S  E   V  O   L   O   D .
     y = [.86,.83,.45,.86,.415,.413,.415,.68,0];

%% Isejimo sluoksnio skaiciavimas
% Generuojami atsitiktiniai svoriai ir neuronu pralaidos koeficientai
% Isejimo neurono aktyvavimo funkcija -- tiesine
yy = v;


%% Klaidos skaiciavimas
e = y - yy

e(1)

%% Perceptrono mokymas, vykdomas kiekvienai isejimo vertei atskirai

    miu = 0.00001;
 

    
%% Dydziai mokymui viena raide
y = y(1);
x = x(:, 1);
w = w(:, 1);
b = b(1);
c = c(:, 1);
r = r(:, 1);
e = e(1);

while 1
   %% Mokymas viena raide
    for n = 1 : length(x) %% Kiekvieno iejimo stebejimas
        %% Atnaujiname svorius
        w(n) = w(n) + miu*e*x(n);
    end
    b = b + miu*b;
    for n = 1 : length(x) %% Kiekvieno iejimo stebejimas
        %% RBF aktyvavimas
        f(n) = rbf_aktyvavimas(x(n), c, r);
    end
    
    for n = 1 : length(x) %% Kiekvieno iejimo stebejimas
        %% Sinaptiniu svorio skaiciavimas
        if n == 1
           v = w(n)*f(n); 
        end
        v = v + w(n)*f(n);
    end
    v = v + b;
    
    %% Tisinis isejimas
    yy = v;
    
    %% Kalidos skaiciavimas
    e = y - yy
    
        if e < 1e-4;
            if e > -1e-4
                break
            end
        end
    
end
    
%% Testavimas
    for n = 1 : length(x) %% Kiekvieno iejimo stebejimas
        %% RBF aktyvavimas
        f(n) = rbf_aktyvavimas(x(n), c, r);
    end
    
    for n = 1 : length(x) %% Kiekvieno iejimo stebejimas
        %% Sinaptiniu svorio skaiciavimas
        if n == 1
           v = w(n)*f(n); 
        end
        v = v + w(n)*f(n);
    end
    v = v + b;
    
    %% Tisinis isejimas
    yy = v;
   
    yy
    
    if yy > 0.8556 %% Nepatenka i ASCII DEC=85 sryti
       if yy < 0.8649 %% Nepatenka i ASCII DEC=87 sryti
          fprintf("Pirma raide yra V\n") 
       end
    else
         fprintf("Kazkas ne taip, viskas blogai!!!")
    end
    
%% Programos galas

%% Spindulio tipo funkcija
function f = rbf_aktyvavimas(x, c, r)
% x - iejimo vektorius
% c - centro parametras
% r - spindulys
    f = exp(-(x-c).^2/(2*r.^2));
end


function pozymiai = pozymiai_raidems_atpazinti(pave, eil)
   
    % nuskaitomas ir atvaizduojamas paveiksliukas su ranka rasytais
    % skaiciais
    pav = imread(pave);
    
    % Iskirpiamos raides ir sudeliojimas i kintamojo celes
    pav_pustonis = im2gray(pav);
    figure(1);imshow(pav_pustonis);
    
    % Vaizdo keitimu dvejetaine slenkties reiksmes paieska
    slenkstis = graythresh(pav_pustonis);
    
    % Pustonio vaizdo keitimas svejetainiu
    pav_bin = im2bw(pav_pustonis, slenkstis);
    imshow(pav_bin);
    
    % Vaizde esanciu konturu paieska
    pav_konturais = edge(uint8(pav_bin));
    figure(2);imshow(pav_konturais);
    
    % Objektu konturu uzpildymas
    se = strel('square', 7);
    pav_uzpildytas = imdilate(pav_konturais, se)
    figure(3);imshow(pav_uzpildytas);
    
    % Tustumu objekto viduje uzpildymas
    pav_vientisi = imfill(pav_uzpildytas,'holes');
    figure(4); imshow(pav_vientisi);
    
    % Vientesu objektu numeravimas dvejatainiame vaizde
    [O_suzymeti Skaicius] = bwlabel(pav_vientisi)
    
    % Apskaiciojami pozymiai
    O_pozymiai = regionprops(O_suzymeti);
    
    % Nuskaitomos pozymiu - objektu ribu koordinaciu - reiksmes
    O_ribos = [O_pozymiai.BoundingBox];
    
    %kintamasis skaicius - objektu skaicius
    
    O_ribos = reshape(O_ribos,[4 Skaicius]); % Skaicius - objektø skaièius
% nuskaitomos poþymiø - objektø masës centro koordinaèiø - reikðmës
% reag the mass center coordinate
O_centras = [O_pozymiai.Centroid];
% kadangi centrà nusako 2 koordinatës, pergrupuojame reikðmes
% group center coordinate values
O_centras = reshape(O_centras,[2 Skaicius]);
O_centras = O_centras';
% pridedamas kiekvienam objektui vaize numeris (treèias stulpelis ðalia koordinaèiø)
% set the label/number for each object in the image
O_centras(:,3) = 1:Skaicius;
% surûðiojami objektai pagal x koordinatæ - stulpelá
% arrange objects according to the column number
O_centras = sortrows(O_centras,2);
% rûðiojama atsiþvelgiant á pavyzdþiø eiluèiø ir raidþiø skaièiø
% sort accordign to the number of rows and number of symbols in the row
raidziu_sk = Skaicius/eil;
for k = 1:eil
    O_centras((k-1)*raidziu_sk+1:k*raidziu_sk,:) = ...
        sortrows(O_centras((k-1)*raidziu_sk+1:k*raidziu_sk,:),3);
end
% ið dvejetainio vaizdo pagal objektø ribas iðkerpami vaizdo fragmentai
% cut the symbol from initial image according to the bounding box estimated in binary image
for k = 1:Skaicius
    objektai{k} = imcrop(pav_bin,O_ribos(:,O_centras(k,3)));
end
% vieno ið vaizdo fragmentø atvaizdavimas
% show one of the symbol's image
figure(5),
for k = 1:Skaicius
    subplot(eil,raidziu_sk,k), imshow(objektai{k})
end
% vaizdo fragmentai apkerpami, panaikinant fonà ið kraðtø (pagal staèiakampá)
% image segments are cutt off
for k = 1:Skaicius % Skaicius = 88, jei yra 88 raidës
    V_fragmentas = objektai{k};
    % nustatomas kiekvieno vaizdo fragmento dydis
    % estimate the size of each segment
    [aukstis, plotis] = size(V_fragmentas);
    
    % 1. Baltø stulpeliø naikinimas
    % eliminate white spaces
    % apskaièiuokime kiekvieno stulpelio sumà
    stulpeliu_sumos = sum(V_fragmentas,1);
    % naikiname tuos stulpelius, kur suma lygi aukðèiui
    V_fragmentas(:,stulpeliu_sumos == aukstis) = [];
    % perskaièiuojamas objekto dydis
    [aukstis, plotis] = size(V_fragmentas);
    % 2. Baltø eiluèiø naikinimas
    % apskaièiuokime kiekvienos seilutës sumà
    eiluciu_sumos = sum(V_fragmentas,2);
    % naikiname tas eilutes, kur suma lygi ploèiui
    V_fragmentas(eiluciu_sumos == plotis,:) = [];
    objektai{k}=V_fragmentas;% áraðome vietoje neapkarpyto
end
% vieno ið vaizdo fragmentø atvaizdavimas
% show the segment
figure(6),
for k = 1:Skaicius
    subplot(eil,raidziu_sk,k), imshow(objektai{k})

end
%%
%% Suvienodiname vaizdo fragmentø dydþius iki 70x50
%% Make all segments of the same size 70x50
for k=1:Skaicius
    V_fragmentas=objektai{k};
    V_fragmentas_7050=imresize(V_fragmentas,[70,50]);
    % padalinkime vaizdo fragmentà á 10x10 dydþio dalis
    % divide each image into 10x10 size segments
    for m=1:7
        for n=1:5
            % apskaièiuokime kiekvienos dalies vidutiná ðviesumà 
            % calculate an average intensity for each 10x10 segment
            Vid_sviesumas_eilutese=sum(V_fragmentas_7050((m*10-9:m*10),(n*10-9:n*10)));
            Vid_sviesumas((m-1)*5+n)=sum(Vid_sviesumas_eilutese);
        end
    end
    % 10x10 dydþio dalyje maksimali ðviesumo galima reikðmë yra 100
    % normuokime ðviesumo reikðmes intervale [0, 1]
    % perform normalization
    Vid_sviesumas = ((100-Vid_sviesumas)/100);
    % rezultatà (poþmius) neuronø tinklui patogiau pateikti stulpeliu
    % transform features into column-vector
    Vid_sviesumas = Vid_sviesumas(:);
    % iðsaugome apskaièiuotus poþymius á bendrà kintamàjá
    % save all fratures into single variable
    pozymiai{k} = Vid_sviesumas;
end 
end
