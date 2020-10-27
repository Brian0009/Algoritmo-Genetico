clc
clear
%ALGORITMO GENETIVO V1.02
%Brian Giovanny Alfonso Rodriguez 20151020600
%----------------------------------------------------------------------------------------------------
%Variables x y
syms x y w r sigma c
%definimos la funcion objetivo Z%
%z = 2*(  180*(sqrt((6-x)^2+(12-y)^2)) + 140*(sqrt((11-x)^2+(7-y)^2)) + 250*(sqrt((x)^2+(8-y)^2)) + 310*(sqrt((9-x)^2+(2-y)^2)) );
%z =  ((0.2*x^2)+(0.08*y^2)+(0.18*w^2)+(0.1*x*y)+(0.04*x*w)+(0.06*y*w))+(r*((-0.14*x)-(0.11*y)-(0.1*w+120))+c*(abs(x+y+w-1000)-sigma));
z =  (0.2*x^2)+(0.08*y^2)+(0.18*w^2)+(0.1*x*y)+(0.04*x*w)+(0.06*y*w)+(r*(-0.14*x-0.11*y-0.1*w+120)+c*abs(x+y+w-1000)-sigma);
%fi = {@(x,y,z) (0.2*x^2)+(0.08*y^2)+(0.18*z^2)+(0.1*x*y)+(0.04*x*z)+(0.06*y*z)+(r*(-0.14*x-0.11*y-0.1*z+120)+c*abs(x+y+z-1000)-sigma)};

%otras variables:
numCrom = 100 ;         %Numero poblacón
nMinValCrom = 0;            %Minimo valor de cromosma
nMaxValCrom = 500;           %Maximo valor de cromosma
r = 80000;
c = 50000;
Pc =0.3;                    %Probabilidad de cruce
Pm =0.2;                   %Probabilidad de mutacion
sigma = 10;
nIteraciones = 50;
%phii= eval(z)
%Definimos los cromosomas
cromosoma=[numCrom,3];
i = 1;      %iterador
%valores aleatorios de los cromosmomas
while i <numCrom+1
    random=aleatorios();
    cromosoma(i,1)= random(1);
    cromosoma(i,2)= random(2);
    cromosoma(i,3)= random(3);
    i= i+1;
    
end
%{
   cromosoma(4,1)=380.952;
   cromosoma(4,2)=476.1904;
   cromosoma(4,3)=142.8571;
%}

GmejorCormosoma=[0,0,0];
GZauxMin=9999999999;
%Mapeado del ejemplo
%{
cromosoma(1,1)=8; cromosoma(1,2)=10;
cromosoma(2,1)=13; cromosoma(2,2)=2;
cromosoma(3,1)=5; cromosoma(3,2)=11;
cromosoma(4,1)=7; cromosoma(4,2)=1 ;
cromosoma(5,1)=0; cromosoma(5,2)=1 ;
cromosoma(6,1)=2; cromosoma(6,2)=10 ;
cromosoma(7,1)=10; cromosoma(7,2)=5;
cromosoma(8,1)=5; cromosoma(8,2)=10;
cromosoma(9,1)=14; cromosoma(9,2)=11;
cromosoma(10,1)=1; cromosoma(10,2)=9;
%}
it = 0;
while it < nIteraciones
    
    %---------------EVALUCACIÔN------------------------------------------
    i = 1;          %iterador
    evalCrom=[numCrom];  %funcion objetivo evaluada
    %se evalua la funcion con los cromosamas obtenidos
    
    x=cromosoma(1,1);
    y=cromosoma(1,2);
    w=cromosoma(1,3);
    peorCromosoma=0;
    mejorCromosoma=1;
    zAuxMin =99999999;
    zAuxMax = 0;
    while i < numCrom +1
        
        x=cromosoma(i,1);
        y=cromosoma(i,2);
        w=cromosoma(i,3);
        evalCrom(i)=eval(z);
       % if evalCrom(i)> 0
            
        if (evalCrom(i)<zAuxMin)
            if evalCrom(i) > 0
              mejorCromosoma = i;
              zAuxMin = evalCrom(i);
            end
        end
        if evalCrom(i)>zAuxMax
            peorCromosoma = i;
            zAuxMax = evalCrom(i);
        end
      %  end
        i= i+1;
    end
    if(2<it)
        cromosoma(peorCromosoma,1)=GmejorCormosoma(1);
        cromosoma(peorCromosoma,2)=GmejorCormosoma(2);
        cromosoma(peorCromosoma,3)=GmejorCormosoma(3);
    end
    if zAuxMin<GZauxMin
        GmejorCormosoma(1)=cromosoma(mejorCromosoma,1);
        GmejorCormosoma(2)=cromosoma(mejorCromosoma,2);
        GmejorCormosoma(3)=cromosoma(mejorCromosoma,3);
        GZauxMin = zAuxMin;
        GZauxMin + "#It" + it
    end
    GZauxMin;
    %                                               -----------------SELECCIÔN------------------------------------------
    fit = [];
    fit = fitness(evalCrom,numCrom);        %se usa la funcion fitness para encontrar arreglo de Z' evaluados
    total = 0;      %sumatoria de todos los fitnessclc
    i = 1;          %iterador
    %suma de todos los fitness
    while i < numCrom +1
        total = total + fit(i);
        i= i+1;
    end
    total;
    p=[numCrom,2];% p(i) probabilidades de seleccion
    %    format long pp=[];
    pAnt=0;  %probailidad anterior --inicia en 0
    i = 1;          %iterador
    maxPP=0;
    posPP=0;
    while i < numCrom +1
        p(i,1)=pAnt;
        p(i,2) = pAnt+(fit(i) / total);
        pp(i)=(fit(i) / total);
        if pp(i)>maxPP
            maxPP=pp(i);
            posPP=i;
        end
        pAnt = p(i,2);
        i= i+1;
    end
    
    cromosoma = Ngen(cromosoma, p,numCrom);
    %                                               -----------------CRUCE---------------------------------------------
    cromosoma= PCruce(cromosoma, Pc,numCrom);
    %                                               -----------------MUTACIÔN------------------------------------------
    cromosoma=PMutacion(cromosoma, Pm, nMinValCrom, nMaxValCrom,numCrom);
    it =it+1;
    if(mod(it,100)==0)
        it
    end
end
"x1 =" + GmejorCormosoma(1)+ " x2=" + GmejorCormosoma(2) + " x3=" + GmejorCormosoma(3) +" phi=" + GZauxMin

% hacer que se varian n numeros
%probar que procd


%-------------------------------------------------------------------------
%                                   FUNCIONES
%-------------------------------------------------------------------------
%Fitness:
function f=fitness(evalCrom, numCrom)
    i = 1;
    f =[numCrom];
    while i < numCrom+1
        f(i)=1/(1+evalCrom(i));
        i= i+1;
    end
end
%-------------------------------------------------------------------------
%Nueva generacion
function n=Ngen(viejaGen, rangos,numCrom)
    i = 1;  %iterador
    n =[numCrom,3];
    while i < numCrom+1
        r = 0 + (1)*rand();
        %   r= rand(1,double)
        j=1;
        while r < rangos(j,1) || r > rangos(j,2)
            j = j+1;
        end
        n(i,1) = viejaGen(j,1);
        n(i,2) = viejaGen(j,2);
        n(i,3) = viejaGen(j,3);
        i= i+1;
    end
end
%-------------------------------------------------------------------------
%Proceso de cruce:
function copiaGens=PCruce(viejaGen, Pc,numCrom)
    i = 1;
    j = 1;
    c=[];
    while i < numCrom+1
        r = rand(1);
        if r <= Pc
            c(j,1)=viejaGen(i,1);
            c(j,2)=viejaGen(i,2);
            c(j,3)=viejaGen(i,3);
            c(j,4)=i;
            j = j+1;
        end
        i= i+1;
    end
    copiaGenss=[];
    i = 1;
    while i < numCrom+1
        copiaGenss(i,1)=viejaGen(i,1);
        copiaGenss(i,2)=viejaGen(i,2);
        copiaGenss(i,3)=viejaGen(i,3);
        i = i+1;
    end
    i = 1;
    k = 1;
    
    while i < length(c) && j > 4
        copiaGenss(c(i,4),1) = c(i,1);
        copiaGenss(c(i,4),2) = c(k+1,2);
        copiaGenss(c(i,4),3) = c(k+1,3);
        i = i+1;
        k = k+1;
        if k > j
            k = 1;
        end
    end
    copiaGens = copiaGenss;
end
%-------------------------------------------------------------------------
%Proceso de mutaciuon:
function cromosoma=PMutacion(cromosom, Pm, nMinValCrom, nMaxValCrom, numCrom)
    nMut = round((length(cromosom)*3)*Pm);
    i = 1;
    while i <= nMut
        %{
        r1= randi(3);
         r2 = randi(800);
         cromosom(r2,r1)=rand()*500;
        i= i+1;
        %}
        
        %revisar
        
        aux = randi([1 , length(cromosom)]);
        random=aleatorios();
        cromosom(aux,1)= random(1);
        cromosom(aux,2)= random(2);
        cromosom(aux,3)= random(3);
        i = i+1;
        
        %
        
    end
    cromosoma=cromosom;
end
%----------------------------------------------------------+---------------
%generacion aleatorios
function a=aleatorios()
    x=rand()*1000;
    y=rand()*(1000-x);
    z=1000-(x+y);
    if (0.14*x + 0.11*y + 0.1*z) <= 120
        a=aleatorios();
    else
        a=[x,y,z]  ;
    end 
end