SETS
        h horas /h1, h2, h3, h4/
        t ciudades /G, CATA, M, V, E, CAST/
        s escenario /s1, s2, s3, s4/
        alias(s,ss)
         ;

PARAMETERS
         prob(s) probabilidad escenario s /s1 0.36, s2 0.24, s3 0.24, s4 0.16/

         d(h) demanda termica en la hora h [MW] /h1 1000, h2 1400, h3 2300, h4 1500/
         al(t) termino lineal coste combustible de t /G 4, CATA 4, M 4, V 4, E 2, CAST 7/
         b(t) termino fijo cost de comb de t  /G 50, CATA 30, M 30, V 25, E 80, CAST 70/
         ca(t) coste de arranque de t  /G 10, CATA 20, M 10, V 15, E 20, CAST 15/
         cp(t) coste de parada de t  /G 5, CATA 10, M 5, V 10, E 15, CAST 10/
         pmax(t) potencia maxima de t  /G 400, CATA 500, M 700, V 400, E 1000, CAST 800/
         pmin(t) potencia minima de t /G 100, CATA 150, M 150, V 50, E 450, CAST 200/
         rs(t) rampa de subida de t  /G 200, CATA 300, M 500, V 300, E 600, CAST 400/
         rb(t) rampa de bajada de t   /G 300, CATA 300, M 200, V 100, E 600, CAST 400/

         r(h) porcentaje de r
         rsto(s,h) porcentaje de r estocastico
;
         r(h) = 0.2*d(h);

TABLE dems(s,h)  demanda termica en la hora h escenario s [MW]
                 h1      h2      h3      h4
         s1      1000    1400    2760    1800
         s2      0       0       0       1050
         s3      0       0       1610    1800
         s4      0       0       0       1050

;
rsto(s,h) = 0.2*dems(s,h);

TABLE arbol(s,h)
                 h1      h2      h3      h4
         s1      1       1       1       1
         s2      0       0       0       1
         s3      0       0       1       1
         s4      0       0       0       1
;

TABLE arbol1(s,h)
                 h1      h2      h3      h4
         s1      1       1       1       1
         s2      1       1       1       2
         s3      1       1       3       3
         s4      1       1       3       4
;

VARIABLES
PS(t,h,s) potencia t en h y s
AS(t,h,s) acoplamiento t en h y s
ARS(t,h,s) arranque t en h y s
*PR(t,h,s) parade t en h y s
Coste Coste total
;
positive variable PS
binary variable AS, ARS
*binary variable PR
;

EQUATIONS
FOBJECTIVO restriccion de la funcion objectivo
maxprod maximo de producion
rampasub rampa subida
rampabaja rampa bajada
reserva restircion de la reserva
demanda
minpot la potencia debe ser mas que la potencia minima
*rel1  relacion entre AR y A
*rel2   relacion entre PR y A
;

*FOBJECTIVO.. Coste =E= SUM((t,h), al(t)*P(t,h) + b(t)*A(t,h)+ ca(t)*AR(t,h) + cp(t)*PR(t,h));
FOBJECTIVO.. Coste =E= SUM((t,h,s,ss)$(ord(ss) = arbol1(s,h)), prob(s)*(al(t)*PS(t,h,ss) + b(t)*AS(t,h,ss) + ca(t)*ARS(t,h,ss) + cp(t)*(ARS(t,h,ss) + AS(t,h-1,ss) -AS(t,h,ss))));

demanda(h,s)$(ord(s) = arbol1(s,h)).. SUM(t,PS(t,h,s)) =E= dems(s,h);
reserva(h,s)$(ord(s) = arbol1(s,h)).. SUM(t, pmax(t)*AS(t,h,s)-PS(t,h,s)) =G= rsto(s,h);
maxprod(t,h,s)$(ord(s) = arbol1(s,h)).. PS(t,h,s) =L= pmax(t)*AS(t,h,s);
minpot(t,h,s)$(ord(s) = arbol1(s,h)).. PS(t,h,s) =G= pmin(t)*AS(t,h,s);
rampasub(t,h,s,ss)$(ord(s) = arbol1(s,h) and ord(ss) = arbol1(s,h-1)).. PS(t,h,s) - PS(t,h-1,ss) =L= rs(t);
rampabaja(t,h,s,ss)$(ord(s) = arbol1(s,h) and ord(ss) = arbol1(s,h-1)).. PS(t,h-1,ss) - PS(t,h,s) =L= rb(t);

*Example;
*rest(h,s,ss)$(ord(s) = arbol1(s,h) and ord(ss) = arbol1(s,h-1)) .. ;


*maxprod(t,h).. PS(t,h) =L= pmax(t)*AS(t,h);
*rampasub(t,h).. PS(t,h) - PS(t,h-1) =L= rs(t);
*rampabaja(t,h).. PS(t,h-1) - PS(t,h) =L= rb(t);
*minpot(t,h).. PS(t,h) =G= pmin(t)*AS(t,h);
*
*reserva(h).. sum(t, pmax(t)*A(t,h) - P(t,h)) =G= r(h);
*demanda(h).. sum(t, P(t,h)) =E= d(h);
*
*rel1(t,h).. AR(t,h) =G= A(t,h) - A(t,h-1);
*rel2(t,h).. PR(t,h) =G= A(t,h-1) - A(t,h);

option optcr = 0;
model pepito /all/;
solve pepito minimize Coste using MIP;
