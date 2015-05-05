$TITLE Planificaci�n �ptima de la expansi�n de la generaci�n
SETS I generadores / gen-1 * gen-4 /
J periodos / per-1 * per-3 /
S escenarios de demanda / s-1 * s-3 /
PARAMETERS
         F(i) coste fijo de inversi�n [u.m.] / gen-1 10
                                               gen-2 7
                                               gen-3 16
                                               gen-4 6 /
         PROB(s) probabilidad cada escenario [p.u.] / s-1 0.2
                                                      s-2 0.5
                                                      s-3 0.3 /
         POTMIN potencia m�nima a instalar [MW] / 12 /
         PRSPTO l�mite presupuestario [u.m.] / 120 /
         DEM(j) demanda para un escenario [MW]
;

TABLE V(i,j) coste variable operaci�n [u.m.per MW]
                                                    per-1 per-2 per-3
                                             gen-1  40    24    4
                                             gen-2  45    27    4.5
                                             gen-3  32    19.2  3.2
                                             gen-4  55    33    5.5
;
TABLE DEMS(s,j) demanda estoc�stica [MW]
                                           per-1 per-2 per-3
                                       s-1 3     3     2
                                       s-2 5     3     2
                                       s-3 7     3     2
;
VARIABLES
X(i) potencia a instalar [MW]
Y(j,i) potencia en operaci�n [MW]
YS(s,j,i) potencia en operaci�n estoc�stica [MW]
COSTE coste total
;
POSITIVE VARIABLES X, Y, YS
;
EQUATIONS
COST coste total [pta]
COSTS coste total estoc�stico [pta]
PRESUP limitaci�n presupuestaria [pta]
INSMIN potencia m�nima instalada [MW]
BALPOT potencia en operaci�n menor que instalada [MW]
BALPOTS potencia en operaci�n menor que instalada estoc�stica [MW]
BALDEM balance de demanda [MW]
BALDEMS balance de demanda estoc�stico [MW] ;


COST .. COSTE =E= SUM(i, F(i) * X(i)) + SUM((j,i), V(i,j) * Y(j,i)) ;
COSTS .. COSTE =E= SUM(i, F(i) * X(i)) + SUM((s,j,i), PROB(s) * V(i,j) * YS(s,j,i)) ;
PRESUP .. SUM(i, F(i) * X(i)) =L= PRSPTO ;
INSMIN .. SUM(i, X(i)) =G= POTMIN ;
BALPOT(j,i).. Y(j,i) =L= X(i) ;
BALPOTS(s,j,i) .. YS(s,j,i) =L= X(i) ;
BALDEM(j).. SUM(i, Y(j,i)) =G= DEM(j) ;
BALDEMS(s,j) .. SUM(i, YS(s,j,i)) =G= DEMS(s,j) ;

MODEL DETERM / COST, INSMIN, PRESUP, BALPOT, BALDEM / ;
MODEL ESTOCA / COSTS, INSMIN, PRESUP, BALPOTS, BALDEMS /;
* cada escenario determinista por separado
LOOP (s, DEM(j) = DEMS(s,j) ;
         SOLVE DETERM MINIMIZING COSTE USING LP ;) ;
* escenario de demanda media:
*DEM(j) = SUM(s, PROB(s) * DEMS(s,j)) ;
*SOLVE DETERM MINIMIZING COSTE USING LP ;
* problema estoc�stico:
*SOLVE ESTOCA MINIMIZING COSTE USING LP
