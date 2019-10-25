***List Sets which will be the variables that go in the paramaters and bigger variables

sets
         m                       Index of energy storage use cases 1 to M(3)   /1, 2, 3/
         s                       Index of energy storage subsystems 1 to S
*t To be found with data
         t                       Index of time periods 0 to T
*To be found with optimization?
         t0(m)                   Set of time periods during which use case m is deployed
*Will have the same range as s
         s_prime                 Index of energy storage subsystem that charges another subsystem or who is being charged by another subsystem


*List parameters with an explanation of what they are so we don't get confused

Parameters
         a(t)                    Price of energy in time period t ($ per kWh)
         C_eff(s)                Charging efficiency of storage subsystem s (p.u.)
         D_eff(s)                Discharging efficiency of storage subsystem s (p.u.)
         delta                   Time step durantion (h)
         lambda(t,m)             Value of energy storage use case m in time period t ($ per kW)
         e_0(s)                  Initial state of charge of storage subsystem s (kWh)
         a_g(t)                  Actual renewable energy generation in time period t (kW)
         t_g(t)                  Target renewable energy generation in time period t (kW)
         beta(s)                 Cost of power capacity for storage subsystem technology s ($ per kWh)
         gamma(s)                Cost of energy storage capacity for storage subsystem technology s ($ per kW)
         rho(s)                  Minimum energy-to-power ratio of energy storage subsystem s (kWh per kW)
         epsylon(s)              Coefficient relating the maximum and minimum state of charge of storage subsystem s (p.u.)
         Fvv(s)                  Coefficient relating the maximum and minimum power of storage subsystem s (p.u.)
         theta_max               Investment budget ($)

*SOMETHING ABOUT OPTIONS!

*Include the main if we decide to do one: $Include "Nameofthefile";

*List all the variables
Variables
         p(s,t,m)                Net power of storage subsystem s during time period t allocated to storage use case m 9kW)
         p_c(s, s_prime,t)       Charging power of storage subsystem s from storage subsystem s_prime during time period t (kW) if s_prime = s then we are chargin from the grid
         p_d(s,s_prime,t)        Discharging power of storage subsystem s to storage subsystem s_prime during time period  (kW) if s_prime = s then discharging to the grid
         e(s,t)                  State of charge of storage subsystem s in period t (kWh)
         e_max(s)                Maximum state of charge of storage subsystem s(kWh)
         e_min(s)                Minimum state of charge of storage subsystem s (kWh)
         p_max(s)                Maximum power of state subsystem s (kWh)
         p_min(s)                Minimum power of storage subsystem s (kWh)
         u(s,t)                  Charging mode (1) or discharging mode (0) of the subsystem s in period t for idle mode doesn't matter since powers at 0

*Variable type
positive variables p, p_c, p_d, e, e_max,e_min, p_max;
Negative variable p_min;
binary variables u;


Equations
*Benefits for each case
         Peak_shaving_m1         Benefits of peak shaving
         Balancing_m2            Benefits of shaping and firming renewable generation (balancing)
         Price_arbitrage_m3      Benefits of selling stored energy at a higher price than when purchased
*Constraints Operational Model
         e_balance               Energy balance of the each energy subsystem
         p_balance               Power balance of each energy subsystem
*         e_lo                    Lower bound of the state of charge of each subsystem
*         e_up                    Higher bound of the state of charge of each subsystem
*         box_e_min               Condition of lower state of charge bound
*         box_e_max               Condition of higher state of charge bound
         e_max                   Equation to determine the maximum state of charge
         p_min                   Lower bound of the power of each subsystem
         p_max                   Higher bound of the power of each subsystem
         p_min_up                Condition of lower power bound
         p_max_lo                Condition of higher power bound
         Coeff_p                 Relationship between min_max power
         Coeff_e                 Relationship between min_max state of charge
*         p_d_lo                  Lower bound of discharging power
*         p_d_up                  Higher bound of discharging power
*         p_c_lo                  Lower bound of chargin power
*         p_c_up                  Higher bound of charging power
         epsylon_up              Epsylon's upper boundary
         Fvv_lo                  Power coefficient's lower boundary
         Cost                    Evaluation of costs
         lambda_lo_m1            Lower bound for Lambda for case 1
         lambda_lo_m2            Lower bound for Lambda for case 2
         lambda_lo_m3            Lower bound for Lambda for case 3
*         beta_con                Constraints on beta
*         lamda_con               Constraints on lamba
*         Cost_boundary           Investment budget restriction

*Equation to be maximised/minimised
          Obj_F                  Maximization function of the benefits of each storage subsystem
          Obj_C                  Summations of the costs of the storage system
          Obj_G                  Cost - Benefits function. Minimization of the costs vs benefits
 ;

*Constraints Operational Model
****Check if need to replace things for the boundaries of teh summations!
e_balance..              e(s,t) =e= e(s,t-1)+ delta*(C_eff(s)*sum(s_prime, p_c(s,s_prime,t))-sum(s_prime, p_d(s,s_prime,t));

p_balance..              sum(m, p(s,t,m)) =e= sum(s_prime,D_eff(s_prime)*p_d(s,s_prime,t)-p_c(s,s_prime,t));

e_min.lo(s) = 0;
*box_e_min..             e_min(s) =l= 0;

e_max.lo(s) = e_min(s);
*box_e_max..             e_max(s) > e_min(s);

e.lo(s,t) = e_min(s);
*e_lo..                  e_min(s) =l= e(s,t);

e.up(s,t) = e_max(s);
*e_up..                  e(s,t) =l= e_max(s);

e_max..                  e_max(s) =g= rho(s)*p_max(s);

p_min_up..               p_min(s) < 0;

p_max_lo..               p_max(s) > 0;

p_min.                   p_min(s) =l= sum(m,p(s,t,m));

p_max..                  sum(m,p(s,t,m)) =l= p_max(s);

Coeff_p..                p_min(s) =e= -Fvv(s)*p_max(s);

Coeff_e..                e_min(s) =e= epsylon*e_max(s);

p_d.lo(s,s_prime,t) = 0;
*p_d_lo..                0 =l= p_d(s,s_prime,t);

p_d.up(s,s_prime,t) = p_min(s)*(u(s,t)-1);
*p_d_up..                p_d(s,s_prime,t) =l= p_min(s)*(u(s,t)-1);

p_c.lo(s,s_prime,t)=0;
*p_c_lo..                0 =l= p_c(s,s_prime,t);

p_c.up(s,s_prime,t) = p_max(s)*u(s,t);
*p_c_up..                p_c(s,s_prime,t) =l= p_max(s)*u(s,t);

epsylon.lo(s) = 0;
*eps_cond_min..           0 =l= epsylon(s);

epsylon_up..             epsylon(s) < 1;

Fvv_lo..                0 < Fvv(s);

Fvv.up(s) = 1;
*Fvv_cond_max..           Fvv(s) =l= 1;

beta.lo(s) = 0;
*beta_con..               0 =l= beta(s);

gamma.lo(s) = 0;
*gamma_con..              0 =l= gamma(s);


*Boundaries of lambda for each cases of the three cases since they varie
lambda_lo_m1..            0 < lambda(t,1)

lambda_lo_m2..            0 =l= lambda(t,2)

lambda_lo_m3..            0 =e= lambda(t,3)

*Cost_boundary..

*Equation to be maximised/minimised
*Max
Ob

*Min
Obj_Gadd al

*Benefits for each case

Peak_shaving_m1

Balancing_m2

Price_arbitrage_m3

Costs



case1 storage_ben_fnt1.. =e= lambda(t,1)*p(s,t,1) + alpha(t)*p(s,t,1)
case2 storage_ben_fnt2.. =e= lambda(t,2)*abs(g(t)+p(s,t,1)-g(t))+alpha(t)*p(s,t,1)
case3 storage_ben_fnt2.. =e= alpha(t)*p(s,t,3)

res1.. =e= sum(t, sum(s, storage_ben_fnt1))
res2.. =e= sum(t, sum(s, storage_ben_fnt2))
res2.. =e= sum(t, sum(s, storage_ben_fnt3))

func.. =e= res1+res2+res3

*Il faut trouver combien de subsystem (soit le set s) que nous allons avoir! Bouffard avait dit qu'on devait optimiser ceci


*Do the equations
*1. name all of the equations with a description, so we don't get lost
*2. EquationName.. the equation itself


**REST TO BE DETERMINED
