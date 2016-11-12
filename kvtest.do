/*******************************************************************************
 KVREG Test
*******************************************************************************/

set more off
clear all

discard
do lnonparam.do 
* set seed 4904934
* set seed 548527343 
* local seed = c(seed)

    clear
    set obs 500
         
    gen vstar   = rnormal(0,1)
    gen zstar   = rnormal(0,1)
    gen ustar   = .33*vstar + zstar
    gen x1      = rnormal(0,1)
    gen x2      = rnormal(0,1)

    gen u =  exp(0.2*x1 + 0.6*x2)*ustar
    gen v =  exp(0.6*x1 + 0.2*x2)*vstar

    gen y2 = 1 + x1 + x2 + v
    gen y1 = 1 + x1 + x2 + y2 + u


timer clear
timer on 1
kvreg y1 y2 x1 x2 
timer off 1

timer list

ereturn list

