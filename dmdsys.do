* Subject: demandsys webinar *
* Chris Cheng *

clear all

**# data and sample
webuse food_consumption
describe w* p* exp n_adults n_kids 

egen w_total = rowtotal(w*)
format %1.0f w_total
list w* in 1/5, abb(10)

list p* in 1/5, abb(10)

list exp in 1/5

list n* in 1/5

qui dtable w* p* exp n_*, sample(,place(seplabels)) sformat("(N=%s)" frequency)
collect style cell var[p_dairy expfd n_adults], border(top)
collect preview

**# aids
demandsys aids w_dairy w_proteins w_fruitveg w_flours w_misc, ///
	  prices(p_dairy p_proteins p_fruitveg p_flours p_misc)   ///
	  expenditure(expfd) nolog
	  
**# quaids
generate lnexp = ln(expfd)
twoway (qfitci w_dairy lnexp, fcolor(%40) acolor(%40)), ///
  ylabel(,format(%5.2f)) legend(off) /// 
  ytitle("Budget share") xtitle("Log of expenditure") ///                       
  title("Nonlinear relationship for dairy products") name(dairy, replace)	
twoway (qfitci w_proteins lnexp, fcolor(%40) acolor(%40)), ///
  ylabel(,format(%5.2f)) legend(off) /// 
  ytitle("Budget share") xtitle("Log of expenditure") ///                       
  title("Nonlinear relationship for fruits and vegetables") name(protein, replace)	
graph combine dairy protein

**## quaids estimation
demandsys quaids w_dairy w_proteins w_fruitveg w_flours w_misc, ///
        prices(p_dairy p_proteins p_fruitveg p_flours p_misc)   ///
	    expenditure(expfd) nolog

**## joint test	
test [lambda]1.Good [lambda]2.Good [lambda]3.Good [lambda]4.Good

**## demographics translating
demandsys quaids w_dairy w_proteins w_fruitveg w_flours w_misc, ///
        prices(p_dairy p_proteins p_fruitveg p_flours p_misc)   ///
	    expenditure(expfd) nolog demographics(n_kids n_adults)
				
**## demographics scaling
demandsys quaids w_dairy w_proteins w_fruitveg w_flours w_misc, ///
        prices(p_dairy p_proteins p_fruitveg p_flours p_misc)   ///
	    expenditure(expfd) nolog demographics(n_kids n_adults, scaling)		
		
**# postestimation
**## price sensitivities
estat elasticities, compensated
estat elasticities, uncompensated

**# make own and cross price tables
collect clear

* we need to compute the p-value as a matrix again since it is not returned
mata: 
V  = st_matrix("r(V)")
se = sqrt(diagonal(V))
st_matrix("SE",se)
end

matrix pvalue = r(b)
matrix se = SE'

forvalues i = 1/25 {
matrix pvalue[1,`i']= 2*(1-(normal(abs(r(b)[1,`i']/se[1,`i']))))
}

matlist pvalue

collect get C = r(b) P = pvalue

collect stars P 0.01 *** 0.05 ** 0.1 *, attach(C) shownote

collect title "Table 1. Own- and cross-price elasticities"
collect notes "See dmdsys.do for how to create this table using the -collect- system"
collect style title, font(,bold)

collect style header result, level(hide)

collect label levels coleq `"Good 1"' "Dairy"    ///
                           `"Good 2"' "Proteins" ///
						   `"Good 3"' "Fruits"   ///
						   `"Good 4"' "Flours"   ///
						   `"Good 5"' "Misc" 
collect label levels Good 1 "Dairy"    ///
                          2 "Proteins" ///
						  3 "Fruits"   ///
						  4 "Flours"   ///
						  5 "Misc"
						  
collect style cell coleq, halign(right)
collect style column, width(equal)
collect style cell Good, nformat(%5.3f) halign(center)

collect layout (coleq#result[C]) (Good)

**# make demand interrelationship table
matrix A = r(b)[1,1..5] \ r(b)[1,6..10] \ r(b)[1,11..15] \ r(b)[1,16..20] \ r(b)[1,21..25]

forvalues i = 1/5{
forvalues j = 1/5{
	if `j' == `i'{
		collect get rel = "--", tags(coleq[`"Good `i'"'] Good[`j'])
	}
	else{
	if  A[`i',`j'] > 0{
		collect get rel = "Subs", tags(coleq[`"Good `i'"'] Good[`j'])
	}
	else{
		collect get rel = "Comp", tags(coleq[`"Good `i'"'] Good[`j'])
	}
	}
}
}

collect title "Table 2. Demand interrelationship"
collect notes, clear
collect notes "Notes: Subs -- substitutes; Comp -- complements"
collect notes "See dmdsys.do for how to create this table using the -collect- system"

forvalues i = 1/5{
forvalues j = 1/5{
if `j' == `i'{
collect get emp = " ", tags(coleq[`"Good `i'"'] Good[`j'])
collect remap result[stars]=result[emp], fortags(coleq[`"Good `i'"'] Good[`j'])

}
}
}

collect composite define fin = rel emp, trim
collect style column, extraspace(2)

collect layout (coleq#result[fin]) (Good)

**## expenditure sensitivities
estat elasticities if n_kids == 2, expenditure atmeans


**# Program your own demand system - AIDS
clear all

capture program drop nlsurmyaids
program nlsurmyaids
version 18.0
syntax varlist(min=10 max=10) if, at(name)
tokenize `varlist'
args w1 w2 w3 w4 p1 p2 p3 p4 p5 exp
//defining alphas and adding-up
tempname a1 a2 a3 a4 a5
scalar `a1' = `at'[1,1]
scalar `a2' = `at'[1,2]
scalar `a3' = `at'[1,3]
scalar `a4' = `at'[1,4]
scalar `a5' = 1-`a1'-`a2'-`a3'-`a4'
//defining betas and adding-up
tempname b1 b2 b3 b4
scalar `b1' = `at'[1,5]
scalar `b2' = `at'[1,6]
scalar `b3' = `at'[1,7]
scalar `b4' = `at'[1,8]
//defining gammas, adding-up, and symmetry
tempname g11 g12 g13 g14 g15
tempname g21 g22 g23 g24 g25
tempname g31 g32 g33 g34 g35
tempname g41 g42 g43 g44 g45
tempname g51 g52 g53 g54 g55
scalar `g11' = `at'[1,9]
scalar `g12' = `at'[1,10]
scalar `g13' = `at'[1,11]
scalar `g14' = `at'[1,12]
scalar `g15' = -`g11'-`g12'-`g13'-`g14'
scalar `g21' = `g12'
scalar `g22' = `at'[1,13]
scalar `g23' = `at'[1,14]
scalar `g24' = `at'[1,15]
scalar `g25' = -`g21'-`g22'-`g23'-`g24'
scalar `g31' = `g13'
scalar `g32' = `g23'
scalar `g33' = `at'[1,16]
scalar `g34' = `at'[1,17]
scalar `g35' = -`g31'-`g32'-`g33'-`g34'
scalar `g41' = `g14'
scalar `g42' = `g24'
scalar `g43' = `g34'
scalar `g44' = `at'[1,18]
scalar `g45' = -`g41'-`g42'-`g43'-`g44'
scalar `g51' = `g15'
scalar `g52' = `g25'
scalar `g53' = `g35'
scalar `g54' = `g45'
scalar `g55' = -`g51'-`g52'-`g53'-`g54'
//price indexes/aggregators and budget share functions
quietly {
tempvar lnpindex
generate double `lnpindex' = 0 + `a1'*ln(`p1') + `a2'*ln(`p2') + ///
`a3'*ln(`p3') + `a4'*ln(`p4') + `a5'*ln(`p5')
forvalues i = 1/5 {
forvalues j = 1/5 {
replace `lnpindex' = `lnpindex' + ///
0.5*`g`i'`j''*ln(`p`i'')*ln(`p`j'')
}
}
replace `w1' = `a1' + `g11'*ln(`p1') + `g12'*ln(`p2') + ///
`g13'*ln(`p3') + `g14'*ln(`p4') + `g15'*ln(`p5') + ///
`b1'*(ln(`exp') - `lnpindex')
replace `w2' = `a2' + `g21'*ln(`p1') + `g22'*ln(`p2') + ///
`g23'*ln(`p3') + `g24'*ln(`p4') + `g25'*ln(`p5') + ///
`b2'*(ln(`exp') - `lnpindex')
replace `w3' = `a3' + `g31'*ln(`p1') + `g32'*ln(`p2') + ///
`g33'*ln(`p3') + `g34'*ln(`p4') + `g35'*ln(`p5') + ///
`b3'*(ln(`exp') - `lnpindex')
replace `w4' = `a4' + `g41'*ln(`p1') + `g42'*ln(`p2') + ///
`g43'*ln(`p3') + `g44'*ln(`p4') + `g45'*ln(`p5') + ///
`b4'*(ln(`exp') - `lnpindex')
}
end

webuse food_consumption, clear

nlsur myaids @ w_dairy w_proteins w_fruitveg w_flours                ///
               p_dairy p_proteins p_fruitveg p_flours p_misc expfd,  ///
             parameters(a1 a2 a3 a4 b1 b2 b3 b4 g11 g12 g13          ///
			 g14 g22 g23 g24 g33 g34 g44) neq(4) ifgnls

demandsys aids w_dairy w_proteins w_fruitveg w_flours w_misc, ///
	  prices(p_dairy p_proteins p_fruitveg p_flours p_misc)   ///
	  expenditure(expfd) piconstant(0)

**# Program your own demand system - QUAIDS
capture program drop nlsurmyquaids
program nlsurmyquaids
version 18.0
syntax varlist(min=10 max=10) if, at(name)
tokenize `varlist'
args w1 w2 w3 w4 p1 p2 p3 p4 p5 exp
//defining alphas and adding-up
tempname a1 a2 a3 a4 a5
scalar `a1' = `at'[1,1]
scalar `a2' = `at'[1,2]
scalar `a3' = `at'[1,3]
scalar `a4' = `at'[1,4]
scalar `a5' = 1-`a1'-`a2'-`a3'-`a4'
//defining betas and adding-up
tempname b1 b2 b3 b4 b5
scalar `b1' = `at'[1,5]
scalar `b2' = `at'[1,6]
scalar `b3' = `at'[1,7]
scalar `b4' = `at'[1,8]
scalar `b5' = -`b1'-`b2'-`b3'-`b4'
//defining gammas, adding-up, and symmetry
tempname g11 g12 g13 g14 g15
tempname g21 g22 g23 g24 g25
tempname g31 g32 g33 g34 g35
tempname g41 g42 g43 g44 g45
tempname g51 g52 g53 g54 g55
scalar `g11' = `at'[1,9]
scalar `g12' = `at'[1,10]
scalar `g13' = `at'[1,11]
scalar `g14' = `at'[1,12]
scalar `g15' = -`g11'-`g12'-`g13'-`g14'
scalar `g21' = `g12'
scalar `g22' = `at'[1,13]
scalar `g23' = `at'[1,14]
scalar `g24' = `at'[1,15]
scalar `g25' = -`g21'-`g22'-`g23'-`g24'
scalar `g31' = `g13'
scalar `g32' = `g23'
scalar `g33' = `at'[1,16]
scalar `g34' = `at'[1,17]
scalar `g35' = -`g31'-`g32'-`g33'-`g34'
scalar `g41' = `g14'
scalar `g42' = `g24'
scalar `g43' = `g34'
scalar `g44' = `at'[1,18]
scalar `g45' = -`g41'-`g42'-`g43'-`g44'
scalar `g51' = `g15'
scalar `g52' = `g25'
scalar `g53' = `g35'
scalar `g54' = `g45'
scalar `g55' = -`g51'-`g52'-`g53'-`g54'
//defning lambdas and adding-up
tempname l1 l2 l3 l4
scalar `l1' = `at'[1,19]
scalar `l2' = `at'[1,20]
scalar `l3' = `at'[1,21]
scalar `l4' = `at'[1,22]

//price indexes/aggregators and budget share functions
quietly {
tempvar lnpindex
generate double `lnpindex' = 0 + `a1'*ln(`p1') + `a2'*ln(`p2') + ///
`a3'*ln(`p3') + `a4'*ln(`p4') + `a5'*ln(`p5')
forvalues i = 1/5 {
forvalues j = 1/5 {
replace `lnpindex' = `lnpindex' + ///
0.5*`g`i'`j''*ln(`p`i'')*ln(`p`j'')
}
}
tempvar bp
generate double `bp' = ///
(`p1')^(`b1')*(`p2')^(`b2')*(`p3')^(`b3')*(`p4')^(`b4')*(`p5')^(`b5')

replace `w1' = `a1' + `g11'*ln(`p1') + `g12'*ln(`p2') + ///
`g13'*ln(`p3') + `g14'*ln(`p4') + `g15'*ln(`p5') + ///
`b1'*(ln(`exp') - `lnpindex') + `l1'/`bp'*(ln(`exp')-`lnpindex')^2
replace `w2' = `a2' + `g21'*ln(`p1') + `g22'*ln(`p2') + ///
`g23'*ln(`p3') + `g24'*ln(`p4') + `g25'*ln(`p5') + ///
`b2'*(ln(`exp') - `lnpindex') + `l2'/`bp'*(ln(`exp')-`lnpindex')^2
replace `w3' = `a3' + `g31'*ln(`p1') + `g32'*ln(`p2') + ///
`g33'*ln(`p3') + `g34'*ln(`p4') + `g35'*ln(`p5') + ///
`b3'*(ln(`exp') - `lnpindex') + `l3'/`bp'*(ln(`exp')-`lnpindex')^2
replace `w4' = `a4' + `g41'*ln(`p1') + `g42'*ln(`p2') + ///
`g43'*ln(`p3') + `g44'*ln(`p4') + `g45'*ln(`p5') + ///
`b4'*(ln(`exp') - `lnpindex') + `l4'/`bp'*(ln(`exp')-`lnpindex')^2
}
end

webuse food_consumption, clear

nlsur myquaids @ w_dairy w_proteins w_fruitveg w_flours                 ///
                 p_dairy p_proteins p_fruitveg p_flours p_misc expfd,   ///
               parameters(a1 a2 a3 a4 b1 b2 b3 b4 g11 g12 g13 g14 g22   ///
			   g23 g24 g33 g34 g44 l1 l2 l3 l4) neq(4) ifgnls

demandsys quaids w_dairy w_proteins w_fruitveg w_flours w_misc, ///
        prices(p_dairy p_proteins p_fruitveg p_flours p_misc)   ///
	    expenditure(expfd) piconstant(0)
