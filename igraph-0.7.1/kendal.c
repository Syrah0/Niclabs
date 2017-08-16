// KENDALL
// Xi are the node from one metric and Yi are the node from the other metrica
// Xi and Yi are differents between them

// Calculate the number of concordant pairs and discordant pairs
// Concordant:
	if Xi < Xj and Yi < Yj => (Xi, Yi) and (Xj, Yj) concordants pairs
	if Xi > Xj and Yi > Yj => (Xi, Yi) and (Xj, Yj) concordants pairs
	if Xi == Xj or Yi == Yj => any classification
	in other case, they are discordant

	=> r = (concordant - discordant)/ [n(n-1)/2] 
	=> -1 <= r <= 1 ... if X and Y are independent => r ~ 0

COMBINATION BETWEEN ALL PAIRS -> X1,Y1 ~ X2,Y2 - X3,Y3 - ... - XN,YN

=> N = TOTAL OF PAIRS

double kendallCoef(int* X, int* Y){
	
}