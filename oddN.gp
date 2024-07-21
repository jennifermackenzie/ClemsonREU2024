\\###############################################################
\\#-------------------------\r oN.gp
\\ odd period polynomial of level N
\\ 0, special zeros, circle 1/sqrt(N).
\\---------------------------------------------
\\ 
\\
\\ Date : 6/21/2024
\\
\\##############################################################

default(parisize, 90000000);
M = 20;
N = 5;

opzs = List(); \\ zeros on the unit circle for level 1

print(1/sqrt(N)); \\ output expected norm for all zeros

{for (k = 9, M, k++; 
    chi = Mod(1,N);
    nk = mfdim([N,k,chi],0);
    printf("N is %d, k is %d, the dimension is %d\n", N, k, nk);
	if (nk < 1, next);
	
	\\---------L-functions of new forms
	mf = mfinit([N,k,chi],0);
	lf = mfeigenbasis(mf);
	La = 1;
	for(a = 1, length(lf),
		f = lf[a];
		LL = lfunmf(mf, f);
		\\ check whether LL is one function or an array of functions?
		\\ Warning: This method maybe wrong for large weight!!!!!!!!!!
		Lcnt = length(LL);   \\-----------number of L-functions.
		if(Lcnt == 6, Lcnt = 1);
		\\ If Lcnt is 6, we think LL contains only one function!!!!!!
		
		for(Ln = 1, Lcnt,
			L = 0;
			if(Lcnt > 1, 
				L = LL[Ln],
				L = LL
			);
		
			print("eigenform ", La);
			La++;
			
		\\ construct the odd period polynomial	
			opp = 0;
			for (n = 1, k-2,
				t = (-1)^((n-1)/2) / factorial(n)* lfun(L,k-n-1);
				opp = opp + t * (2*Pi*x)^n;
				n++;
			);
			pr = polroots(opp);
			opk = List([N,k,a]);
			opkz = List();
			print(length(pr));
			for(nz = 1, length(pr),
				print(abs(pr[nz])," & ", pr[nz]);
			);				
			listput(opk, vecsort(opkz));

			listput(opzs, opk);
		);
	);
)}











