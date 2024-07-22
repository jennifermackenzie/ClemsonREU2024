\\###############################################################
\\#-------------------------\r evenN.gp
\\ even period polynomial of level N
\\ 0, special zeros, circle 1/sqrt(N).
\\---------------------------------------------
\\ 
\\
\\ Date : 7/21/2024
\\
\\##############################################################

default(parisize, 90000000);
M = 20;
N = 17;

epzs = List(); \\ zeros on the circle |X|=1/sqrt(N) for level N

print(1/sqrt(N)); \\ output expected norm for all zeros

{for (k = 13, M, k++; 
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
			
		\\ construct the even period polynomial	
			epp = 0;
			for (n = 0, (k-2)/2,
				t = (-1)^n / factorial(2*n) * lfun(L,k-1-(2*n));
				epp = epp + t * (2*Pi*x)^(2*n);
			);
			pr = polroots(epp);
			epk = List([N,k,a]);
			epkz = List();
			print(length(pr));
			for(nz = 1, length(pr),
				print(abs(pr[nz])," & ", pr[nz]);
			);				
			listput(epk, vecsort(epkz));

			listput(epzs, epk);
		);
	);
)}











