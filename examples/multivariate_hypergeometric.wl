(* ::Package:: *)

(* ::Text:: *)
(*This example computes the reduction to master integrals (MIs) and the differential equations (DEs) satisfied by the hypergeometric functions (n+1)F(n)*)


<<CALICO`


(* ::Text:: *)
(*You can insert the number n identifying the hypergeometric function here.*)


n=4


(* ::Text:: *)
(*Define the number of integration variables and the twist eq.(5.12)*)


zs = z/@Range[n];


as = Join[Symbol["a"<>ToString[#]]&/@Range[n],{Symbol["a"<>ToString[n+1]]}];
bs = Symbol["b"<>ToString[#]]&/@Range[n];
polynomials = Flatten[{Table[{z[ii],(1-z[ii])},{ii,n}],(1-x Product[z[ii],{ii,n}])}];
exponents = Flatten[{Table[{(as[[ii]]-1),(-as[[ii]]+bs[[ii]]-1)},{ii,n}],(-as[[-1]])}];
twist = CATMultiBaikov[polynomials,exponents,zs];


(* ::Text:: *)
(*Computing  annihilators:*)


ann = CATAnnihilator[twist,zs,n+3,1];


Length/@ann[[1]]


(* ::Text:: *)
(*Higher order annihilators exist for n>1, but they are not needed for the reduction. Un-comment the following block to include them anyway.*)


(*
annextra = CATAnnihilator[twist,zs,n+2,2,"KnownSolutions"->ann];
ann = CATAnnihilatorMerge[ann,annextra];
*)


(* ::Text:: *)
(*Write  template identities from annihilators*)


tmpeqs = CATMonIdsFromAnnihilators["family",ann,zs];


(* ::Text:: *)
(*Generate a list of seeds*)


seedmaxdeg = 4;
seeds = Join@@(CATExponentList[#,n]&/@Range[0,seedmaxdeg]);


(* ::Text:: *)
(*Generate the linear system*)


eqs = Join@@(tmpeqs@@#&/@seeds);
eqs = DeleteCases[eqs,0];


(* ::Text:: *)
(*Choice of a weight *)


weight = {-Plus@@#[[2]],-Max[#[[2]]],#[[2]]}&;


(* ::Text:: *)
(*Define a list of integrals that we need to reduce to master integrals (i.e. having maximum degree `maxneededdeg`) .*)


ints = CATIntCases[eqs,weight];
maxneededdeg = 3;
neededints = Select[ints,
  Plus@@#[[2]] <= maxneededdeg&
];


(* ::Text:: *)
(*Finding the master integrals (MIs) as the independent unknowns of the system*)


mis = SortBy[FFSparseSolve[Thread[eqs==0],ints,"NeededVars"->neededints,"IndepVarsOnly"->True],weight]


(* ::Text:: *)
(*Finding the differential equations by deriving the MIs with respect to the external parameter x.*)


dop = CATDiffOperator[twist,z/@Range[n],{x},n+1,1];
dexop = CATMonIdsFromDiffOperators["family",dop,z/@Range[n]];
derivs = Collect[(dexop@@#)[[1]]&/@( mis[[;;,2]] ),_CATInt,Factor];


(* ::Text:: *)
(*Reducing to MIs the integrals appearing in the derivatives*)


reductions = FFSparseSolve[Thread[eqs==0],ints,"NeededVars"->CATIntCases[derivs,weight]];


(* ::Text:: *)
(*DE with respect to x*)


de = Coefficient[#,mis]&/@(derivs /. Dispatch[reductions])//Factor;
