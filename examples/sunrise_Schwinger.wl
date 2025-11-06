(* ::Package:: *)

(* ::Text:: *)
(*This example computes differential equations (DEs) for the 2-loop unequal mass banana integral family, using the Schwinger representation.*)


<<CALICO`


L=2;


(* ::Text:: *)
(*Define loop denominators and the twist:*)


loopmoms = k/@Range[L];
extmoms = p;
dens = Join[ Table[{loopmoms[[ii]],m2[ii]},{ii,1,L}], {{-Sum[loopmoms[[ii]],{ii,1,L}]+extmoms,m2[L+1]}} ];


{u,f,g,zs} = CATUFGPolys[
  dens,
  k/@Range[L],
  {p^2 -> s},
  z
];


twist = CATSchwinger[u,f,d,zs];


(* ::Text:: *)
(*Subsectors factorize as a product of 1-loop tapole integrals. It is more efficient to derive and use 1-loop integral relations for those. Hence we define the denominators and the twist of a one-loop family.*)


dens1L = {{k,m2}};
{u1L,f1L,g1L,zs1L} = CATUFGPolys[
  dens1L,
  {k},
  {},
  z
];


twist1L = CATSchwinger[u1L,f1L,d,zs1L];


(* ::Text:: *)
(*Computing annihilators and template identities:*)


ann = CATAnnihilator[twist,zs,L+2(*+1*),1]; //AbsoluteTiming


(* ::Text:: *)
(*Heuristically, we find that Schwinger parametrization needs more frequently annihilators of order 2 to perform a complete reduction.  We compute annihilators of order 2 that are independent modulo linear combinations and differentiation of the ones find at order 1 using the option "KnownSolutions".*)


ann2 = CATAnnihilator[twist,zs,2,2,"MinOrder"->1,"KnownSolutions"->ann]; //AbsoluteTiming


(* ::Text:: *)
(*We merge ann and ann2*)


ann = CATAnnihilatorMerge[ann,ann2];


tmpeq = CATSchwingerIdsFromAnnihilators["family",ann,zs];//AbsoluteTiming


(* ::Text:: *)
(*Computing annihilators and template identities for the 1-loop family:*)


ann1L = CATAnnihilator[twist1L,zs1L,4,1];


tmpeq1L = CATSchwingerIdsFromAnnihilators["f1L",ann1L,zs1L];


(* ::Text:: *)
(*Convert one-loop identities into 2-loop identities. We can build 3 of these, i.e. one for each position of the denominator.*)
(*These identities are only valid when the integrand factorizes, i.e. when the numerator is trivial (this is handled when seeds are selected).*)


(*
   Convert 1L ids into L-loop id.s
*)
Module[{aleft,aright,a},
  tmpeqs1Lfactors = Join@@Table[
    aleft = a/@Range[1,ii-1];
    aright = a/@Range[ii+1,L+1];
    HeavisideTheta[a[ii]-1/2]tmpeq1L[a[ii]]/. CATInt->(CATInt["family",Join[aleft,#2,aright]]&)/.m2->m2[ii],
  {ii,1,L+1}];
  tmpeqs1Lfactors = Evaluate[tmpeqs1Lfactors/.a->Slot]&;
]


(* ::Text:: *)
(*We have one top-sector and 3 subsectors.*)


topsect = ConstantArray[1,L+1];
lowersects = Table[Module[{sect=topsect}, sect[[ii]]=0; sect],{ii,L+1}];


(* ::Text:: *)
(*We use seed integrals up to 3 dots*)


ndots = 3;
seeds = Join[CATGenerateSeeds[topsect,{0,ndots},{0,0}],Join@@(CATGenerateSeeds[#,{0,ndots},{0,1}]&/@lowersects)];
seeds1L = Join@@(CATGenerateSeeds[#,{0,ndots+L},{0,(*1*)0},{MemberQ[#,0,1]&}]&/@lowersects);


(* ::Text:: *)
(*The function simplifyeqs, simplifies the identities by removing zero integrals, which are all integrals with fewer than 3 denominators.*)


simplifyeqs[eqs_]:=Module[{zeroes},
  zeroes = Select[CATUnorderedIntCases[eqs],Length[Select[#[[2]], #>0&]]<L&];
  Collect[eqs /. Dispatch[#->0&/@zeroes],_CATInt,Expand]
]


(* ::Text:: *)
(*Generate the linear system*)


eqs = DeleteCases[Join[
  Join@@simplifyeqs[tmpeqs1Lfactors@@#&/@seeds1L],
  Join@@simplifyeqs[tmpeq@@#&/@seeds]
  ]
,0];


(* ::Text:: *)
(*Find an overcomplete set of integrals appearing in the system.*)


ints = CATIntCases[{eqs,CATInt["family",#]&/@seeds,CATInt["family",#]&/@seeds1L},CATIntMellinDefaultWeight];


(* ::Text:: *)
(*Define a list of integrals that we need to reduce to master integrals (from their t,s,u). This is a list that is likely to be a superset of those needed for the DEs.*)


intT[a_]:=Length[Select[a,#>0&]];
intDots[a_]:=(Plus@@# - Length[#])&[Select[a,#>0&]];
intS[a_]:=-Plus@@Select[a,#<0&];


neededints = Select[ints, intDots[#[[2]]]<=ndots && intS[#[[2]]]<=0&];


(* ::Text:: *)
(*Free parameters*)


masses = m2/@Range[L+1];


params = Join[{d,s},masses]


(* ::Text:: *)
(*Here we use FiniteFlow to solve the system, following section 6.1 of https://arxiv.org/pdf/1905.08019. This example is simple enough that we can also easily find reduction tables with*)
(*    FFSparseSolve[Thread[eqs==0],ints,"NeededVars"->neededints]*)
(*but the overall strategy below is generally more efficient, especially in more complex cases.*)


FFNewGraph["g","in",params]


FFAlgSparseSolver["g","ibps",{"in"},params,Thread[eqs==0],ints,
"NeededVars"->neededints]


FFSolverOnlyHomogeneous["g","ibps"]


FFGraphOutput["g","ibps"]


ibplearn = FFSparseSolverLearn["g",ints];


(* ::Text:: *)
(*The independent unknowns of the system are the master integrals (MIs)*)


mis = "IndepVars"/.ibplearn


(* ::Text:: *)
(*Having found the MIs, we differentiate them with respect to s and m2. Here we focus on the derivatives of MIs belonging to the top sector, since those of lower sectors are single-scale integrals whose DE are uniquely fixed by dimensional analysis.*)


mistop = Select[mis,intT[#[[2]]]==L+1&];


deop = CATDiffOperator[twist,zs,params[[2;;]],4,1,"KnownAnnihilatorSolutions"->ann];
deriv = CATSchwingerIdsFromDiffOperators["family",deop,zs];
misderivatives = Collect[simplifyeqs[Transpose[deriv@@#[[2]]&/@mistop]],_CATInt,Together];


(* ::Text:: *)
(*List the integrals appearing in the derivatives and compute the DEs with FiniteFlow*)


newneeded = CATUnorderedIntCases[{misderivatives,mis}];


If[!SubsetQ[neededints,newneeded],
  Print["Missing needed integrals: ",Complement[newneeded,neededints]];
]


FFSolverResetNeededVars["g","ibps",ints,newneeded]


ibplearn = FFSparseSolverLearn["g",ints];


FFSparseSolverMarkAndSweepEqs["g","ibps"]


FFSparseSolverDeleteUnneededEqs["g","ibps"]


mis = "IndepVars"/.ibplearn;
ibpints = Join["DepVars"/.ibplearn,"IndepVars"/.ibplearn];


decoeffs = Map[Coefficient[#,ibpints]&,misderivatives,{2}];


FFAlgRatNumEval["g","misred",Join@@IdentityMatrix[Length[mis]]]


FFAlgChain["g","ibpsfull",{"ibps","misred"}]


FFAlgRatFunEval["g","deunreduced",{"in"},params,Flatten[decoeffs]]


FFAlgMatMul["g","de",{"deunreduced","ibpsfull"},
  Length[params[[2;;]]]*Length[mistop],Length[ibpints],Length[mis]
]


FFGraphOutput["g","de"]


rec=FFReconstructFunction["g",params];


(* ::Text:: *)
(*Here res[[i]] is the differential equation matrix w.r.t. params[[1+i]] (i.e. s for i=1, m2[1] for i=2, etc ). Each of them is the full DE for the MIs of the top sector, i.e. a {Length[mistop],Length[mis]}-dimensional matrix.*)


res = ArrayReshape[rec,{Length[params[[2;;]]],Length[mistop],Length[mis]}];
