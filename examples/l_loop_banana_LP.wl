(* ::Package:: *)

(* ::Text:: *)
(*This example computes differential equations (DEs) for the L-loop equal mass banana integral family, using the Lee-Pomeransky representation.*)


<<CALICO`


(* ::Text:: *)
(*You can set the number of loops here. The example has been tested up to L=6.*)


L=4


(* ::Text:: *)
(*Define loop denominators and the twist:*)


loopmoms = k/@Range[L];
dens = Join[ Table[{loopmoms[[ii]],m2},{ii,1,L}], {{Sum[loopmoms[[ii]],{ii,1,L}]-p,m2}} ]


{u,f,g,zs} = CATUFGPolys[
  dens,
  k/@Range[L],
  {p^2 -> s},
  z
];


twist = CATLP[g,d,zs];


(* ::Text:: *)
(*Subsectors factorize as a product of 1-loop tapole integrals. It is more efficient to derive and use 1-loop integral relations for those. Hence we define the denominators and the twist of a one-loop family.*)


dens1L = {{k,m2}};
{u1L,f1L,g1L,zs1L} = CATUFGPolys[
  dens1L,
  {k},
  {},
  z
]


twist1L = CATLP[g1L,d,zs1L];


(* ::Text:: *)
(*Computing annihilators and template identities:*)


ann = CATAnnihilator[twist,zs,L+1(*+1*),1]; //AbsoluteTiming


tmpeq = CATLPIdsFromAnnihilators["family",ann,zs,L,d];//AbsoluteTiming


(* ::Text:: *)
(*Computing annihilators and template identities for the 1-loop family:*)


ann1L = CATAnnihilator[twist1L,zs1L,4,1];


tmpeq1L = CATLPIdsFromAnnihilators["f1L",ann1L,zs1L,1,d];


(* ::Text:: *)
(*Convert one-loop identities into L-loop identities. We can build L+1 of these, i.e. one for each position of the denominator.*)
(*These identities are only valid when the integrand factorizes, i.e. when the numerator is trivial (this is handled when seeds are selected).*)


(*
   Convert 1L ids into L-loop id.s
*)
Module[{aleft,aright,a},
  tmpeqs1Lfactors = Join@@Table[
    aleft = a/@Range[1,ii-1];
    aright = a/@Range[ii+1,L+1];
    HeavisideTheta[a[ii]-1/2]tmpeq1L[a[ii]]/. CATInt->(CATInt["family",Join[aleft,#2,aright]]&),
  {ii,1,L+1}];
  tmpeqs1Lfactors = Evaluate[tmpeqs1Lfactors/.a->Slot]&;
]


(* ::Text:: *)
(*We have one top-sector and L+1 subsectors. However, we only write identities for one of the subsectors, since the other ones are related by a symmetry relation, namely a symmetry under any permutation of the exponents.*)


topsect = ConstantArray[1,L+1];
uniquesubsector = Join[ConstantArray[1,L],{0}];


(* ::Text:: *)
(*We use seed integrals up to L dots. Due to the symmetry of this family, we select only seeds that are independent modulo permutations the exponents, by making sure the exponents are sorted (from higher to lower).*)


filterseeds[seeds_List]:=Select[seeds, OrderedQ[#,-Order[##]&]&]


ndots = L; (* seed up to L dots at L loops *)
seeds = filterseeds[ CATGenerateSeeds[topsect,{0,ndots},{0,0}] ];
seeds1L = filterseeds[ CATGenerateSeeds[uniquesubsector,{0,ndots+L},{0,0}] ];


(* ::Text:: *)
(*The function simplifyeqs, simplifies the identities by removing zero integrals, which are all integrals with fewer than L denominators here, and it implements the symmetry relation, by sorting the exponents of all integrals in an expression.*)


intmap[expr_]:=expr/.CATInt->(CATInt[#1,Reverse[Sort[#2]]]&)


simplifyeqs[eqs_]:=Module[{zeroes},
  zeroes = Select[CATUnorderedIntCases[eqs],Length[Select[#[[2]], #>0&]]<L&];
  Collect[intmap[eqs] /. Dispatch[#->0&/@zeroes],_CATInt,Expand]
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


neededints = Select[ints, intDots[#[[2]]]<=ndots && intS[#[[2]]]==0&];


(* ::Text:: *)
(*Free parameters*)


params = {d,s,m2}


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
deriv = CATLPIdsFromDiffOperators["family",deop,zs,L,d];
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
(*Here res[[i]] is the differential equation matrix w.r.t. params[[1+i]] (i.e. s for i=1 and m2 for i=2). Each of them is the full DE for the MIs of the top sector, i.e. a {Length[mistop],Length[mis]}-dimensional matrix.*)


res = ArrayReshape[rec,{Length[params[[2;;]]],Length[mistop],Length[mis]}];
