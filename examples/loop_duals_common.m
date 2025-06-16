(* ::Package:: *)

(* ::Text:: *)
(*Routine to compute the Omega matrix, common to all examples involving duals of loop integrals.*)


ComputeOmega[ii_]:=Module[
  {ann,fam,tmpids,seeds,eqs,ints,weight,diffop,dz,dmis,sol},
  
  (* computing annihilators in the first "ii" variables *)
  ann = CATAnnihilator[twist,zs[[1;;ii]],2,1];
  
  (* define regulated denominators - "regulated" is defined in each example *)
  CATDualRegulated["fam"] = regulated[[1;;ii]];
  
  (* compute template identities *)
  tmpids = CATDualInvMonIdsFromAnnihilators["fam",ann,zs[[1;;ii]]];
  
  (* seed up to rank 2 and 2 dots and generate system of equations *)
  seeds = Join@@(CATGenerateSeeds[#,{0,2},{0,2}]&/@sectors[[ii]]);
  eqs = Collect[Join@@(tmpids@@#&/@seeds),_CATInt,Expand];
  
  (* integral weight and list of integrals *)
  weight = {MemberQ[basis[[ii]],#]&,CATIntDualInvMonDefaultWeight};
  ints = CATIntCases[{eqs,CATInt["fam",#]&/@seeds},weight];
  
  (* derivatives of MIs from differential operators *)
  diffop = CATDiffOperator[twist,zs[[1;;ii]],zs[[ii+1;;ii+1]],1,1];
  dz = CATDualInvMonIdsFromDiffOperators["fam",diffop,zs[[1;;ii]]];
  dmis = (dz@@#)[[1]]&/@basis[[ii]][[;;,2]];
  
  (* solve for the integrals appearing in unreduced derivatives of MIs *)
  sol = FFSparseSolve[Thread[eqs==0],ints,
        "NeededVars"->CATIntCases[dmis,weight]
  ];
  dmis = dmis /. sol;
  If[!(CATIntCases[dmis,weight]===CATIntCases[basis[[ii]],weight]),
    Print["Failed to reduce to MIs"];
    Print["- expected MIs: ",CATIntCases[dmis,weight]];
    Print["- found MIs   : ",CATIntCases[basis[[ii]],weight]];
    Return[$Failed];
  ];
  Factor[Coefficient[#,basis[[ii]]]&/@(dmis)]
];
