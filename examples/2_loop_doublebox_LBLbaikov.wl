(* ::Package:: *)

<<CALICO`


loopmoms = {k1,k2};
extmoms = {p1,p2,p3}; (* p4 = -p1-p2-p3 *)
dens = {
k1,
k1+p1,
k1+p1+p2,
k1+k2,
k2,
k2-p1-p2-p3,
k2-p1-p2,
(*isps*)
k1+p1+p2+p3(*,
k2-p1*)  (* second ISP not needed *)
};
repl = {p3^2 -> 0, p2^2 -> 0, p1^2 -> 0, p2*p3 -> s23/2, 
        p1*p2 -> s12/2, p1*p3 -> -1/2*s12 - s23/2};


(* ::Text:: *)
(*Polynomials and exponents of the loop-by-loop (LBL) Baikov parametrization. These do not include, yet, the ones for the external prefactor which only depends on kinematic variables since they don't contribute to reduction identities (but they will be included later when deriving differential equations).*)


polys = {s12^2 s23^2-2 s12 s23^2 z[1]+s23^2 z[1]^2-2 s12^2 s23 z[2]+2 s12 s23 z[1] z[2]+s12^2 z[2]^2-2 s12 s23^2 z[3]-4 s12 s23 z[1] z[3]-2 s23^2 z[1] z[3]+2 s12 s23 z[2] z[3]+s23^2 z[3]^2-2 s12^2 s23 z[8]+2 s12 s23 z[1] z[8]-2 s12^2 z[2] z[8]-4 s12 s23 z[2] z[8]+2 s12 s23 z[3] z[8]+s12^2 z[8]^2,-s12 (z[1] z[3]+s12 z[8]-z[1] z[8]-z[3] z[8]+z[8]^2),s12^2 z[4]^2-2 s12 z[3] z[4] z[5]+z[3]^2 z[5]^2-4 s12 z[1] z[3] z[6]-2 s12^2 z[4] z[6]+2 s12 z[1] z[4] z[6]+2 s12 z[3] z[4] z[6]+2 s12 z[3] z[5] z[6]+2 z[1] z[3] z[5] z[6]-2 z[3]^2 z[5] z[6]+s12^2 z[6]^2-2 s12 z[1] z[6]^2+z[1]^2 z[6]^2-2 s12 z[3] z[6]^2-2 z[1] z[3] z[6]^2+z[3]^2 z[6]^2-2 s12 z[1] z[4] z[7]-2 z[1] z[3] z[5] z[7]+2 s12 z[1] z[6] z[7]-2 z[1]^2 z[6] z[7]+2 z[1] z[3] z[6] z[7]+z[1]^2 z[7]^2-2 s12^2 z[4] z[8]+2 s12 z[3] z[5] z[8]+2 s12 z[4] z[5] z[8]-2 z[3] z[5]^2 z[8]-2 s12^2 z[6] z[8]+2 s12 z[1] z[6] z[8]+2 s12 z[3] z[6] z[8]-4 s12 z[4] z[6] z[8]+2 s12 z[5] z[6] z[8]-2 z[1] z[5] z[6] z[8]+2 z[3] z[5] z[6] z[8]+2 s12 z[1] z[7] z[8]+2 s12 z[4] z[7] z[8]-4 s12 z[5] z[7] z[8]+2 z[1] z[5] z[7] z[8]+2 z[3] z[5] z[7] z[8]+2 s12 z[6] z[7] z[8]+2 z[1] z[6] z[7] z[8]-2 z[3] z[6] z[7] z[8]-2 z[1] z[7]^2 z[8]+s12^2 z[8]^2-2 s12 z[5] z[8]^2+z[5]^2 z[8]^2-2 s12 z[7] z[8]^2-2 z[5] z[7] z[8]^2+z[7]^2 z[8]^2};
exponents = {-5/2 + d/2, 2 - d/2, -5/2 + d/2};
zs = {z[1],z[2],z[3],z[4],z[5],z[6],z[7],z[8]};


(* ::Text:: *)
(*The input above can be generated using the "BaikovPackage.m" by Hjalte Frellesvig: https://github.com/HjalteFrellesvig/BaikovPackage*)
(*You may do that by un-commenting the following block, if you have the package installed.*)


(*
Get["BaikovPackage.m"];
Internal = loopmoms;
External = extmoms;
Propagators = dens[[1;;7]]^2;
PropagatorsExtra = {dens[[8]]^2};
Replacements = repl;
baikovcomp = BaikovLBL[] /. x->z;

(* note: first polynomial yields a common loop-independent prefactor, hence we remove it *)
polys = baikovcomp[[1]][[2;;]];
exponents = baikovcomp[[2]][[2;;]];
zs = Union[Cases[baikovcomp[[1]],_x,Infinity]] /. x->z;
*)


(* ::Text:: *)
(*Define twist and compute annihilators*)


twist = CATMultiBaikov[polys,exponents,zs];


ann = CATAnnihilator[twist,zs,3,1];


Length/@ann[[1]]


(* ::Text:: *)
(*There are second order annihilators of degree 1, but they are not needed for the reduction below. To include them in the set of identities, un-comment the following block.*)


(*
annext = CATAnnihilator[twist,zs,1,2,"MinOrder"->2,"KnownSolutions"->ann];
ann = CATAnnihilatorMerge[ann,annext];
*)


(* ::Text:: *)
(*Compute template identities from annihilators.*)


tmpids = CATInvMonIdsFromAnnihilators["fam",ann,zs];


(* ::Text:: *)
(*Zero and non-zero sectors*)


nonzerosects = {{0,0,1,1,1,0,0,0},{0,0,1,1,1,0,1,0},{0,0,1,1,1,1,0,0},{0,0,1,1,1,1,1,0},{0,1,0,1,0,1,0,0},{0,1,0,1,0,1,1,0},{0,1,0,1,1,0,1,0},{0,1,0,1,1,1,0,0},{0,1,0,1,1,1,1,0},{0,1,1,1,0,1,0,0},{0,1,1,1,0,1,1,0},{0,1,1,1,1,0,0,0},{0,1,1,1,1,0,1,0},{0,1,1,1,1,1,0,0},{0,1,1,1,1,1,1,0},{1,0,0,1,0,0,1,0},{1,0,0,1,0,1,1,0},{1,0,0,1,1,0,1,0},{1,0,0,1,1,1,1,0},{1,0,1,0,1,0,1,0},{1,0,1,0,1,1,1,0},{1,0,1,1,0,0,1,0},{1,0,1,1,0,1,0,0},{1,0,1,1,0,1,1,0},{1,0,1,1,1,0,0,0},{1,0,1,1,1,0,1,0},{1,0,1,1,1,1,0,0},{1,0,1,1,1,1,1,0},{1,1,0,1,0,0,1,0},{1,1,0,1,0,1,0,0},{1,1,0,1,0,1,1,0},{1,1,0,1,1,0,1,0},{1,1,0,1,1,1,0,0},{1,1,0,1,1,1,1,0},{1,1,1,0,1,0,1,0},{1,1,1,0,1,1,1,0},{1,1,1,1,0,0,1,0},{1,1,1,1,0,1,0,0},{1,1,1,1,0,1,1,0},{1,1,1,1,1,0,0,0},{1,1,1,1,1,0,1,0},{1,1,1,1,1,1,0,0},{1,1,1,1,1,1,1,0}};
zerosects = {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,0,0,0,1,0,0,0},{0,0,0,1,0,0,0,0},{0,0,1,0,0,0,0,0},{0,1,0,0,0,0,0,0},{1,0,0,0,0,0,0,0},{0,0,0,0,0,1,1,0},{0,0,0,0,1,0,1,0},{0,0,0,0,1,1,0,0},{0,0,0,1,0,0,1,0},{0,0,0,1,0,1,0,0},{0,0,0,1,1,0,0,0},{0,0,1,0,0,0,1,0},{0,0,1,0,0,1,0,0},{0,0,1,0,1,0,0,0},{0,0,1,1,0,0,0,0},{0,1,0,0,0,0,1,0},{0,1,0,0,0,1,0,0},{0,1,0,0,1,0,0,0},{0,1,0,1,0,0,0,0},{0,1,1,0,0,0,0,0},{1,0,0,0,0,0,1,0},{1,0,0,0,0,1,0,0},{1,0,0,0,1,0,0,0},{1,0,0,1,0,0,0,0},{1,0,1,0,0,0,0,0},{1,1,0,0,0,0,0,0},{0,0,0,0,1,1,1,0},{0,0,0,1,0,1,1,0},{0,0,0,1,1,0,1,0},{0,0,0,1,1,1,0,0},{0,0,1,0,0,1,1,0},{0,0,1,0,1,0,1,0},{0,0,1,0,1,1,0,0},{0,0,1,1,0,0,1,0},{0,0,1,1,0,1,0,0},{0,1,0,0,0,1,1,0},{0,1,0,0,1,0,1,0},{0,1,0,0,1,1,0,0},{0,1,0,1,0,0,1,0},{0,1,0,1,1,0,0,0},{0,1,1,0,0,0,1,0},{0,1,1,0,0,1,0,0},{0,1,1,0,1,0,0,0},{0,1,1,1,0,0,0,0},{1,0,0,0,0,1,1,0},{1,0,0,0,1,0,1,0},{1,0,0,0,1,1,0,0},{1,0,0,1,0,1,0,0},{1,0,0,1,1,0,0,0},{1,0,1,0,0,0,1,0},{1,0,1,0,0,1,0,0},{1,0,1,0,1,0,0,0},{1,0,1,1,0,0,0,0},{1,1,0,0,0,0,1,0},{1,1,0,0,0,1,0,0},{1,1,0,0,1,0,0,0},{1,1,0,1,0,0,0,0},{1,1,1,0,0,0,0,0},{0,0,0,1,1,1,1,0},{0,0,1,0,1,1,1,0},{0,0,1,1,0,1,1,0},{0,1,0,0,1,1,1,0},{0,1,1,0,0,1,1,0},{0,1,1,0,1,0,1,0},{0,1,1,0,1,1,0,0},{0,1,1,1,0,0,1,0},{1,0,0,0,1,1,1,0},{1,0,0,1,1,1,0,0},{1,0,1,0,0,1,1,0},{1,0,1,0,1,1,0,0},{1,1,0,0,0,1,1,0},{1,1,0,0,1,0,1,0},{1,1,0,0,1,1,0,0},{1,1,0,1,1,0,0,0},{1,1,1,0,0,0,1,0},{1,1,1,0,0,1,0,0},{1,1,1,0,1,0,0,0},{1,1,1,1,0,0,0,0},{0,1,1,0,1,1,1,0},{1,1,0,0,1,1,1,0},{1,1,1,0,0,1,1,0},{1,1,1,0,1,1,0,0}};


(* ::Text:: *)
(*Seed the reduction.*)


maxrank=2;
maxdots=2;
seeds = Join@@(CATGenerateSeeds[#,{0,maxdots},{0,maxrank}]&/@nonzerosects);


(* ::Text:: *)
(*Utilities to simplify the equations*)


sectorof = (HeavisideTheta[#-1/2]&/@#[[2]])&


simplifyeqs[expr_]:=Module[{zeroes},
  zeroes = Select[CATUnorderedIntCases[expr],MemberQ[zerosects,sectorof[#]]&];
  Collect[expr /. Dispatch[#->0&/@zeroes], _CATInt]
]


(* ::Text:: *)
(*Generate the system*)


equations = simplifyeqs[ Flatten[tmpids@@#&/@seeds] ];


(* ::Text:: *)
(*This family has non-trivial symmetries, which have to be found via different methods. The following have been found using FFIntRed (to be released) but in principle any software capable of finding sector mappings and sector symmetries would be able to generate these equations.*)


mappings = {-CATInt["fam",{0,0,1,1,1,0,0,0}]+CATInt["fam",{1,0,0,1,0,0,1,0}],3 CATInt["fam",{0,1,0,1,0,1,0,-1}]-s23 CATInt["fam",{0,1,0,1,0,1,0,0}],CATInt["fam",{0,1,0,1,0,1,-1,0}]+CATInt["fam",{0,1,0,1,0,1,0,-1}],CATInt["fam",{0,1,0,1,-1,1,0,0}]+CATInt["fam",{0,1,0,1,0,1,0,-1}],CATInt["fam",{0,1,-1,1,0,1,0,0}]+CATInt["fam",{0,1,0,1,0,1,0,-1}],CATInt["fam",{0,0,1,1,1,0,-1,0}]+CATInt["fam",{0,0,1,1,1,0,0,-1}],CATInt["fam",{0,0,1,1,1,0,-1,0}]-2 CATInt["fam",{0,0,1,1,1,0,0,-1}]-s12 CATInt["fam",{0,0,1,1,1,0,0,0}],CATInt["fam",{0,0,1,1,1,-1,0,0}]-CATInt["fam",{0,0,1,1,1,0,0,-1}],CATInt["fam",{0,-1,1,1,1,0,0,0}]-CATInt["fam",{0,0,1,1,1,0,0,-1}],CATInt["fam",{-1,1,0,1,0,1,0,0}]+CATInt["fam",{0,1,0,1,0,1,0,-1}],CATInt["fam",{-1,0,1,1,1,0,0,0}]+CATInt["fam",{0,0,1,1,1,0,0,-1}],-CATInt["fam",{-1,0,1,1,1,0,0,0}]+CATInt["fam",{0,0,1,1,1,0,0,-1}]+s12 CATInt["fam",{0,0,1,1,1,0,0,0}]+CATInt["fam",{1,0,0,1,0,0,1,-1}],CATInt["fam",{0,0,1,1,1,-1,0,0}]-CATInt["fam",{0,0,1,1,1,0,-1,0}]+s12 CATInt["fam",{0,0,1,1,1,0,0,0}]+CATInt["fam",{1,0,0,1,0,-1,1,0}],-CATInt["fam",{0,0,1,1,1,0,-1,0}]+CATInt["fam",{1,0,0,1,-1,0,1,0}],-CATInt["fam",{-1,0,1,1,1,0,0,0}]+CATInt["fam",{1,0,-1,1,0,0,1,0}],-CATInt["fam",{-1,0,1,1,1,0,0,0}]+CATInt["fam",{0,-1,1,1,1,0,0,0}]+s12 CATInt["fam",{0,0,1,1,1,0,0,0}]+CATInt["fam",{1,-1,0,1,0,0,1,0}],-CATInt["fam",{0,1,0,1,0,1,1,0}]+CATInt["fam",{1,1,0,1,0,1,0,0}],-CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{1,1,0,1,0,0,1,0}],-CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,1,0,0,0}],-CATInt["fam",{0,1,0,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,0,1,0,0}],-CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,0,0,1,0}],-CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,0,1,1,0,1,0}],-CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{1,0,0,1,0,1,1,0}],-CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{0,1,1,1,1,0,0,0}],-CATInt["fam",{0,1,0,1,0,1,1,0}]+CATInt["fam",{0,1,1,1,0,1,0,0}],-CATInt["fam",{0,1,0,1,0,1,1,0}]+CATInt["fam",{0,1,0,1,1,1,0,0}],2 CATInt["fam",{1,0,1,0,1,0,1,-1}]+s12 CATInt["fam",{1,0,1,0,1,0,1,0}],2 CATInt["fam",{1,0,1,0,1,-1,1,0}]+s12 CATInt["fam",{1,0,1,0,1,0,1,0}],2 CATInt["fam",{1,0,1,-1,1,0,1,0}]+s12 CATInt["fam",{1,0,1,0,1,0,1,0}],2 CATInt["fam",{1,-1,1,0,1,0,1,0}]+s12 CATInt["fam",{1,0,1,0,1,0,1,0}],2 CATInt["fam",{0,1,-1,1,1,0,1,0}]+CATInt["fam",{0,1,0,1,1,-1,1,0}]-2 CATInt["fam",{0,1,0,1,1,0,1,-1}]+s23 CATInt["fam",{0,1,0,1,1,0,1,0}],-CATInt["fam",{0,0,1,1,1,0,0,0}]+CATInt["fam",{0,0,1,1,1,1,-1,0}]+2 CATInt["fam",{0,0,1,1,1,1,0,-1}],2 CATInt["fam",{0,0,1,1,1,-1,1,0}]-CATInt["fam",{0,0,1,1,1,0,0,0}]+s12 CATInt["fam",{0,0,1,1,1,0,1,0}],CATInt["fam",{0,0,1,1,1,-1,1,0}]-CATInt["fam",{0,0,1,1,1,0,0,0}]-2 CATInt["fam",{0,0,1,1,1,0,1,-1}],CATInt["fam",{0,-1,1,1,1,0,1,0}]-CATInt["fam",{0,0,1,1,1,0,1,-1}],CATInt["fam",{-1,1,0,1,1,0,1,0}]-CATInt["fam",{0,1,-1,1,1,0,1,0}],2 CATInt["fam",{-1,0,1,1,1,1,0,0}]+CATInt["fam",{0,0,1,1,1,1,-1,0}]-s12 CATInt["fam",{0,0,1,1,1,1,0,0}],CATInt["fam",{-1,0,1,1,1,0,1,0}]+CATInt["fam",{0,0,1,1,1,-1,1,0}],2 CATInt["fam",{0,1,0,1,0,1,1,-1}]-s23 CATInt["fam",{0,1,0,1,0,1,1,0}]+CATInt["fam",{1,1,0,1,0,1,0,-1}],2 CATInt["fam",{0,1,-1,1,0,1,1,0}]-CATInt["fam",{0,1,0,1,0,1,0,0}]+CATInt["fam",{1,1,0,1,0,1,0,-1}],2 CATInt["fam",{-1,1,0,1,0,1,1,0}]-CATInt["fam",{0,1,0,1,-1,1,1,0}]+CATInt["fam",{1,1,0,1,0,1,0,-1}],CATInt["fam",{-1,1,0,1,0,1,1,0}]-CATInt["fam",{0,1,0,1,-1,1,1,0}]+CATInt["fam",{1,1,0,1,0,1,-1,0}]+CATInt["fam",{1,1,0,1,0,1,0,-1}],-2 CATInt["fam",{0,-1,1,1,1,1,0,0}]-CATInt["fam",{0,0,1,1,1,1,-1,0}]+CATInt["fam",{1,1,0,1,0,0,1,-1}],CATInt["fam",{0,-1,1,1,1,1,0,0}]+CATInt["fam",{0,0,1,1,1,1,-1,0}]+CATInt["fam",{1,1,0,1,0,-1,1,0}]-CATInt["fam",{1,1,0,1,0,0,1,-1}],CATInt["fam",{0,1,-1,1,0,1,1,0}]-CATInt["fam",{0,1,0,1,0,1,0,0}]+CATInt["fam",{1,1,0,1,-1,1,0,0}]+CATInt["fam",{1,1,0,1,0,1,0,-1}],CATInt["fam",{-1,0,1,1,1,1,0,0}]+CATInt["fam",{0,0,1,1,1,1,-1,0}]-s12 CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{1,1,0,1,-1,0,1,0}],-CATInt["fam",{0,1,0,1,-1,1,1,0}]+CATInt["fam",{1,1,-1,1,0,1,0,0}],-CATInt["fam",{0,0,1,1,1,1,-1,0}]+CATInt["fam",{1,1,-1,1,0,0,1,0}],CATInt["fam",{0,0,1,1,1,-1,1,0}]-CATInt["fam",{0,0,1,1,1,0,0,0}]+s12 CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,1,0,0,-1}],CATInt["fam",{-1,0,1,1,1,0,1,0}]+CATInt["fam",{0,0,1,1,1,0,0,0}]-s12 CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,1,0,-1,0}],CATInt["fam",{-1,0,1,1,1,0,1,0}]+CATInt["fam",{0,0,1,1,1,-1,1,0}]-CATInt["fam",{0,0,1,1,1,0,1,-1}]+CATInt["fam",{1,0,1,1,1,-1,0,0}],2 CATInt["fam",{0,1,-1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,0,1,0,-1}],CATInt["fam",{-1,1,0,1,1,0,1,0}]-2 CATInt["fam",{0,1,-1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,0,1,-1,0}],-CATInt["fam",{0,0,1,1,1,-1,1,0}]+CATInt["fam",{1,0,1,1,0,0,1,-1}],-CATInt["fam",{0,0,1,1,1,-1,1,0}]+CATInt["fam",{0,0,1,1,1,0,0,0}]+CATInt["fam",{0,0,1,1,1,0,1,-1}]+CATInt["fam",{1,0,1,1,0,-1,1,0}],-CATInt["fam",{0,1,-1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,-1,1,0,0}],CATInt["fam",{-1,0,1,1,1,0,1,0}]+CATInt["fam",{0,0,1,1,1,0,0,0}]-s12 CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,-1,0,1,0}],-CATInt["fam",{-1,0,1,1,1,0,1,0}]+CATInt["fam",{0,0,1,1,1,0,1,-1}]+s12 CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,0,1,1,0,1,-1}],CATInt["fam",{0,0,1,1,1,-1,1,0}]-CATInt["fam",{0,0,1,1,1,0,0,0}]+s12 CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,0,1,1,-1,1,0}],-CATInt["fam",{0,0,1,1,1,1,0,-1}]+CATInt["fam",{1,0,0,1,0,1,1,-1}],-CATInt["fam",{0,0,1,1,1,1,-1,0}]+CATInt["fam",{1,0,0,1,-1,1,1,0}],-CATInt["fam",{-1,0,1,1,1,0,1,0}]+CATInt["fam",{1,0,-1,1,1,0,1,0}],-CATInt["fam",{-1,0,1,1,1,1,0,0}]+CATInt["fam",{1,0,-1,1,0,1,1,0}],CATInt["fam",{0,0,1,1,1,-1,1,0}]-CATInt["fam",{0,0,1,1,1,0,0,0}]+s12 CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,-1,1,1,1,0,0,0}],-CATInt["fam",{0,1,0,1,1,-1,1,0}]+CATInt["fam",{1,-1,1,1,0,1,0,0}],-CATInt["fam",{0,0,1,1,1,-1,1,0}]+CATInt["fam",{1,-1,1,1,0,0,1,0}],-CATInt["fam",{-1,0,1,1,1,0,1,0}]+CATInt["fam",{0,-1,1,1,1,0,1,0}]+s12 CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{1,-1,0,1,1,0,1,0}],-CATInt["fam",{0,-1,1,1,1,1,0,0}]+CATInt["fam",{1,-1,0,1,0,1,1,0}],CATInt["fam",{0,1,1,1,1,0,0,-1}]-CATInt["fam",{1,1,0,1,0,0,1,-1}],CATInt["fam",{-1,0,1,1,1,1,0,0}]+CATInt["fam",{0,0,1,1,1,1,-1,0}]-s12 CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{0,1,1,1,1,0,-1,0}],CATInt["fam",{0,-1,1,1,1,1,0,0}]+CATInt["fam",{0,0,1,1,1,1,-1,0}]+CATInt["fam",{0,1,1,1,1,-1,0,0}]-CATInt["fam",{1,1,0,1,0,0,1,-1}],CATInt["fam",{0,1,1,1,0,1,0,-1}]-CATInt["fam",{1,1,0,1,0,1,0,-1}],CATInt["fam",{0,1,-1,1,0,1,1,0}]-CATInt["fam",{0,1,0,1,0,1,0,0}]+CATInt["fam",{0,1,1,1,0,1,-1,0}]+CATInt["fam",{1,1,0,1,0,1,0,-1}],CATInt["fam",{-1,1,0,1,0,1,1,0}]-CATInt["fam",{0,1,0,1,-1,1,1,0}]+CATInt["fam",{0,1,1,1,-1,1,0,0}]+CATInt["fam",{1,1,0,1,0,1,0,-1}],-CATInt["fam",{0,1,0,1,0,1,1,-1}]+CATInt["fam",{0,1,0,1,1,1,0,-1}],-CATInt["fam",{0,1,0,1,-1,1,1,0}]+CATInt["fam",{0,1,0,1,1,1,-1,0}],-CATInt["fam",{-1,1,0,1,0,1,1,0}]+CATInt["fam",{0,1,-1,1,1,1,0,0}],CATInt["fam",{-1,1,1,1,1,0,0,0}]-CATInt["fam",{0,0,1,1,1,1,-1,0}],CATInt["fam",{-1,1,1,1,0,1,0,0}]-CATInt["fam",{0,1,0,1,-1,1,1,0}],CATInt["fam",{-1,1,0,1,1,1,0,0}]-CATInt["fam",{0,1,-1,1,0,1,1,0}],-CATInt["fam",{0,0,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,1,0,0,0}],-CATInt["fam",{0,1,0,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,0,1,0,0}],-CATInt["fam",{0,0,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,0,0,1,0}],-CATInt["fam",{1,0,1,0,1,1,1,0}]+CATInt["fam",{1,1,1,0,1,0,1,0}],-CATInt["fam",{0,1,1,1,0,1,1,0}]+CATInt["fam",{1,1,0,1,1,1,0,0}],-CATInt["fam",{0,1,1,1,1,0,1,0}]+CATInt["fam",{1,1,0,1,1,0,1,0}],-CATInt["fam",{0,1,1,1,1,1,0,0}]+CATInt["fam",{1,1,0,1,0,1,1,0}],-CATInt["fam",{0,1,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,1,1,0,0}],-CATInt["fam",{0,1,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,0,1,1,0}],-CATInt["fam",{0,0,1,1,1,1,1,0}]+CATInt["fam",{1,0,0,1,1,1,1,0}],-CATInt["fam",{1,0,1,1,0,0,1,0}]-CATInt["fam",{1,0,1,1,1,0,0,0}]+2 CATInt["fam",{1,0,1,1,1,0,1,-1}]+s12 CATInt["fam",{1,0,1,1,1,0,1,0}],CATInt["fam",{1,0,1,1,1,-1,1,0}]-CATInt["fam",{1,0,1,1,1,0,1,-1}],2 CATInt["fam",{1,0,1,0,1,1,1,-1}]+s12 CATInt["fam",{1,0,1,0,1,1,1,0}],2 CATInt["fam",{1,0,1,-1,1,1,1,0}]+s12 CATInt["fam",{1,0,1,0,1,1,1,0}],CATInt["fam",{1,-1,1,1,1,0,1,0}]-CATInt["fam",{1,0,1,1,1,0,1,-1}],2 CATInt["fam",{1,-1,1,0,1,1,1,0}]+s12 CATInt["fam",{1,0,1,0,1,1,1,0}],-CATInt["fam",{0,0,1,1,1,0,1,0}]+CATInt["fam",{0,0,1,1,1,1,0,0}]+2 CATInt["fam",{0,0,1,1,1,1,1,-1}],CATInt["fam",{-1,1,1,1,1,1,0,0}]-CATInt["fam",{0,1,1,1,1,1,-1,0}],CATInt["fam",{-1,1,1,1,0,1,1,0}]-CATInt["fam",{0,1,1,1,-1,1,1,0}],CATInt["fam",{-1,1,0,1,1,1,1,0}]-CATInt["fam",{0,1,-1,1,1,1,1,0}],2 CATInt["fam",{-1,0,1,1,1,1,1,0}]+CATInt["fam",{0,0,1,1,1,1,0,0}]-s12 CATInt["fam",{0,0,1,1,1,1,1,0}],-2 CATInt["fam",{0,-1,1,1,1,1,1,0}]-CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{1,1,1,1,1,0,0,-1}],CATInt["fam",{-1,0,1,1,1,1,1,0}]+CATInt["fam",{0,0,1,1,1,1,0,0}]-s12 CATInt["fam",{0,0,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,1,0,-1,0}],CATInt["fam",{0,-1,1,1,1,1,1,0}]+CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{1,1,1,1,1,-1,0,0}]-CATInt["fam",{1,1,1,1,1,0,0,-1}],2 CATInt["fam",{0,1,-1,1,1,1,1,0}]-CATInt["fam",{0,1,0,1,1,1,0,0}]+CATInt["fam",{1,1,1,1,0,1,0,-1}],-CATInt["fam",{0,1,0,1,1,0,1,0}]+2 CATInt["fam",{0,1,0,1,1,1,1,-1}]-s23 CATInt["fam",{0,1,0,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,0,1,0,-1}],CATInt["fam",{-1,1,0,1,1,1,1,0}]-CATInt["fam",{0,1,0,1,0,1,1,0}]+CATInt["fam",{1,1,1,1,0,1,-1,0}]+CATInt["fam",{1,1,1,1,0,1,0,-1}],CATInt["fam",{1,1,1,1,0,0,1,-1}]-CATInt["fam",{1,1,1,1,1,0,0,-1}],CATInt["fam",{0,-1,1,1,1,1,1,0}]+CATInt["fam",{0,0,1,1,1,1,0,0}]+CATInt["fam",{1,1,1,1,0,-1,1,0}]-CATInt["fam",{1,1,1,1,1,0,0,-1}],CATInt["fam",{0,1,-1,1,1,1,1,0}]-CATInt["fam",{0,1,0,1,1,1,0,0}]+CATInt["fam",{1,1,1,1,-1,1,0,0}]+CATInt["fam",{1,1,1,1,0,1,0,-1}],CATInt["fam",{-1,0,1,1,1,1,1,0}]+CATInt["fam",{0,0,1,1,1,1,0,0}]-s12 CATInt["fam",{0,0,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,-1,0,1,0}],CATInt["fam",{1,-1,1,0,1,1,1,0}]+s12 CATInt["fam",{1,0,1,0,1,1,1,0}]+CATInt["fam",{1,1,1,0,1,-1,1,0}],CATInt["fam",{1,0,1,-1,1,1,1,0}]+s12 CATInt["fam",{1,0,1,0,1,1,1,0}]+CATInt["fam",{1,1,1,-1,1,0,1,0}],-CATInt["fam",{0,1,1,1,0,1,1,-1}]+CATInt["fam",{1,1,0,1,1,1,0,-1}],-CATInt["fam",{0,1,1,1,-1,1,1,0}]+CATInt["fam",{1,1,0,1,1,1,-1,0}],-CATInt["fam",{0,1,1,1,1,0,1,-1}]+CATInt["fam",{1,1,0,1,1,0,1,-1}],-CATInt["fam",{0,1,1,1,1,-1,1,0}]+CATInt["fam",{1,1,0,1,1,-1,1,0}],-CATInt["fam",{0,1,1,1,1,1,0,-1}]+CATInt["fam",{1,1,0,1,0,1,1,-1}],-CATInt["fam",{0,1,1,1,1,1,-1,0}]+CATInt["fam",{1,1,0,1,-1,1,1,0}],-CATInt["fam",{-1,1,1,1,0,1,1,0}]+CATInt["fam",{1,1,-1,1,1,1,0,0}],-CATInt["fam",{-1,1,1,1,1,0,1,0}]+CATInt["fam",{1,1,-1,1,1,0,1,0}],-CATInt["fam",{-1,1,1,1,1,1,0,0}]+CATInt["fam",{1,1,-1,1,0,1,1,0}],-CATInt["fam",{-1,1,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,1,1,-1,0}],CATInt["fam",{1,0,1,1,0,1,1,-1}]-CATInt["fam",{1,0,1,1,1,1,0,-1}],-CATInt["fam",{-1,1,1,1,1,0,1,0}]+CATInt["fam",{1,0,1,1,-1,1,1,0}],-CATInt["fam",{0,0,1,1,1,1,1,-1}]+CATInt["fam",{1,0,0,1,1,1,1,-1}],-CATInt["fam",{-1,0,1,1,1,1,1,0}]+CATInt["fam",{1,0,-1,1,1,1,1,0}],-CATInt["fam",{0,1,1,1,1,-1,1,0}]+CATInt["fam",{1,-1,1,1,1,1,0,0}],-CATInt["fam",{0,1,1,1,1,-1,1,0}]+CATInt["fam",{1,-1,1,1,0,1,1,0}],-CATInt["fam",{0,-1,1,1,1,1,1,0}]+CATInt["fam",{1,-1,0,1,1,1,1,0}],-CATInt["fam",{0,1,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,1,1,0,0}],-CATInt["fam",{1,0,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,1,0,1,0}],-CATInt["fam",{0,1,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,0,1,1,0}],-CATInt["fam",{0,1,1,1,1,1,1,0}]+CATInt["fam",{1,1,0,1,1,1,1,0}],-CATInt["fam",{-1,1,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,1,1,-1,0}],-CATInt["fam",{1,-1,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,1,-1,1,0}],CATInt["fam",{1,1,1,1,0,1,1,-1}]-CATInt["fam",{1,1,1,1,1,1,0,-1}],-CATInt["fam",{-1,1,1,1,1,1,1,0}]+CATInt["fam",{1,1,1,1,-1,1,1,0}],-CATInt["fam",{0,1,1,1,1,1,1,-1}]+CATInt["fam",{1,1,0,1,1,1,1,-1}],-CATInt["fam",{-1,1,1,1,1,1,1,0}]+CATInt["fam",{1,1,-1,1,1,1,1,0}]};


(* ::Text:: *)
(*Collect all equations and integrals*)


allequations = Join[equations,mappings];


intcases = CATIntCases[#,CATIntInvMonDefaultWeight]&;


ints = intcases[allequations];


(* ::Text:: *)
(*Define s,t,u as in section 6...*)


intS = -Plus@@Select[#[[2]],#<0&]&;
intT = Plus@@sectorof[#]&;
intU = ((Plus@@Select[#[[2]],#>0&]) - intT[#])&;


(* ::Text:: *)
(*...and use them to define a set of needed integrals.*)


needed = Select[ints,intU[#]<=1 && intS[#]<=2 &];


(* ::Text:: *)
(*Solve the system and compute DEs with FiniteFlow*)


extparams = {s12,s23};
params = Join[{d},extparams];


FFNewGraph["g","in",params]


FFAlgSparseSolver["g","ls",{"in"},params,Thread[allequations==0],ints,"NeededVars"->needed]


FFSolverSparseOutput["g","ls"]


FFGraphOutput["g","ls"]


learn = FFSparseSolverLearn["g",ints];


mis = {Complement[needed,"DepVars"/.learn],"IndepVars"/.learn}//intcases


mis//Length


(* ::Text:: *)
(*Having found a list of MIs, we need to differentiate them. In order to compute the differential operators, we need the full twist, including the external factor which we previously neglected, since the latter depends on the external invariants.*)


extpoly = -s12 s23 (s12+s23);
extexponent = 2-d/2;
fulltwist = CATMultiBaikov[Join[{extpoly},polys],Join[{extexponent},exponents],zs];


(* ::Text:: *)
(*When using the BaikovPackage, the external factor can be obtained by un-commenting the following block*)


(*
extpoly = baikovcomp[[1,1]];
extexponent = baikovcomp[[2,1]];
fulltwist = CATMultiBaikov[Join[{extpoly},polys],Join[{extexponent},exponents],zs];
*)


diffop = CATDiffOperator[fulltwist,zs,extparams,2,1,"KnownAnnihilatorSolutions"->ann];


(* ::Text:: *)
(*and we turn them into template identities*)


derivof = CATInvMonIdsFromDiffOperators["fam",diffop,zs];


(* ::Text:: *)
(*which we apply to the MIs*)


misderivs = Transpose[simplifyeqs[ Collect[derivof@@#[[2]]&/@mis,_CATInt,Together] ]];


(* ::Text:: *)
(*Now we follow the method 6.1 of section of https://arxiv.org/pdf/1905.08019 to directly reconstruct differential equations:*)


newneeded = intcases[misderivs];


FFSolverResetNeededVars["g","ls",ints,newneeded]


learn = FFSparseSolverLearn["g",ints];//AbsoluteTiming


newmis = Complement[newneeded,"DepVars"/.learn]//intcases


newmis == mis


FFSparseSolverMarkAndSweepEqs["g","ls"]


FFSparseSolverDeleteUnneededEqs["g","ls"]


(* ::Text:: *)
(*Note: here we could easily reconstruct reduction tables here (un-comment the following line to do that), but for complex examples we recommend reconstructing DEs directly.*)


(*reductions = FFSparseSolverSol[FFReconstructFunction["g",params],learn];*)


depints = "DepVars"/.learn;
alldints = Join[depints,mis];


dcoeffs = FFSparseRowRules[#,alldints]&/@(Join@@misderivs);


FFAlgRatFunEval["g","dcoeffs",{"in"},params,(Join@@dcoeffs)[[;;,2]]]


FFAlgRatNumEval["g","misred",ConstantArray[1,Length[mis]]]


FFAlgChain["g","fullred",{"ls","misred"}]


misno = Association[{}];
Do[misno[mis[[ii]]]=ii;,{ii,Length[mis]}];
redcols = Join[
misno/@#&/@("IndepVars"/.learn),
{#}&/@Range[Length[mis]]
];
Clear[misno];


FFAlgSparseMatMul["g","de",{"dcoeffs","fullred"},Length[dcoeffs],Length[alldints],Length[mis],
dcoeffs[[;;,;;,1]],
redcols
]


FFGraphOutput["g","de"]


rec = FFReconstructFunction["g",params,"PrintDebugInfo"->1];


des = ArrayReshape[rec,{Length[params]-1, Length[mis], Length[mis]}];


(* ::Text:: *)
(*DE matrix w.r.t. to s12*)


des[[1]]//Factor//MatrixForm


(* ::Text:: *)
(*DE matrix w.r.t. s23*)


des[[2]]//Factor//MatrixForm


(* ::Text:: *)
(*Check integrability conditions*)


D[des[[1]],s23]-D[des[[2]],s12] + (des[[1]] . des[[2]] - des[[2]] . des[[1]])//Factor
