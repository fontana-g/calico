(* ::Package:: *)

(* ::Text:: *)
(*This example computes the connection matrices (Omega) satisfied by the dual integrals in the different layers of the calculation of intersection numbers for the family of 2 - loop massless kite.*)


<<CALICO`


(* ::Text:: *)
(*Define loop denominators and the Baikov polynomial:*)


loopmoms = {k1,k2};
extmoms = {p};
dens = {{k1,0},{k1-p,0},{k1-k2,0},{k2,0},{k2-p,0}};
repl = {p^2 -> s};

{bpoly,zs} = CATBaikovPoly[
dens,
loopmoms,
extmoms,
repl,
z
];


(* ::Text:: *)
(*Define which denominators have to be regulated and compute the twist*)


regulated = {1,1,1,1,1};


CATDualRegulated["fam"] = regulated;


twist = CATBaikovDual[bpoly,d,Length[loopmoms],Length[extmoms],zs];


(* ::Text:: *)
(*Alternative way of computing the twist using CATMultiBaikov*)


(*(* alternative way *)
gamma = CATBaikovExponent[d,Length[loopmoms],Length[extmoms]];
twist = CATMultiBaikov[{bpoly},{-gamma},zs];*)


(* ::Text:: *)
(*Compute annihilators and template equations*)


ann = CATAnnihilator[twist,zs,2,1];//AbsoluteTiming


tmpids = CATDualInvMonIdsFromAnnihilators["fam",ann,zs];//AbsoluteTiming


(* ::Text:: *)
(*Choose the basis of the dual integrals for each layer as in ref . arXiv : 2304.14336*)


basis = {
{CATInt["fam",{1}],CATInt["fam",{0}]},
{CATInt["fam",{1,1}],CATInt["fam",{1,0}],CATInt["fam",{0,1}],CATInt["fam",{0,0}]},
{CATInt["fam",{1,1,1}],CATInt["fam",{0,1,1}],CATInt["fam",{1,0,1}],CATInt["fam",{1,1,0}]},
{CATInt["fam",{0,1,1,1}],CATInt["fam",{1,1,0,1}],CATInt["fam",{0,1,1,0}],CATInt["fam",{1,0,1,0}],CATInt["fam",{1,1,0,0}],CATInt["fam",{0,0,1,0}]}
};


sectors = {
{{0},{1}},
{{0, 0}, {0, 1}, {1, 0}, {1, 1}},
{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1,
   1, 0}, {1, 1, 1}},
{{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 1, 1}, {0, 1, 0, 
  0}, {0, 1, 0, 1}, {0, 1, 1, 0}, {0, 1, 1, 1}, {1, 0, 0, 0}, {1, 0, 
  0, 1}, {1, 0, 1, 0}, {1, 0, 1, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 
  1, 1, 0}, {1, 1, 1, 1}}
};


(* ::Text:: *)
(*Compute Omega*)


If[$FrontEnd=!=Null, SetDirectory[NotebookDirectory[]]];
Get["loop_duals_common.m"];


Omega1 = ComputeOmega[1];


Omega2 = ComputeOmega[2];


Omega3 = ComputeOmega[3];


Omega4 = ComputeOmega[4];


MatrixForm/@{
Omega1,
Omega2,
Omega3,
Omega4
}
