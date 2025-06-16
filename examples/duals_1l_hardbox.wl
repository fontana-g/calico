(* ::Package:: *)

(* ::Text:: *)
(*This example computes the connection matrices (Omega) satisfied by the dual integrals in the different layers of the calculation of intersection numbers for the family of 1-loop box with two off-shell legs.*)


<<CALICO`


(* ::Text:: *)
(*Define loop denominators and the Baikov polynomial :*)


loopmoms = {k1};
extmoms = {p1,p2,p3};
dens = {{k1,0},{k1-p1,0},{k1-p1-p2,0},{k1-p1-p2-p3,0}};
repl = {p3^2 -> p3sq, p2^2 -> 0, p1^2 -> 0, p1 p3 -> -(p3sq/2) + t/2, 
 p1 p2 -> s/2, p2 p3 -> p4sq/2 - s/2 - t/2};
 
{bpoly,zs} = CATBaikovPoly[
dens,
loopmoms,
extmoms,
repl,
z
];


(* ::Text:: *)
(*Define which denominators have to be regulated and compute the twist*)


regulated = {1,1,1,1};


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


tmpeqs = CATDualInvMonIdsFromAnnihilators["fam",ann,zs];//AbsoluteTiming


(* ::Text:: *)
(*Choose the basis of the dual integrals for each layer as in ref. arXiv:2304.14336*)


basis = {
{CATInt["fam",{1}],CATInt["fam",{0}]},
{CATInt["fam",{1,1}],CATInt["fam",{0,1}],CATInt["fam",{1,0}],CATInt["fam",{0,0}]},
{CATInt["fam",{1,1,1}],CATInt["fam",{1,0,1}],CATInt["fam",{0,0,1}],CATInt["fam",{0,1,0}],CATInt["fam",{1,0,0}],CATInt["fam",{0,0,0}]}
};


sectors = {
{{0},{1}},
{{0,0},{0,1},{1,0},{1,1}},
{{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}}
};


(* ::Text:: *)
(*ComputeOmega*)


If[$FrontEnd=!=Null, SetDirectory[NotebookDirectory[]]];
Get["loop_duals_common.m"];


Omega1 = ComputeOmega[1];


Omega2 = ComputeOmega[2];


Omega3 = ComputeOmega[3];


MatrixForm/@{
Omega1,
Omega2,
Omega3
}
