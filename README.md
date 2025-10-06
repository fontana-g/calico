# CALICO

CALICO is a Mathematica package based on the following publication:
- Giuseppe Bertolini, Gaia Fontana, Tiziano Peraro, *CALICO: Computing Annihilators from Linear Identities Constraining (differential) Operators*, [JHEP 10 (2025) 018](https://link.springer.com/article/10.1007/JHEP10(2025)018) [[arXiv:2506.13653]](https://arxiv.org/abs/2506.13653).

CALICO contains various routines for computing parametric annihilators and other differential operators, and using them to find identities between special functions with a suitable integral representation. It contains specialized routines for various representations of loop integrals, as well as more general routines which can be adapted to custom problems. Functions for solving syzygy equations and polynomial decomposition problems are also provided.

<p align="center"> <img src="https://raw.githubusercontent.com/fontana-g/calico-logo/refs/heads/main/logo.png" width="350" alt="CALICO logo"> </p>

## Installation

CALICO is a [Mathematica](https://www.wolfram.com/mathematica/) package and depends on the latest version of [FiniteFlow](https://github.com/peraro/finiteflow), which must be installed or updated first.

To install the package, either copy the file `CALICO.m` in a directory in your Mathematica `$Path` or add the following to your `init.m` file
```
$CalicoPath = "/path/to/calico"; (* directory containing CALICO.m )
If[Not[MemberQ[$Path,$CalicoPath]], $Path = Flatten[{$Path, $CalicoPath}]];
```

## Usage

The most important functions of the package are documented in section 10 of [the paper](https://arxiv.org/abs/2506.13653). Examples can be found in the [examples/](examples/) subdirectory.


## Reporting issues

The preferred way of reporting issues is by contacting the maintainers Gaia Fontana and Tiziano Peraro via email. Before reporting an issue, please make sure that your CALICO and FiniteFlow installations are both up to date.

## Credits

 Package authors:
 - Giuseppe Bertolini
 - Gaia Fontana
 - Tiziano Peraro

 Logo:
 - Gaia Fontana - [@qftoons](https://www.instagram.com/qftoons/)
