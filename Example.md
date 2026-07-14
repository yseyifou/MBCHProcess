# MBCHProcess Example

This example demonstrates a complete, real-world workflow using `MBCHProcess` to identify the hypergeometric structure of a Feynman integral. We compute the two-loop sunset diagram $H_{\{1,2,3\}}(m, M, M; m^2)$ with $p^2 = m^2$.

## Prerequisites

Ensure that `AMBRE`, `MultivariateResidues`, `MBConicHulls`, `HypExp`, and `MBCHProcess` are loaded into your Mathematica environment. For reproducibility, the loading sequence and the diagram definition are shown below.

## Mathematica Workflow

Open a Mathematica notebook and evaluate the following cells:

### 1. Load Packages

Load the packages in the correct order to ensure all dependencies are met. (Adjust the file paths according to your local directory structure).

```mathematica
Get["/path/to/AMBREv2.1.1.m"];
Get["/path/to/MultivariateResidues.m"];
Get["/path/to/MBConicHulls.wl"];
Get["/path/to/HypExp.m"];
Get["/path/to/MBCHProcess.m"];
```

### 2. Define the Diagram and Generate MB Representation

Use `AMBRE` to construct the Mellin-Barnes representation for the sunset diagram.

```mathematica
(* Define the sunset diagram using AMBRE *)
(* H{1,2,3}(m,M,M;m^2) with p^2=m^2 *)
Propagators = {PR[k1 - p, m, a1]*PR[k2, M, a2]*PR[k1 - k2, M, a3]};
invariants = {p^2 -> m^2};
PreFactor = {1};

(* Generate the raw Mellin-Barnes representation *)
mbRep = MBrepr[PreFactor, Propagators, {k2, k1}, Text -> False, BarnesLemma -> True];

(* Substitute specific indices and simplify *)
Simplify[mbRep /. {a1 -> 1, a2 -> 2, a3 -> 3, eps -> e}]
```

### 3. Resolve the MB Representation

Format the output for `MBConicHulls` and use `ResolveMB` to geometrically identify the series expansions.

```mathematica
(* Format the representation for MBConicHulls *)
Rep = MBRep[
  1/2 (m^2)^(-2 e - 1), 
  {z1}, 
  {M^2/m^2}, 
  {{-3 - 4 e - 2 z1, -1 - e - z1, -e - z1, -z1, 3 + e + z1, 2 + 2 e + z1}, 
   {-1 - 2 e - 2 z1, -3 e - z1}}
];

(* Resolve the MB representation *)
test = ResolveMB[Rep]
```

*Note: `ResolveMB` will output details about the conic hulls and the series solutions found.*

### 4. Process the Series with MBCHProcess

Finally, pass the resolved representation to `ProcessSeries` to convert the raw Gamma function sums into standard special functions.

```mathematica
(* Process the first series solution *)
ProcessSeries[test, 1]
```

## Expected Output

`MBCHProcess` will evaluate the series and return a compact analytical expression. For this specific sunset configuration, the output is a combination of generalized hypergeometric functions (${_4F_3}$) with their respective prefactors:

```mathematica
-(((-1 + e) e (-1 + 2 e) m^(-2 - 2 e) (1/M^2)^(3 + e) Gamma[5 - 2 e] Gamma[-1 + e] Gamma[3 + e] 
    (m^2 M^4 - e m^2 M^4 - 4 M^6 - 4 e M^6 + 
     4 M^6 (1 - e) HypergeometricPFQ[{1, 1, 2 - 2 e, e}, {4, -1 - e, -2 e}, m^2/M^2] - 
     2 m^2 M^4 (1 - e) HypergeometricPFQ[{1, 1, 3 - 2 e, 1 + e}, {5, 1 - 2 e, -e}, m^2/M^2])
   )/(48 (-2 + e) (1 + e) (2 + e) (-3 + 2 e) Gamma[3 - 2 e])) 
+ 
(m^2 (1/M^4)^(1 + e) Gamma[1 - e] Gamma[1 + e] Gamma[2 + e] Gamma[2 + 2 e] 
  HypergeometricPFQ[{3, 1 + e, 2 + e, 2 + 2 e}, {2 - e, e, 5 + 2 e}, m^2/M^2]
 )/(Gamma[2 - e] Gamma[5 + 2 e])
```
