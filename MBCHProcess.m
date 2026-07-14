(* :Title: MBCHProcess *)
(* :Context: MBCHProcess *)
(* :Author: YS, QU, EM *)
(* :Summary: Post-processing of MBConicHulls output (EvaluateSeries function): automatically identifies, term by term, the expression in terms of known special functions (HypergeometricPFQ, AppellF1, generic Horn functions, Lauricella A/B/C/D, Kampe de Feriet). *)
(* :Requires: MBConicHulls must be loaded before this package (the EvaluateSeries function must be present in the Global or MBConicHulls context). *)

BeginPackage["MBCHProcess`"]

Unprotect[ProcessSeries, KampeDeFeriet, LauricellaA, LauricellaB, 
  LauricellaC, LauricellaD, HornFunction];

Clear[ProcessSeries, KampeDeFeriet, LauricellaA, LauricellaB, 
  LauricellaC, LauricellaD, HornFunction];

ProcessSeries::usage = 
  "ProcessSeries[rep, n, Evaluate :> False, OptionsPattern[]] processes \
the resolved MB representation rep (from ResolveMB) and identifies \
the n-th series term, returning an expression in terms of known \
special functions. Options: Verbose -> True/False (default False).";
ProcessSeries::noeval = 
  "EvaluateSeries not found. Please load MBConicHulls first.";

KampeDeFeriet::usage = 
  "KampeDeFeriet[a, b, bp, c, d, dp, x, y] represents the Kampe de Feriet \
function.";

LauricellaA::usage = 
  "LauricellaA[a, b, c, vars] represents the Lauricella A function.";

LauricellaB::usage = 
  "LauricellaB[a, b, c, vars] represents the Lauricella B function.";

LauricellaC::usage = 
  "LauricellaC[a, b, c, vars] represents the Lauricella C function.";

LauricellaD::usage = 
  "LauricellaD[a, b, c, vars] represents the Lauricella D function.";

HornFunction::usage = 
  "HornFunction[name, P, Q, vars, numIndices, demIndices] represents \
a generic Horn series.";

Options[ProcessSeries] = {Verbose -> False};

Begin["`Private`"]

(* ==========Import confirmation message==========*)
Print["========================================="];
Print["  MBCHProcess v1.0 "];
Print["  ProcessSeries: Identify special functions from MB \
representations"];
Print["  Usage: ProcessSeries[rep, n, Evaluate :> False, \
Verbose -> True/False]"];
Print["  Dependencies: MBConicHulls must be loaded first."];
Print["========================================="];

(* ==========Check that EvaluateSeries is available==========*)
If[! MemberQ[Names["MBConicHulls`*"], "MBConicHulls`EvaluateSeries"], 
  Message[ProcessSeries::noeval]];

(* ==========Special functions: display and numerical evaluation==========*)

(*KampeDeFeriet*)
KampeDeFeriet /: 
  MakeBoxes[
   KampeDeFeriet[a_List, b_List, bPrime_List, c_List, d_List, 
    dPrime_List, x_, y_], StandardForm] := 
  Module[{superscript, subscript, aBox, bBox, bPrimeBox, cBox, dBox, 
    dPrimeBox, numRow, denRow, matrix, varBox}, 
   superscript = 
    ToString[Length[a]] <> ":" <> ToString[Length[b]] <> ";" <> 
      ToString[Length[bPrime]];
   subscript = 
    ToString[Length[c]] <> ":" <> ToString[Length[d]] <> ";" <> 
      ToString[Length[dPrime]];
   aBox = If[a === {}, "\[Dash]", RowBox[Riffle[Map[MakeBoxes, a], ","]]];
   bBox = If[b === {}, "\[Dash]", RowBox[Riffle[Map[MakeBoxes, b], ","]]];
   bPrimeBox = 
    If[bPrime === {}, "\[Dash]", 
     RowBox[Riffle[Map[MakeBoxes, bPrime], ","]]];
   cBox = If[c === {}, "\[Dash]", RowBox[Riffle[Map[MakeBoxes, c], ","]]];
   dBox = If[d === {}, "\[Dash]", RowBox[Riffle[Map[MakeBoxes, d], ","]]];
   dPrimeBox = 
    If[dPrime === {}, "\[Dash]", 
     RowBox[Riffle[Map[MakeBoxes, dPrime], ","]]];
   numRow = RowBox[{aBox, " : ", bBox, " ; ", bPrimeBox}];
   denRow = RowBox[{cBox, " : ", dBox, " ; ", dPrimeBox}];
   matrix = 
    GridBox[{{numRow}, {denRow}}, ColumnAlignments -> Left, 
     RowSpacings -> 0, ColumnSpacings -> 1];
   varBox = RowBox[{MakeBoxes[x], " , ", MakeBoxes[y]}];
   With[{sub = subscript, sup = superscript}, 
    RowBox[{SubsuperscriptBox["F", sub, sup], "[", matrix, 
      "\[VerticalLine]", varBox, "]"}]]];

KampeDeFeriet /: 
  N[KampeDeFeriet[a_List, b_List, bPrime_List, c_List, d_List, 
    dPrime_List, x_, y_]] := 
  Module[{p, q, r, s, t, u, maxM = 30, maxN = 30, term, sum = 0., 
    tol = 10^-12}, p = Length[a]; q = Length[b]; r = Length[bPrime];
   s = Length[c]; t = Length[d]; u = Length[dPrime];
   For[m = 0, m <= maxM, m++, 
    For[n = 0, n <= maxN, n++, 
     term = N[(Times @@ Pochhammer[a, m + n])*(Times @@ 
          Pochhammer[b, m])*(Times @@ 
           Pochhammer[bPrime, 
            n])/((Times @@ Pochhammer[c, m + n])*(Times @@ 
             Pochhammer[d, m])*(Times @@ Pochhammer[dPrime, n]))*(x^m/
          m!)*(y^n/n!)];
     sum += term;
     If[Abs[term] < tol && m > 5 && n > 5, Break[]];];
    If[Abs[term] < tol && m > 5, Break[]];];
   sum];

(*Lauricella A*)
LauricellaA /: 
  MakeBoxes[LauricellaA[a_, b_List, c_List, vars_List], 
   StandardForm] := 
  Module[{bArgs, cArgs, varArgs}, 
   bArgs = If[b === {}, "_", RowBox[Riffle[Map[MakeBoxes, b], ","]]];
   cArgs = If[c === {}, "_", RowBox[Riffle[Map[MakeBoxes, c], ","]]];
   varArgs = 
    If[vars === {}, "_", RowBox[Riffle[Map[MakeBoxes, vars], ","]]];
   RowBox[{"F_A", "[", MakeBoxes[a], "; ", bArgs, "; ", cArgs, "; ", 
     varArgs, "]"}]];

LauricellaA /: N[LauricellaA[a_, b_List, c_List, vars_List]] := 
  Module[{n = Length[vars], ks}, ks = Array[k, n];
   NSum[Pochhammer[a, Total[ks]]*
     Product[Pochhammer[b[[i]], k[i]]/Pochhammer[c[[i]], k[i]]*
       vars[[i]]^k[i]/k[i]!, {i, n}], 
    Evaluate[Sequence @@ Table[{k[i], 0, Infinity}, {i, n}]], 
    Method -> "WynnEpsilon", NSumTerms -> 50]];

(*Lauricella B*)
LauricellaB /: 
  MakeBoxes[LauricellaB[a_List, b_List, c_, vars_List], 
   StandardForm] := 
  Module[{aArgs, bArgs, varArgs}, 
   aArgs = If[a === {}, "_", RowBox[Riffle[Map[MakeBoxes, a], ","]]];
   bArgs = If[b === {}, "_", RowBox[Riffle[Map[MakeBoxes, b], ","]]];
   varArgs = 
    If[vars === {}, "_", RowBox[Riffle[Map[MakeBoxes, vars], ","]]];
   RowBox[{"F_B", "[", aArgs, "; ", bArgs, "; ", MakeBoxes[c], "; ", 
     varArgs, "]"}]];

LauricellaB /: N[LauricellaB[a_List, b_List, c_, vars_List]] := 
  Module[{n = Length[vars], ks}, ks = Array[k, n];
   NSum[Product[
      Pochhammer[a[[i]], k[i]]*Pochhammer[b[[i]], k[i]]*
       vars[[i]]^k[i]/k[i]!, {i, n}]/Pochhammer[c, Total[ks]], 
    Evaluate[Sequence @@ Table[{k[i], 0, Infinity}, {i, n}]], 
    Method -> "WynnEpsilon", NSumTerms -> 50]];

(*Lauricella C*)
LauricellaC /: 
  MakeBoxes[LauricellaC[a_, b_, c_List, vars_List], StandardForm] := 
  Module[{numRow, denRow, matrix, varBox}, 
   numRow = RowBox[{MakeBoxes[a], ",", MakeBoxes[b]}];
   denRow = If[c === {}, "\[Dash]", RowBox[Riffle[Map[MakeBoxes, c], ","]]];
   matrix = 
    GridBox[{{numRow}, {denRow}}, ColumnAlignments -> Left, 
     RowSpacings -> 0, ColumnSpacings -> 1];
   varBox = 
    If[vars === {}, "\[Dash]", RowBox[Riffle[Map[MakeBoxes, vars], ","]]];
   RowBox[{SubscriptBox["F", "C"], "[", matrix, "\[VerticalLine]", 
     varBox, "]"}]];

LauricellaC /: N[LauricellaC[a_, b_, c_List, vars_List]] := 
  Module[{n = Length[vars], ks}, ks = Array[k, n];
   NSum[Pochhammer[a, Total[ks]]*Pochhammer[b, Total[ks]]*
     Product[1/Pochhammer[c[[i]], k[i]]*vars[[i]]^k[i]/k[i]!, {i, n}], 
    Evaluate[Sequence @@ Table[{k[i], 0, Infinity}, {i, n}]], 
    Method -> "WynnEpsilon", NSumTerms -> 50]];

(*Lauricella D*)
LauricellaD /: 
  MakeBoxes[LauricellaD[a_, b_List, c_, vars_List], StandardForm] := 
  Module[{bArgs, varArgs}, 
   bArgs = If[b === {}, "_", RowBox[Riffle[Map[MakeBoxes, b], ","]]];
   varArgs = 
    If[vars === {}, "_", RowBox[Riffle[Map[MakeBoxes, vars], ","]]];
   RowBox[{"F_D", "[", MakeBoxes[a], "; ", bArgs, "; ", MakeBoxes[c], 
     "; ", varArgs, "]"}]];

LauricellaD /: N[LauricellaD[a_, b_List, c_, vars_List]] := 
  Module[{n = Length[vars], ks}, ks = Array[k, n];
   NSum[Pochhammer[a, Total[ks]]*
     Product[Pochhammer[b[[i]], k[i]]*vars[[i]]^k[i]/k[i]!, {i, n}]/
      Pochhammer[c, Total[ks]], 
    Evaluate[Sequence @@ Table[{k[i], 0, Infinity}, {i, n}]], 
    Method -> "WynnEpsilon", NSumTerms -> 50]];

(*Generic HornFunction*)
HornFunction /: 
  MakeBoxes[HornFunction[name_String, P_List, Q_List, vars_List, ___], 
    StandardForm] := 
  Module[{pBox, qBox, argsBox}, 
   pBox = If[P === {}, "_", RowBox[Riffle[Map[MakeBoxes, P], ","]]];
   qBox = If[Q === {}, "_", RowBox[Riffle[Map[MakeBoxes, Q], ","]]];
   argsBox = 
    If[vars === {}, "_", RowBox[Riffle[Map[MakeBoxes, vars], ","]]];
   RowBox[{name, "[", pBox, " ; ", qBox, " | ", argsBox, "]"}]];

HornFunction /: 
  N[HornFunction[name_String, P_List, Q_List, vars_List, 
    numIndices_List, demIndices_List]] := 
  Module[{x, y, rep, numTerm, denTerm, term}, {x, y} = vars;
   rep = {Subscript[n, 1] -> m, Subscript[n, 2] -> n};
   numTerm = 
    If[Length[P] > 0, 
     Times @@ 
      Table[Pochhammer[P[[i]], (numIndices[[i]] /. rep)], {i, 
        Length[P]}], 1];
   denTerm = 
    If[demIndices === {0} || Length[Q] == 0, 1, 
     Times @@ 
      Table[Pochhammer[Q[[i]], (demIndices[[i]] /. rep)], {i, 
        Length[Q]}]];
   term = numTerm/denTerm*x^m/m!*y^n/n!;
   NSum[term, {m, 0, Infinity}, {n, 0, Infinity}, 
    Method -> "WynnEpsilon", NSumTerms -> 50]];

(*Special case: Appell F1*)
HornFunction /: 
  N[HornFunction["F1", P_List, Q_List, vars_List, numIndices_List, 
    demIndices_List]] := 
  N[AppellF1[P[[1]], P[[2]], P[[3]], Q[[1]], vars[[1]], vars[[2]]]];

HornFunction::nimp = 
  "Numerical evaluation of Horn series `` is not yet implemented.";

(* ==========Auxiliary functions for ProcessSeries==========*)

(*Factorization in Gamma*)
FactorGammaPositive[expr_] := 
  expr /. Gamma[w_] :> 
    Module[{factored, term}, factored = FactorTerms[w];
     If[MatchQ[factored, Times[_Integer, _]], term = factored[[1]];
      If[term < 0, Gamma[-term*(-factored[[2]])], Gamma[factored]], 
      Gamma[w]]];

(*Gauss multiplication formula*)
GMultiplication[expr_] := 
  expr /. Gamma[
     Times[n_Integer?Positive, 
      w_]] :> (2 Pi)^((1 - n)/2) n^(n w - 1/2) Product[
      Gamma[w + k/n], {k, 0, n - 1}];

(*Reflection formula*)
Reflection[expr_] := 
  expr /. Gamma[z_ - Subscript[n_, i_Integer]] :> 
    Gamma[z] Gamma[1 - z] (-1)^Subscript[n, i]/
      Gamma[Subscript[n, i] + 1 - z];

(*Gamma to Pochhammer replacement*)
ReplacePochhammer[expr_] := 
  expr /. Gamma[args___] :> 
    With[{nList = Cases[{args}, Subscript[n, _Integer], Infinity], 
      restList = 
        DeleteCases[{args}, Subscript[n, _Integer], Infinity]}, {x, 
       nSum} = {Total[restList], Total[nList]};
      If[nList != {}, Gamma[x]*Pochhammer[x, nSum], Gamma[args]]];

(*Separation of exponents according to indices*)
SplitExponent[expr_, k_Integer, IndexList_] := 
  Module[{base, exp, terms, const, indices = IndexList[[k]], coeff}, 
   If[MatchQ[expr, Power[_, _]], base = expr[[1]];
    exp = Expand[expr[[2]]];
    terms = If[Head[exp] === Plus, List @@ exp, {exp}];
    const = Select[terms, FreeQ[Alternatives @@ indices]];
    Prepend[
      Table[coeff = Plus @@ (Coefficient[#, x] & /@ terms); 
       Inactive[Power][base, coeff*x], {x, indices}], 
      Inactive[Power][base, Plus @@ const]], 
    Prepend[Table[1, {Length[indices]}], Inactive[Power][expr, 1]]]];

(*Gamma formatting*)
FormatGamma[expr_, k_Integer, IndexList_] := 
  Module[{L, Start, result, withN, withnoN, final}, 
   L = Apply[List, expr];
   Start = Table[1, {Length[IndexList[[k]]] + 1}];
   result = 
    Simplify[
      PowerExpand[
       Activate[Fold[#1*SplitExponent[#2, k, IndexList] &, Start, L]]]];
   withN = Select[result[[1]], ! FreeQ[#, Subscript[n, _]] &];
   withnoN = Select[result[[1]], FreeQ[#, Subscript[n, _]] &];
   (*final = Join[{withnoN}, {withN}, Rest[result]];*)
   final = {withnoN, withN, Rest[result]};
   Return[final]];

(*Extraction of (a, n) pairs*)
ExtractPochhammerCouples[expr_] := 
  Flatten[Cases[
    expr, (Pochhammer[d_, f_]^n_.) :> Table[{d, f}, n], {0, 1}], 1];

(*Grouping by index combination*)
IdentificationByIndex[expr_, IndexList_, k_Integer] := 
  Module[{adjustedExpr, numCouples, denCouples, allIndexCombos, 
    result}, 
   adjustedExpr = expr*Product[Pochhammer[1, ni], {ni, IndexList[[k]]}];
   numCouples = ExtractPochhammerCouples[Numerator[adjustedExpr]];
   denCouples = ExtractPochhammerCouples[Denominator[adjustedExpr]];
   allIndexCombos = 
    DeleteDuplicates[Join[numCouples[[All, 2]], denCouples[[All, 2]]]];
   result = 
    Association[
      Table[combo -> <|
         "NumeratorArgs" -> Cases[numCouples, {a_, combo} :> a], 
         "DenominatorArgs" -> 
          Cases[denCouples, {a_, combo} :> a]|>, {combo, 
        allIndexCombos}]];
   Return[result]];

(*Function recognition*)
Recognize[rawNum_, rawDem_, IndexList_, k_] := 
  Module[{Nv, standardize, nNum, nDem, uniqueNum, uniqueDem, singles, 
    total}, Nv = Length[IndexList[[k]]];
   standardize = 
    Thread[IndexList[[k]] -> Array[Symbol["n" <> ToString[#]] &, Nv]];
   nNum = Sort[rawNum /. standardize];
   nDem = Sort[DeleteCases[rawDem /. standardize, 0]];
   uniqueNum = Union[nNum];
   uniqueDem = Union[nDem];
   singles = Array[Symbol["n" <> ToString[#]] &, Nv];
   total = Total[singles];
   Which[(*Horn with 2 variables*)
    nNum === Sort[{n1, n2, n1 + n2}] && nDem === Sort[{n1 + n2}], "F1", 
     nNum === Sort[{n1, n2, n1 + n2}] && nDem === Sort[{n1, n2}], "F2", 
     nNum === Sort[{n1, n1, n2, n2}] && nDem === Sort[{n1 + n2}], "F3", 
     nNum === Sort[{n1 + n2, n1 + n2}] && nDem === Sort[{n1, n2}], 
    "F4", nNum === Sort[{n1 + n2, n2 - n1, n1 - n2}] && nDem === {}, 
    "G1", nNum === Sort[{n1, n2, n2 - n1, n1 - n2}] && nDem === {}, 
    "G2", nNum === Sort[{2 n2 - n1, 2 n1 - n2}] && nDem === {}, "G3", 
    nNum === Sort[{n1 - n2, n1 + n2, n2}] && nDem === Sort[{n1}], "H1", 
     nNum === Sort[{n1 - n2, n1, n2, n2}] && nDem === Sort[{n1}], "H2", 
     nNum === Sort[{2 n1 + n2, n2}] && nDem === Sort[{n1 + n2}], "H3", 
     nNum === Sort[{2 n1 + n2, n2}] && nDem === Sort[{n1, n2}], "H4", 
    nNum === Sort[{2 n1 + n2, n2 - n1}] && nDem === Sort[{n2}], "H5", 
    nNum === Sort[{2 n1 - n2, n2 - n1, n2}] && nDem === {}, "H6", 
    nNum === Sort[{2 n1 - n2, n2, n2}] && nDem === Sort[{n1}], 
    "H7",(*Lauricella N-variables*)
    nNum === Sort[Append[singles, total]] && nDem === Sort[{total}], 
    "LauricellaD", 
    nNum === Sort[Append[singles, total]] && nDem === Sort[singles], 
    "LauricellaA", 
    nNum === Sort[Join[singles, singles]] && nDem === Sort[{total}], 
    "LauricellaB", 
    nNum === Sort[{total, total}] && nDem === Sort[singles], 
    "LauricellaC",(*KampeDeFeriet*)
    SubsetQ[{n1, n2, n1 + n2}, uniqueNum] && 
     SubsetQ[{n1, n2, n1 + n2}, uniqueDem], "KampeDeFeriet", True, 
    "UnknownFunction"]];

(*Association for a term*)
FunctionIdentification[tab_, IndexList_, k_] := 
  Module[{FunctionList, idRes, keys, numIndices, demIndices, pArgs, 
    qArgs, seriesArgs}, idRes = tab[[2]];
   keys = Keys[idRes];
   pArgs = Flatten[Values[idRes[[All, "NumeratorArgs"]]]];
   qArgs = Flatten[Values[idRes[[All, "DenominatorArgs"]]]];
   seriesArgs = 
    Table[(tab[[3]][[i]] /. IndexList[[k, i]] -> 1), {i, 
      Length[tab[[3]]]}];
   FunctionList = <|"Prefactor" -> tab[[1]], "PArgs" -> pArgs, 
     "QArgs" -> qArgs, "IdRes" -> idRes, "Indices" -> IndexList[[k]]|>;
   If[Length[keys] <= 1, FunctionList["SeriesArg"] = seriesArgs[[1]], 
    numIndices = 
     Flatten[Table[
       Table[idx, Length[idRes[idx, "NumeratorArgs"]]], {idx, keys}]];
     demIndices = 
      Flatten[Table[
        Table[idx, Length[idRes[idx, "DenominatorArgs"]]], {idx, keys}]];
     If[demIndices === {}, demIndices = {0}];
     FunctionList["NumIndices"] = numIndices;
     FunctionList["DemIndices"] = demIndices;
     FunctionList["SeriesArgs"] = seriesArgs;
     FunctionList["Name"] = 
      Recognize[numIndices, demIndices, IndexList, k];];
   Return[FunctionList]];

(*Display for 1D series*)
Display1[L_List] := 
  Sum[Module[{coef, P, Q, x}, coef = item["Prefactor"];
    P = Flatten[{item["PArgs"]}];
    Q = Flatten[{item["QArgs"]}];
    x = item["SeriesArg"];
    coef*HypergeometricPFQ[P, Q, x]], {item, L}];

(*Display for 2D series*)
Display2[L_List] := 
  Sum[Module[{coef, name, P, Q, idRes, vars, indices, keysAll, 
      keysSingle, a, b, bList, cList, b1List, b2List, c, keysA, keysB, 
      keysBPrime, aKdF, bKdF, bPrimeKdF, cKdF, dKdF, dPrimeKdF}, 
    indices = item["Indices"];
    coef = item["Prefactor"];
    name = item["Name"];
    P = Flatten[{item["PArgs"]}];
    Q = Flatten[{item["QArgs"]}];
    idRes = item["IdRes"];
    vars = item["SeriesArgs"];
    Which[name === "F1", 
     coef*AppellF1[P[[1]], P[[2]], P[[3]], Q[[1]], vars[[1]], vars[[2]]], 
      name === "KampeDeFeriet", 
     keysA = Select[
       Keys[idRes], ! FreeQ[#, Subscript[n, 1]] && ! 
          FreeQ[#, Subscript[n, 2]] &];
     keysB = 
      Select[Keys[
        idRes], ! FreeQ[#, Subscript[n, 1]] && 
         FreeQ[#, Subscript[n, 2]] &];
     keysBPrime = 
      Select[Keys[idRes], 
       FreeQ[#, Subscript[n, 1]] && ! FreeQ[#, Subscript[n, 2]] &];
     aKdF = 
      Flatten[Table[Lookup[idRes[k], "NumeratorArgs", {}], {k, keysA}]];
     bKdF = 
      Flatten[Table[Lookup[idRes[k], "NumeratorArgs", {}], {k, keysB}]];
     bPrimeKdF = 
      Flatten[Table[
        Lookup[idRes[k], "NumeratorArgs", {}], {k, keysBPrime}]];
     cKdF = 
      Flatten[Table[
        Lookup[idRes[k], "DenominatorArgs", {}], {k, keysA}]];
     dKdF = 
      Flatten[Table[
        Lookup[idRes[k], "DenominatorArgs", {}], {k, keysB}]];
     dPrimeKdF = 
      Flatten[Table[
        Lookup[idRes[k], "DenominatorArgs", {}], {k, keysBPrime}]];
     coef*
      KampeDeFeriet[aKdF, bKdF, bPrimeKdF, cKdF, dKdF, dPrimeKdF, 
       vars[[1]], vars[[2]]], 
     MemberQ[{"LauricellaD", "LauricellaA", "LauricellaB", 
       "LauricellaC"}, name], 
     keysAll = Select[Keys[idRes], # === Total[indices] &];
     keysSingle = 
      SortBy[Select[Keys[idRes], Head[#] === Subscript &], #[[2]] &];
     Which[name === "LauricellaD", 
      a = Flatten[
        Table[Lookup[idRes[k], "NumeratorArgs", {}], {k, keysAll}]];
      bList = 
       Table[Flatten[
         Table[Lookup[idRes[k], 
           "NumeratorArgs", {}], {k, {ks}}]], {ks, keysSingle}];
      c = 
       Flatten[Table[
         Lookup[idRes[k], "DenominatorArgs", {}], {k, keysAll}]];
      coef*LauricellaD[a, bList, c, vars], name === "LauricellaA", 
      a = 
       Flatten[Table[
         Lookup[idRes[k], "NumeratorArgs", {}], {k, keysAll}]];
      bList = 
       Table[Flatten[
         Table[Lookup[idRes[k], 
           "NumeratorArgs", {}], {k, {ks}}]], {ks, keysSingle}];
      cList = 
       Table[Flatten[
         Table[Lookup[idRes[k], 
           "DenominatorArgs", {}], {k, {ks}}]], {ks, keysSingle}];
      coef*LauricellaA[a, bList, cList, vars], name === "LauricellaB", 
       b1List = 
       Table[Flatten[
          Table[Lookup[idRes[k], 
            "NumeratorArgs", {}], {k, {ks}}]][[1]], {ks, keysSingle}];
       b2List = 
       Table[Flatten[
          Table[Lookup[idRes[k], 
            "NumeratorArgs", {}], {k, {ks}}]][[2]], {ks, keysSingle}];
       c = 
       Flatten[Table[
         Lookup[idRes[k], "DenominatorArgs", {}], {k, keysAll}]];
       coef*LauricellaB[b1List, b2List, c, vars], 
      name === "LauricellaC", 
      a = Flatten[
         Table[Lookup[idRes[k], "NumeratorArgs", {}], {k, 
           keysAll}]][[1]];
      b = 
       Flatten[Table[
          Lookup[idRes[k], "NumeratorArgs", {}], {k, keysAll}]][[2]];
      cList = 
       Table[Flatten[
         Table[Lookup[idRes[k], 
           "DenominatorArgs", {}], {k, {ks}}]], {ks, keysSingle}];
      coef*LauricellaC[a, b, cList, vars]], True, 
     With[{nInd = item["NumIndices"], dInd = item["DemIndices"]}, 
      coef*HornFunction[name, P, Q, vars, nInd, dInd]]]], {item, L}];

(*Prefactor simplification*)
ProcessPrefactor[pre_] := Module[{simpler}, simpler = PowerExpand[pre];
   simpler = Simplify[simpler];
   simpler = Simplify[simpler];
   simpler];

(* ==========Main function ProcessSeries==========*)
ProcessSeries[ResolveSunsetMBRepOut_, SeriesNumber_Integer, 
   Evaluate_ : False, opts : OptionsPattern[]] := 
  Module[{EvaluateSeriesOut, GenericSeriesTerm, IndexList = {}, expr, 
    SeriesList = {}, list1D, list2D, finalresult, k, j, 
    verbose = 
     OptionValue[Verbose]},(*Check the presence of EvaluateSeries*)
   If[! MemberQ[Names["Global`*"], "EvaluateSeries"] && ! 
      MemberQ[Names["MBConicHulls`*"], "EvaluateSeries"], 
    Message[ProcessSeries::noeval];
    Return[$Failed];];
   (*Step 1: Evaluation of the generic series*)
   EvaluateSeriesOut = 
    MBConicHulls`EvaluateSeries[ResolveSunsetMBRepOut, {}, 
      SeriesNumber, PrintSeries -> False];
   (*Print["ESO -> ", EvaluateSeriesOut];*)
   GenericSeriesTerm = EvaluateSeriesOut[[2]];
   (*Standardize indices: transform Subscript[_, 
   i_Integer] into Subscript[n, i]*)
   GenericSeriesTerm = 
    GenericSeriesTerm /. Subscript[s_, i_Integer] :> Subscript[n, i];
   (*Print["GST -> ", GenericSeriesTerm];*)(*Construction of IndexList*)
   IndexList = 
    Table[Sort[
      DeleteDuplicates[
        Cases[term, Subscript[n, _Integer], Infinity]]], {term, 
      GenericSeriesTerm}];
   If[verbose, Print["IndexList -> ", IndexList]];
   (*Removal of useless i0*)
   For[j = 1, j <= Length[GenericSeriesTerm], j++, 
    GenericSeriesTerm[[j]] = GenericSeriesTerm[[j]] //. i0 :> 0];
   (*Assumptions for simplification*)$Assumptions = 
    And @@ (Element[#, PositiveIntegers] & /@ Flatten[IndexList]);
   SetOptions[Simplify, 
    TransformationFunctions -> {Automatic, (# /. {a1_^(b1_*
              c1_) :> (a1^b1)^c1} &)}];
   SeriesList = {};
   For[k = 1, k <= Length[GenericSeriesTerm], k++, 
    expr = FactorGammaPositive[GenericSeriesTerm[[k]]];
    If[verbose, Print[k, "_1 -> ", expr]];
    expr = GMultiplication[expr];
    If[verbose, Print[k, "_2 -> ", expr]];
    expr = Reflection[expr];
    If[verbose, Print[k, "_3 -> ", expr]];
    expr = ReplacePochhammer[expr];
    If[verbose, Print[k, "_4 -> ", expr]];
    expr = FormatGamma[expr, k, IndexList];
    If[verbose, Print[k, "_5 -> ", expr]];
    expr = Simplify[expr];
    If[verbose, Print[k, "_6 -> ", expr]];
    expr = {expr[[1]], IdentificationByIndex[expr[[2]], IndexList, k], 
       expr[[3]]};
    If[verbose, Print[k, "_7 -> ", expr]];
    expr[[1]] = ProcessPrefactor[expr[[1]]];
    If[verbose, Print[k, "_8 -> ", expr]];
    expr = FunctionIdentification[expr, IndexList, k];
    If[verbose, Print[k, "_9 -> ", expr]];
    AppendTo[SeriesList, expr];];
   list1D = Select[SeriesList, KeyMemberQ[#, "SeriesArg"] &];
   list2D = Select[SeriesList, KeyMemberQ[#, "SeriesArgs"] &];
   finalresult = Display1[list1D] + Display2[list2D];
   If[Evaluate, Return[InputForm[finalresult]], Return[finalresult]]];

Protect[ProcessSeries, KampeDeFeriet, LauricellaA, LauricellaB, 
  LauricellaC, LauricellaD, HornFunction];

End[]

EndPackage[]
