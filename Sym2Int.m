(* ::Package:: *)

(* ::Input::Initialization:: *)
failedToLoadGroupMath[]:=(
Print[Style["Problem: failed to load the GroupMath package.",{Bold,Darker[Red],FontFamily->"Consolas"}],Style[" Please ensure that ",{GrayLevel[0.5],FontFamily->"Consolas"}],Style[Hyperlink["GroupMath","http://renatofonseca.net/groupmath.php",Appearance->"Palette"],{GrayLevel[0.5],FontFamily->"Consolas"}],Style[" is installed in any of the following paths:\n",{GrayLevel[0.5],FontFamily->"Consolas"}],Grid[{Style[#,{FontFamily->"Consolas",GrayLevel[0.5]}]}&/@$Path,Alignment->{Left,Baseline}]];
Print[Style["Alternatively, you may try first to load GroupMath manually, and then load Sym2Int.",{GrayLevel[0.5],FontFamily->"Consolas"}]];
);
GroupMathAlreadyLoadedQ=Or@@(!StringFreeQ[#,"GroupMath`"]&/@Contexts[]);
If[!GroupMathAlreadyLoadedQ,
Quiet[Check[Get["GroupMath`"],failedToLoadGroupMath[]]];
];


(* ::Input::Initialization:: *)
SMexample="(* ***** Standard Model ***** *)\ngaugeGroup[SM]^={SU3,SU2,U1};\n\nfld1={\"u\",{3,1,2/3},\"R\",\"C\",3};\nfld2={\"d\",{3,1,-1/3},\"R\",\"C\",3};\nfld3={\"Q\",{3,2,1/6},\"L\",\"C\",3};\nfld4={\"e\",{1,1,-1},\"R\",\"C\",3};\nfld5={\"L\",{1,2,-1/2},\"L\",\"C\",3};\nfld6={\"H\",{1,2,1/2},\"S\",\"C\",1};\nfields[SM]^={fld1,fld2,fld3,fld4,fld5,fld6};\n\nGenerateListOfCouplings[SM,MaxOrder\[Rule]4];";

THDMexample="(* ***** Two Higgs doublet model ***** *)\ngaugeGroup[THDM]^={SU3,SU2,U1};\n\nfld1={\"u\",{3,1,2/3},\"R\",\"C\",3};\nfld2={\"d\",{3,1,-1/3},\"R\",\"C\",3};\nfld3={\"Q\",{3,2,1/6},\"L\",\"C\",3};\nfld4={\"e\",{1,1,-1},\"R\",\"C\",3};\nfld5={\"L\",{1,2,-1/2},\"L\",\"C\",3};\nfld6={\"H1\",{1,2,1/2},\"S\",\"C\",1};\nfld7={\"H2\",{1,2,1/2},\"S\",\"C\",1};\nfields[THDM]^={fld1,fld2,fld3,fld4,fld5,fld6,fld7};\n\nGenerateListOfCouplings[THDM];";

SMSeesawIIexample="(* ***** Standard Model plus one complex scalar triplet ***** *)\ngaugeGroup[SMSeesawII]^={SU3,SU2,U1};\n\nfld1={\"u\",{3,1,2/3},\"R\",\"C\",3};\nfld2={\"d\",{3,1,-1/3},\"R\",\"C\",3};\nfld3={\"Q\",{3,2,1/6},\"L\",\"C\",3};\nfld4={\"e\",{1,1,-1},\"R\",\"C\",3};\nfld5={\"L\",{1,2,-1/2},\"L\",\"C\",3};\nfld6={\"H\",{1,2,1/2},\"S\",\"C\",1};\nfld7={\"\[CapitalDelta]\",{1,3,1},\"S\",\"C\",1};\nfields[SMSeesawII]^={fld1,fld2,fld3,fld4,fld5,fld6,fld7};\n\nGenerateListOfCouplings[SMSeesawII,MaxOrder\[Rule]4];";

SVSexample="gaugeGroup[SVS331Model] ^= {SU3, SU3, U1};\n\n\[Psi]l = {\"\[Psi]l\", {1, 3, -1/3}, \"L\", \"C\", 3};\nec = {\"ec\", {1, 1, 1}, \"L\", \"C\", 3};\nQ12 = {\"Q12\", {3, -3, 0}, \"L\", \"C\", 2};\nQ3 = {\"Q3\", {3, 3, 1/3}, \"L\", \"C\", 1};\nuc = {\"uc\", {-3, 1, -2/3}, \"L\", \"C\", 4};\ndc = {\"dc\", {-3, 1, 1/3}, \"L\", \"C\", 5};\n\n\[Phi]1 = {\"\[Phi]1\", {1, 3, 2/3}, \"S\", \"C\", 1};\n\[Phi]23 = {\"\[Phi]23\", {1, 3, -1/3}, \"S\", \"C\", 2}; (* Note that \[Phi]2 and \[Phi]3 are seems as two flavors of \[Phi]23 *)\n\nfields[SVS331Model] ^= {\[Psi]l, ec, Q12, Q3, uc, dc, \[Phi]1, \[Phi]23};\nGenerateListOfCouplings[SVS331Model];";

PPFexample="gaugeGroup[PPF331Model]^={SU3,SU3,U1};\n\n\[Psi]l={\"\[Psi]l\",{1,3,0},\"L\",\"C\",3};\nQ23L={\"Q23L\",{3,-3,-1/3},\"L\",\"C\",2};\nQ1L={\"Q1L\",{3,3,2/3},\"L\",\"C\",1};\nuc={\"uc\",{-3,1,-2/3},\"L\",\"C\",3};\ndc={\"dc\",{-3,1,1/3},\"L\",\"C\",3};\nJ12={\"J12\",{-3,1,4/3},\"L\",\"C\",1};\nJ3={\"J3\",{-3,1,-5/3},\"L\",\"C\",2};\n\n\[Chi]={\"\[Chi]\",{1,3,-1},\"S\",\"C\",1};\n\[Eta]={\"\[Eta]\",{1,3,0},\"S\",\"C\",1};\n\[Rho]={\"\[Rho]\",{1,3,1},\"S\",\"C\",1};\n\nfields[PPF331Model]^={\[Psi]l,Q23L,Q1L,uc,dc,J12,J3,\[Chi],\[Eta],\[Rho]};\n\nGenerateListOfCouplings[PPF331Model];";

SU5example="gaugeGroup[modelSU5]^={SU5};\n\nfld1={\"F\",{-5},\"L\",\"C\",3};\nfld2={\"T\",{10},\"L\",\"C\",3};\nfld3={\"5\",{5},\"S\",\"C\",1};\nfld4={\"24\",{24},\"S\",\"R\",1};\nfields[modelSU5]^={fld1,fld2,fld3,fld4};\n\nGenerateListOfCouplings[modelSU5];";

SMEFTExamplePart1="(* ***** Standard Model Effective Field Theory (SMEFT) up to dimension 10 ***** *)\ngaugeGroup[SM]^={SU3,SU2,U1};\n\nfld1={\"u\",{3,1,2/3},\"R\",\"C\",3};\nfld2={\"d\",{3,1,-1/3},\"R\",\"C\",3};\nfld3={\"Q\",{3,2,1/6},\"L\",\"C\",3};\nfld4={\"e\",{1,1,-1},\"R\",\"C\",3};\nfld5={\"L\",{1,2,-1/2},\"L\",\"C\",3};\nfld6={\"H\",{1,2,1/2},\"S\",\"C\",1};\nfields[SM]^={fld1,fld2,fld3,fld4,fld5,fld6};\n\nsavedResults=GenerateListOfCouplings[SM,MaxOrder\[Rule]10,Verbose\[Rule]\"OnlyStatistics\"];";

SMEFTExamplePart2="(* Convert result into the notation of the datafiles of JHEP 1708 (2017) 016, arXiv:1512.03433 [hep-ph] *)\nConvertSym2IntResult[termAll_]:=Module[{rule,aux,result},\nrule={-10\[Rule]Br,-9\[Rule]Wr,-8\[Rule]Gr,-6\[Rule]Hd,-5\[Rule]Ld,-4\[Rule]e,-3\[Rule]Qd,-2\[Rule]d,-1\[Rule]u,0\[Rule]t,1\[Rule]ud,2\[Rule]dd,3\[Rule]Q,4\[Rule]ed,5\[Rule]L,6\[Rule]H,8\[Rule]Gl,9\[Rule]Wl,10\[Rule]Bl};\n\naux=If[termAll[[4]],{termAll[[2]]},{termAll[[2]],-termAll[[2]]}];\nresult=Expand[termAll[[5]]Total[Times@@@(aux/.rule)]];\n\nReturn[result];\n]\n\nTotal[ConvertSym2IntResult/@savedResults]";

SMEFTExamplePart3="(* Select operators which violate baryon number by +-2 units *)\nBaryonNumber[term_]:=Module[{baryonNumber},\nbaryonNumber=Sign[term[[2]]].(Abs[term[[2]]]/.{0\[Rule]0,1\[Rule]1/3,2\[Rule]1/3,3\[Rule]1/3,4\[Rule]0,5\[Rule]0,6\[Rule]0,8\[Rule]0,9\[Rule]0,10\[Rule]0});\nReturn[baryonNumber];\n]\n\naux=Cases[savedResults,x_/;Abs[BaryonNumber[x]]\[Equal]2\[RuleDelayed]x[[1]]];\nPrintOperatorTable[SM,savedResults[[aux]]]";


Block[{result},
result={};
AppendTo[result,Row[{Style["XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX ",{GrayLevel[0.5]}],Hyperlink[Mouseover[Style["Sym2Int",{GrayLevel[0.5]}],Style["Sym2Int",{Darker[Blue,0.5],Bold}]],"http://renatofonseca.net/sym2int.php"],Style[" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",{GrayLevel[0.5]}]}]];
AppendTo[result,Row[{Style["Version: 2.0; Author: Renato Fonseca.",{GrayLevel[0.5]}]}]];
AppendTo[result,Row[{Style["The Sym2Int program list all possible interactions with a given list of fields. ",{GrayLevel[0.5]}],Style[Hyperlink["This","http://renatofonseca.net/sym2int.php"],{GrayLevel[0.5]}],Style[" webpage explains how to use it.",{GrayLevel[0.5]}]}]];
AppendTo[result,Style["XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",{GrayLevel[0.5]}]];

Print[Row[result,"\n",BaseStyle->(FontFamily->"Consolas")]];
];

Print[Style["Get started quickly by considering some examples (click on one of the models below):",{GrayLevel[0.5],FontFamily->"Consolas"}]];


buttonSM=Button[Style["Standard Model",{Darker[Red],Bold,FontFamily->"Consolas"}],CellPrint[Cell[BoxData[SMexample],"Input",Background->Lighter[Orange,0.9]]],Appearance->"DialogBox"];
buttonTHDM=Button[Style["Two Higgs Doublet Model",{Darker[Red],Bold,FontFamily->"Consolas"}],CellPrint[Cell[BoxData[THDMexample],"Input",Background->Lighter[Orange,0.9]]],Appearance->"DialogBox"];
buttonSMSeesawII=Button[Style["SM+scalar triplet",{Darker[Red],Bold,FontFamily->"Consolas"}],CellPrint[Cell[BoxData[SMSeesawIIexample],"Input",Background->Lighter[Orange,0.9]]],Appearance->"DialogBox"];
buttonPPF=Button[Style["Pisano-Pleitez-Frampton 331 Model",{Darker[Red],Bold,FontFamily->"Consolas"}],CellPrint[Cell[BoxData[PPFexample],"Input",Background->Lighter[Orange,0.9]]],Appearance->"DialogBox"];
buttonSVS=Button[Style["Singer-Valle-Schechter 331 Model",{Darker[Red],Bold,FontFamily->"Consolas"}],CellPrint[Cell[BoxData[SVSexample],"Input",Background->Lighter[Orange,0.9]]],Appearance->"DialogBox"];
buttonSU5=Button[Style["SU(5) Model",{Darker[Red],Bold,FontFamily->"Consolas"}],CellPrint[Cell[BoxData[SU5example],"Input",Background->Lighter[Orange,0.9]]],Appearance->"DialogBox"];
buttonSMEFT=Button[Style["SMEFT",{Darker[Red],Bold,FontFamily->"Consolas"}],CellPrint[Cell[BoxData[SMEFTExamplePart3],"Input",Background->Lighter[Orange,0.9],GeneratedCell->False,CellLabel->"PART3"]];CellPrint[Cell[BoxData[SMEFTExamplePart2],"Input",Background->Lighter[Orange,0.9],GeneratedCell->False,CellLabel->"PART2"]];CellPrint[Cell[BoxData[SMEFTExamplePart1],"Input",Background->Lighter[Orange,0.9],GeneratedCell->False,CellLabel->"PART1"]],Appearance->"DialogBox"];
Print[Grid[{{buttonSM,buttonTHDM,buttonSMSeesawII,buttonPPF,buttonSVS,buttonSU5},{buttonSMEFT}}]];


(* ::Input::Initialization:: *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX GenerateListOfCouplings XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)

(* This funtion generates all terms up to some order allowed by gauge symmetry *)
Options[GenerateListOfCouplings]={Verbose->True,HCTerms->False,CalculateSnSymmetries->True,CalculateInvariants->False,DiscreteSym->{},SymExcludedOps->False,MaxOrder->4,MassDimensions->{},GrassmannCorrection->None,IncludeDerivatives->True,UseEOM->True};
GenerateListOfCouplings[model_,OptionsPattern[]]:=Module[{reps\[UnderBracket]discreteCharges,reps\[UnderBracket]lorentz,order,listOfInteractions,result,operatorTypes,outputData,toPrint,aux,vSpace},

tmp=TimeUsed[];

(* .................................................................................... *)
(* ............................... Initialize variables ............................... *)
(* .................................................................................... *)

(* clear this variable *)
fieldsAUX[model]^={};

(* The discrete charges *)
reps\[UnderBracket]discreteCharges=If[OptionValue[DiscreteSym]==={},ConstantArray[1,Length[fields[model]]],OptionValue[DiscreteSym]];
discSym[model]^=reps\[UnderBracket]discreteCharges;

(* The mass dimensions *)
If[OptionValue[MassDimensions]==={},
reps\[UnderBracket]lorentz=SimpleLorentzInputConversion/@fields[model][[All,3]];
massDs[model]^=If[Mod[Total[reps\[UnderBracket]lorentz[[#]],2],2]==0,1,3/2]&/@Range[Length[fields[model]]];
,
massDs[model]^=OptionValue[MassDimensions];
];

(* Grassmann nature of the fields; True=anti-commutes; False=commutes *)
If[OptionValue[GrassmannCorrection]===None,
grassmannNature[model]^=Mod[Total[#,2],2]==1&/@reps\[UnderBracket]lorentz;
,
grassmannNature[model]^=OptionValue[GrassmannCorrection];
];

(* This is the index reserved to a pure derivative, which later on is added to the list of input fields *)
derivativeIndex[model]^=Length[fields[model]]+1;

order=OptionValue[MaxOrder];

(* .................................................................................... *)
(* ............. Find list of interactions - done by FindAllInteractions .............. *)
(* .................................................................................... *)

(* This step may change: fields[model], fieldsAUX[model] and fieldsORIGINAL[model] *)
listOfInteractions=FindAllInteractions[model,MaxOrder->order,HCTerms->OptionValue[HCTerms],IncludeDerivatives->OptionValue[IncludeDerivatives],UseEOM->OptionValue[UseEOM]];
 (* Print["P4 ",TimeUsed[]-tmp]; *)

(* .................................................................................... *)
(* .............. Analize each interaction - done by AnalizeInteraction ............... *)
(* .................................................................................... *)

(* We would like to add more information to resultRaw. In particular, each entry of result, result[[i]], will be of the form  {resultRaw[[i]],<Sn information>, <invariants>}. The last two entries are only calculated if the user so wishes; otherwise their value is {} *)

(* For now only calculate explicit invariants if CalculateInvariants=True and IncludeDerivatives=False; UseFieldsAUX is doing nothing right now *)

result=AnalizeInteraction[model,#,CalculateSnSymmetries->OptionValue[CalculateSnSymmetries],CalculateInvariants->(OptionValue[CalculateInvariants]&&!OptionValue[IncludeDerivatives]),UseFieldsAUX->True]&/@listOfInteractions;

If[!OptionValue[SymExcludedOps]&&(OptionValue[CalculateSnSymmetries]||OptionValue[CalculateInvariants]),
result=DeleteCases[result,x_/;!x[[4]]];
];
result=result[[All,1;;3]];

(* Print["P5 ",TimeUsed[]-tmp]; *)


(* .................................................................................... *)
(* ..................... Condense the various operator types which .................... *)
(* .............. only differ by the field where derivatives are applied .............. *)
(* .................................................................................... *)

operatorTypes=OperatorsOfEachType[model,result,OptionValue[CalculateSnSymmetries]]//Transpose;

(* Delete those cases where the total number of operators is 0 *)
operatorTypes=DeleteCases[operatorTypes,x_/;x[[5]]===0];

outputData=Table[Prepend[ReorganizeOperatorData[model,operatorTypes[[i]],result[[operatorTypes[[i,2]]]],OptionValue[CalculateSnSymmetries]],i],{i,Length[operatorTypes]}];

(* Print["P6 ",TimeUsed[]-tmp]; *)

(* .................................................................................... *)
(* ...... Last part: make arrangements to produce a final list with all the data ...... *)
(* .................................................................................... *)

Sym2Int\[UnderBracket]Table=Null;
vSpace="";
If[OptionValue[Verbose]=!=False&&OptionValue[Verbose]=!="OnlyStatistics",

If[Length[outputData]<=200,
toPrint=BuildTableToPrint[model,outputData,OptionValue[CalculateSnSymmetries],OptionValue[Verbose]];
Sym2Int\[UnderBracket]OperatorTable=toPrint;
Print[toPrint];
vSpace="\n";
,
Print[Style["WARNING: there are more than 200 types of operators. It is not a good idea to try to print them on Mathematica's front end. Use the raw data output of GenerateListOfCouplings instead.",{Bold,FontFamily->"Consolas",Darker[Red]}]];
vSpace="\n";
];
];

If[OptionValue[Verbose]=!=False,
aux=ModelStatistics[outputData,order];
Sym2Int\[UnderBracket]Statistics=aux[[1]];
Sym2Int\[UnderBracket]StatisticsTable=aux[[2]];
Print[Style[vSpace<>"***************************** Statistics ****************************",{Bold,FontFamily->"Consolas"}]];
Print[aux[[2]]];
];

(* Print["P7 ",TimeUsed[]-tmp]; *)

Return[outputData];
]


(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX FindAllInteractions XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)

(* FindAllInteractions generates list of interactions which are invariant under the gauge, lorentz, and discrete groups. The list is elaborated assuming infinite copies of each field. *)
(* Output of the function is a list {<int1>, <int2>,...} where each intN={<field1>,<field2>,...} gives the list of participating fields, given by their position in the model list (and a -1 sign means conjugation) *)
Options[FindAllInteractions]={HCTerms->False,MaxOrder->4,IncludeDerivatives->False,UseEOM->True};
FindAllInteractions[model_,OptionsPattern[]]:=Module[{nameOfTheFields,reps\[UnderBracket]gaugeIn,reps\[UnderBracket]lorentzIn,reps\[UnderBracket]RI,reps\[UnderBracket]copies,reps\[UnderBracket]gauge,reps\[UnderBracket]discreteCharges,order,reps\[UnderBracket]lorentz,fullGroup,theFullReps,theFullRepsC,tableConjClasses,massDimensions,indicesForMod,U1pos,nonU1pos,U1charges,otherCharges,indicesForModOfNonU1s,maxF,aux,trueFieldIndicesRule,reps,result,rule,ruleForSelfConjFields,includeDerivatives,useEOMQ,uniqueFields,gaugeBosonsIndices,gb,derBounded,derFree,aux2,replacementRule},


(* .................................................................................... *)
(* ............................... Initialize variables ............................... *)
(* .................................................................................... *)

order=OptionValue[MaxOrder];
includeDerivatives=OptionValue[IncludeDerivatives];
{nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,massDimensions,reps\[UnderBracket]discreteCharges}=ExtractModelData[model][[1;;7]];
If[includeDerivatives,
AppendTo[nameOfTheFields,"D"];
AppendTo[reps\[UnderBracket]gauge,SimpleRepInputConversion[gaugeGroup[model],If[#===U1,0,1]&/@gaugeGroup[model]]];
AppendTo[reps\[UnderBracket]lorentz,SimpleLorentzInputConversion["V"]];
AppendTo[reps\[UnderBracket]RI,"R"];
AppendTo[reps\[UnderBracket]copies,1];
AppendTo[massDimensions,1];
AppendTo[reps\[UnderBracket]discreteCharges,1];
];

(* Print["P1 ",TimeUsed[]-tmp]; *)

(* .................................................................................... *)
(* .......... Prepare things for recursive function TryAllFieldCombinations ........... *)
(* .................................................................................... *)

fullGroup=Join[gaugeGroup[model],{SU2,SU2}];

theFullReps=Flatten[#,1]&/@Transpose[{reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz}];
theFullRepsC=Flatten[#,1]&/@Transpose[{ConjugateIrrep[gaugeGroup[model],#]&/@reps\[UnderBracket]gauge,Reverse/@reps\[UnderBracket]lorentz}];

tableConjClasses=Flatten/@Join[Table[Table[ConjugacyClass[fullGroup[[i]],theFullReps[[repI,i]]],{i,Length[fullGroup]}],{repI,Length[theFullReps]}],Table[Table[ConjugacyClass[fullGroup[[i]],theFullRepsC[[repI,i]]],{i,Length[fullGroup]}],{repI,Flatten[Position[reps\[UnderBracket]RI,"C"]]}]];
PrependTo[tableConjClasses,0tableConjClasses[[1]]];

(* NOTE: massDimensions is a special list at this point of the code; it includes an entry for the complex-conjugated fields and also for a dummny 'no field' *)
massDimensions=Join[massDimensions,massDimensions[[Flatten[Position[reps\[UnderBracket]RI,"C"]]]]];
PrependTo[massDimensions,1]; (* index 1 corresponds to no field; ie it is used to get the operators with less than MaxOrder *)

indicesForMod=Flatten[ConjugacyClassGroupModIndices/@fullGroup]; (* -1 means it is a U1 charge *)
U1pos=Flatten[Position[indicesForMod,-1]];
nonU1pos=Complement[Range[Length[indicesForMod]],U1pos];
U1charges=tableConjClasses[[All,U1pos]];
otherCharges=tableConjClasses[[All,nonU1pos]];
indicesForModOfNonU1s=indicesForMod[[nonU1pos]];

maxF=Length[massDimensions];

(* .................................................................................... *)
(* .................................................................................... *)
(* .................................................................................... *)
(* ............ Use conjugacy classes and U1 charges to get a reduced list ............ *)
(* .............. of POTENTIALLY gauge and lorentz invariant interactions ............. *)
(* .................................................................................... *)
(* .................................................................................... *)
(* .................................................................................... *)

(* Print["P2 ",TimeUsed[]-tmp]; *)

ClearAll[massDData];
ClearAll[fieldsWhichCanBeAdded];
Do[massDData[i,j]=i+massDimensions[[j]];,{i,0,order,1/2},{j,Length[massDimensions]}];
Do[fieldsWhichCanBeAdded[i,minF]=DeleteCases[Flatten[Position[massDimensions,x_/;x<=i]],y_/;y<minF],{i,1/2,order,1/2},{minF,1,Length[massDimensions]}];

TryAllFieldCombinations[initCombination_,massD_,minF_]:=Block[{tmp},
If[massD==order,
aux=Total[U1charges[[initCombination]]];

If[(aux==0aux),
aux=Total[otherCharges[[initCombination]]];

If[And@@Divisible[aux,indicesForModOfNonU1s],Sow[initCombination]];
];
,
TryAllFieldCombinations[Append[initCombination,#],massDData[massD,#],#]&/@fieldsWhichCanBeAdded[order-massD,minF];
];
];
aux=Drop[Reap[TryAllFieldCombinations[{},0,1]][[2,1]],1];

(* Print["P3 ",TimeUsed[]-tmp]; *)

(* .................................................................................... *)
(* ................ Some clean up and change of notation is necessary ................. *)
(* .................................................................................... *)

trueFieldIndicesRule=Join[Range[Length[theFullReps]],-Flatten[Position[reps\[UnderBracket]RI,"C"]]];
trueFieldIndicesRule=MapThread[Rule,{1+Range[Length[trueFieldIndicesRule]],trueFieldIndicesRule}];
aux=DeleteCases[aux,1,{2}]/.trueFieldIndicesRule;

(* .................................................................................... *)
(* .... If derivatives are included, then the model to be evaluated must be changed ... *)
(* .................................................................................... *)

If[includeDerivatives,
aux=AddGaugeBosons[model,aux,massDimensions,order,1+Length[fields[model]]];
aux=ApplyDerivativesToFields[aux,1+Length[fields[model]]];

(* Print["P3-b ",TimeUsed[]-tmp]; *)

useEOMQ=True;

(* modify model reps *)
uniqueFields=Flatten[aux,1]//DeleteDuplicates;
uniqueFields=DeleteCases[DeleteDuplicates[Abs[uniqueFields]],x_/;x[[1]]==1+Length[fields[model]]];
uniqueFields=SortBy[uniqueFields,(2#[[2]]Length[uniqueFields]+#[[1]])&];

fieldsAUX[model]^=ModifyModelFieldsToIncludeDerivatives[model,uniqueFields,useEOMQ,1+Length[fields[model]],massDimensions[[2;;-1]],reps\[UnderBracket]discreteCharges];


{nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,massDimensions,reps\[UnderBracket]discreteCharges}=ExtractModelData[model][[1;;7]];
(* PrependTo[massDimensions,1]; (* because massDimensions[[1]] used to be a placeholder for a 'no field' *) *)


(* Change notation of the list of interactions obtained earlier *);
aux2=GatherBy[fieldsAUX[model],#[[2]]&];
replacementRule=Join[Table[aux2[[i,1,2]]->aux2[[i,All,1]],{i,Length[aux2]}],Table[{-1,1}aux2[[i,1,2]]->-aux2[[i,All,1]],{i,Length[aux2]}]];

aux=Flatten[Tuples/@(aux/.Dispatch[replacementRule]),1];
,
(* Even if no derivatives are to be considered, it is useful to populate fieldsAUX[model] as in the derivative case so that later on the same code can be used in all cases *)

aux2=Table[{fldI,{fldI,0},fields[model][[fldI,1]],SimpleRepInputConversion[gaugeGroup[model],fields[model][[fldI,2]]],SimpleLorentzInputConversion[fields[model][[fldI,3]]],fields[model][[fldI,4]],fields[model][[fldI,5]],massDimensions[[1+fldI]],reps\[UnderBracket]discreteCharges[[fldI]],grassmannNature[model][[fldI]]},{fldI,Length[fields[model]]}];

fieldsAUX[model]^=AppendTo[aux2,{Length[fields[model]]+1,{Length[fields[model]]+1,0},"D",Null,Null,Null,Null,Null,Null,Null}];
];

(* .................................................................................... *)
(* ............. Now get the correct list of gauge-invariant interactions ............. *)
(* .................................................................................... *)

reps=MapThread[List,{Range[Length[reps\[UnderBracket]gauge]],reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]discreteCharges,reps\[UnderBracket]RI,reps\[UnderBracket]copies}];
rule=Join[Table[i->reps[[i]],{i,Length[reps\[UnderBracket]gauge]}],Table[-i->ConjugateRepSpecial[gaugeGroup[model],reps[[i]]],{i,1,Length[reps\[UnderBracket]gauge]}]];

result=Cases[aux,x_/;ContainsSingletQ[gaugeGroup[model],x/.rule]];

(* .................................................................................... *)
(* ............................. Maybe cut down h.c. terms ............................ *)
(* .................................................................................... *)

If[!OptionValue[HCTerms],
ruleForSelfConjFields=MapThread[Rule,{-Flatten[Position[reps\[UnderBracket]RI,"R"]],Flatten[Position[reps\[UnderBracket]RI,"R"]]}];

(* Chosen terms are those with less conjugations: SortBy[#,(Count[#,x_/;x<0]&)]& *)
result=DeleteDuplicates[SortBy[#,(Count[#,x_/;x<0]&)]&/@MapThread[List,{Sort/@result,Sort/@((-result)/.ruleForSelfConjFields)}]][[All,1]];
];

(* .................................................................................... *)
(* .... Delete dim 4 or smaller interactions with derivatives and/or gauge bosons ..... *)
(* .................................................................................... *)

If[includeDerivatives,
gaugeBosonsIndices=Length[fields[model]]+1+Range[Length[gaugeGroup[model]]];
gb=If[MemberQ[gaugeBosonsIndices,#[[2,1]]],1,0]&/@fieldsAUX[model];
derBounded=fieldsAUX[model][[All,2,2]];
derFree=UnitVector[Length[fieldsAUX[model]],derivativeIndex[model]];
result=DeleteCases[result,x_/;Total[massDimensions[[Abs[x]]]]<=4&&(Total[gb[[Abs[x]]]]+Total[derBounded[[Abs[x]]]]+Total[derFree[[Abs[x]]]])>0,{1}];
];

(* .................................................................................... *)
(* .................................. Sort results .................................... *)
(* .................................................................................... *)

(* Sort results by operator dimension *)
result=SortBy[result,Total[massDimensions[[Abs[#]]]]&];

(* Now we would like to sort the fields in each interaction: this is done by Abs[<pos of the field given by user>] as the first criteria, and then conjugated fields last (second criteria) *)
result=SortBy[#,({Abs[#],-#}&)]&/@result;

Return[result];
]


(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX AnalizeInteraction XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)


(* AnalizeInteraction considers more information about each interaction: (A) the permutation symmetry, (B) the explicit invariants, and  (C) if the operator exists at all for the given number of flavors [based on information (A)] *)
(* The interaction should have the form {<field1>,<field2>,...}, where each field is identified with its position in the model list (and a -1 sign means conjugation) *)

Options[AnalizeInteraction]={CalculateSnSymmetries->True,CalculateInvariants->False,UseFieldsAUX->False};
AnalizeInteraction[model_,fieldPositions_,OptionsPattern[]]:=Module[{nameOfTheFields,reps\[UnderBracket]gaugeIn,reps\[UnderBracket]RI,reps\[UnderBracket]copies,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,interaction,result,survivingInvariants,nameHeadsToUse,distinguishReps,gReps,lorentzRep,rep,aux,aux2,isFermionQ,needToConsiderGrassmannNature,fermionsPos,nFlavs,nFlavsMin,validInvsPos,reps,grassmannCorrection},

(* .................................................................................... *)
(* .................................. Initializations ................................. *)
(* .................................................................................... *)

{nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,grassmannCorrection}=ExtractModelData[model][[{1,2,3,4,5,8}]];

(* interaction has the format {# of the field,gauge rep (maybe already conjugated),Lorentz rep (maybe already conjugated),discreteSym, is it conjugated?0=no,1=yes, # copies/flavours of the field} *)
interaction=Table[{Abs[el],If[el>0,reps\[UnderBracket]gauge[[Abs[el]]],ConjugateIrrep[gaugeGroup[model],reps\[UnderBracket]gauge[[Abs[el]]]]],If[el>0,reps\[UnderBracket]lorentz[[Abs[el]]],Reverse[reps\[UnderBracket]lorentz[[Abs[el]]]]],"DS",If[el>0,0,1],reps\[UnderBracket]copies[[Abs[el]]]},{el,fieldPositions}];
result={interaction,Null,Null,Null};


(* .................................................................................... *)
(* .......... Calculate the invariants (comes before permutation symmetry) ............ *)
(* .................................................................................... *)

If[OptionValue[CalculateInvariants],
Print[fieldPositions,"  ",reps\[UnderBracket]lorentz];
survivingInvariants=LorentzInvariants[gaugeGroup[model],reps\[UnderBracket]gauge[[interaction[[All,1]]]],reps\[UnderBracket]lorentz[[interaction[[All,1]]]],(interaction[[All,5]]/.{0->False,1->True}),DistinguishFields->(interaction[[All,1]](interaction[[All,5]]/.{0->-1,1->1}))];
(* Note that at this point, these are not really just the surviving invariants, since all invariants are being kept. Below, at the end of PART4, this is fixed. *)

nameHeadsToUse=Table[nameOfTheFields[[interaction[[i,1]]]]<>If[interaction[[i,5]]==0,"","C"],{i,Length[interaction]}];
survivingInvariants[[2,All,All,0]]=Table[ConstantArray[nameHeadsToUse[[elI]],Length[survivingInvariants[[2,elI]]]],{elI,Length[survivingInvariants[[2]]]}];

result[[3]]=survivingInvariants;
];


(* .................................................................................... *)
(* ................ Look at the permutation symmetry of the invariants ................ *)
(* .................................................................................... *)

If[OptionValue[CalculateSnSymmetries],
distinguishReps=interaction[[All,1]](1-2interaction[[All,5]]);
(* gReps=If[#[[5]]\[Equal]0,#[[2]],ConjugateIrrep[gaugeGroup,#[[2]]]]&/@interaction; *)
gReps=#[[2]]&/@interaction;
(* lorentzRep=If[#[[5]]\[Equal]0,#[[3]],Reverse[#[[3]]]]&/@interaction; *)
lorentzRep=#[[3]]&/@interaction;
reps=Flatten[#,1]&/@MapThread[List,{gReps,lorentzRep}];

If[OptionValue[CalculateInvariants],
aux=$GroupMath\[UnderBracket]Invariants\[UnderBracket]Symmetries;
,
aux=PermutationSymmetryOfInvariants[Join[gaugeGroup[model],{SU2,SU2}],reps,DistinguishFields->distinguishReps];
];

(* Correct Grassmann nature of fermions *)
isFermionQ=grassmannCorrection[[Abs[fieldPositions]]];

needToConsiderGrassmannNature=isFermionQ[[#[[1]]]]&&Length[#]>1&/@aux[[1]];
fermionsPos=Flatten[Position[needToConsiderGrassmannNature,True]];

aux[[2,All,1,fermionsPos]]=Map[SnRepTimesTotallyAntisymmetric,aux[[2,All,1,fermionsPos]],{2}];

nFlavs=interaction[[All,6]];
nFlavs=nFlavs[[#]]&/@aux[[1,All,1]];
(* nFlavsMin=If[NumericQ[#],#,1]&/@nFlavs; *)
aux[[2]]=Table[{aux[[2,i,1]],Times@@MapThread[HookContentFormula,{aux[[2,i,1]],nFlavs}],aux[[2,i,2]]},{i,Length[aux[[2]]]}];

(* Delete those cases which cannot be realized with the given number of flavors *)
aux[[2]]=DeleteCases[aux[[2]],x_/;x[[2]]===0];

If[OptionValue[CalculateInvariants],
aux2=Table[{aux[[2,i,3]]Times@@(SnIrrepDim/@aux[[2,i,1]]),aux[[2,i,2]],aux[[2,i,3]]},{i,Length[aux[[2]]]}];
aux2=Range@@@Transpose[{Accumulate[aux2[[All,1]]]-aux2[[All,1]]+1,Accumulate[aux2[[All,1]]]}];
validInvsPos=Flatten[aux2[[Flatten[Position[aux[[2]],x_/;x[[2]]=!=0,{1},Heads->False]]]]];
result[[3]]={survivingInvariants[[1,validInvsPos]],survivingInvariants[[2]]};
];

result[[2]]=aux;
];

(* Can the interaction be realized with the given number of fields? *)
result[[4]]=!(Length[result[[2]]]==2&&result[[2,2]]==={});

(* simplify 1st entry of the output, which identifies the interaction (should be the same as result[[1,All,1]] (1-2result[[1,All,5]])) *)
result[[1]]=fieldPositions;
Return[result];

];

(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX DetailedInteractionInformation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)

(* model=SM for example; interaction should be contain the field positions, with - signs eventually for conjugations example: {-5,-5,4,4,1,1,-3,-3,-2,-2,-2,-2} *)

DetailedInteractionInformation[model_,interaction_]:=Module[{aux,aux1,aux2,nCopies,fieldNames,inputMod,whichFields,fermionQ,productsOfFieldsAux,productsOfFields,result,dataToPresent,posTransition,divsV,divsH},

aux=Tally[interaction][[All,1]];
fieldNames=If[#>0,fields[model][[#,1]]<>"[R]",fields[model][[-#,1]]<>"[C]"]&/@aux;

inputMod=Tally[interaction];
inputMod[[All,1]]=Join@@@PositionToReps[model,inputMod[[All,1]]][[All,2;;3]];

productsOfFields={};
Do[
whichFields=el;
fermionQ=(Mod[whichFields[[1,4,1]]+whichFields[[1,5,1]],2]==1);
aux=PermutationSymmetryOfTensorProductParts[Join[gaugeGroup[model],{SU2,SU2}],ConstantArray[whichFields[[1]],whichFields[[2]]]];
aux=aux[[2]];

If[fermionQ,
aux[[All,1,2]]=IntegerPartitions[Total[#[[1]]]][[Position[DecomposeSnProduct[{#[[1]],ConstantArray[1,Total[#[[1]]]]}],1][[1,1]]]]&/@aux[[All,1,2]],
(* Just remove a parentesis *)
aux[[All,1,2]]=#[[1]]&/@aux[[All,1,2]]; 
];
(* productsOfFieldsAux=If[!fermionQ,Cases[aux[[2,All,1]],x_/;x[[2,1]]\[Equal]{whichFields[[2]]}\[RuleDelayed]x[[1]]],Cases[aux[[2,All,1]],x_/;x[[2,1]]\[Equal]ConstantArray[1,whichFields[[2]]]\[RuleDelayed]x[[1]]]]; *)
AppendTo[productsOfFields,aux];

,{el,inputMod}];

result={Cases[ReduceRepProduct[Join[gaugeGroup[model],{SU2,SU2}],#[[All,1,1]]],x_/;x[[1]]==0x[[1]]:>x[[2]]],#[[All,1,1]],#[[All,1,2]],Times@@#[[All,2]]}&/@Tuples[productsOfFields];
result={#[[2]],#[[3]],#[[1,1]]#[[4]]}&/@DeleteCases[result,x_/;x[[1]]==={}];
result=SortBy[result,#[[2]]&];

(* --- Get things ready to present the data ---- *)

If[Length[result]==0,Return[result]];

dataToPresent=Prepend[result,{Style[fieldNames,{Bold,Darker[Red]}],Style[Tally[interaction][[All,2]],{Bold,Darker[Red]}],Style["Number of invariants",{Bold,Darker[Red]}]}];
posTransition=Position[result[[All,2]],#][[1,1]]&/@(result[[All,2]]//DeleteDuplicates);

divsV={1->{Black,Thick},2->Lighter[Gray,0.8],3->Lighter[Gray,0.8],4->Lighter[Gray,0.8],5->{Black,Thick}};
divsH=Join[{1->{Black,Thick},2->Black,-1->{Black,Thick}},Table[i->{Black,Dashed},{i,Drop[posTransition,1]+1}],Table[i->Lighter[Gray,0.8],{i,Complement[Range[Length[posTransition]],posTransition]+1}]];

nCopies=fields[model][[Abs[Tally[interaction][[All,1]]],5]];
aux1=Times@@MapThread[HookContentFormula,{#,nCopies}]&/@Tally[result[[All,2]]][[All,1]];
aux2=ConstantArray[SpanFromAbove,#-1]&/@Tally[result[[All,2]]][[All,2]];
aux=Prepend[Flatten[MapThread[Prepend,{aux2,If[#==0,Item[#,Background->Lighter[Red,0.85]],Item[#,Background->Lighter[Green,0.85]]]&/@aux1}]],"# Parameters"];
Print[Grid[Transpose[Append[Transpose[dataToPresent],aux]],Frame->All,FrameStyle->LightGray,Spacings->{1,1},Dividers->{divsV,divsH},Background->{None,{None,{Lighter[Gray,0.95],None}}},Alignment->{Center,Center}]];

(* Add number of parameters information *)
result=Transpose[Append[Transpose[result],Times@@MapThread[HookContentFormula,{#,nCopies}]&/@result[[All,2]]]];
Return[result];
]

(* auxiliar function of DetailedInteractionInformation *)
PositionToReps[model_,fieldPositions_]:=Module[{nameOfTheFields,reps\[UnderBracket]gaugeIn,reps\[UnderBracket]lorentzIn,reps\[UnderBracket]RI,discreteSym,reps\[UnderBracket]copies,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,interaction},
{nameOfTheFields,reps\[UnderBracket]gaugeIn,reps\[UnderBracket]lorentzIn,reps\[UnderBracket]RI,reps\[UnderBracket]copies}=Transpose[fields[model]];

(* Just converts the simplified input into Dynkin coefficients *)
reps\[UnderBracket]gauge=SimpleRepInputConversion[gaugeGroup[model],#]&/@reps\[UnderBracket]gaugeIn;
reps\[UnderBracket]lorentz=SimpleLorentzInputConversion/@reps\[UnderBracket]lorentzIn;

(* interaction has the format{# of the field,gauge rep (maybe already conjugated),Lorentz rep (maybe already conjugated),discreteSym, is it conjugated?0=no,1=yes, # copies/flavours of the field} *)
interaction=Table[{Abs[el],If[el>0,reps\[UnderBracket]gauge[[Abs[el]]],ConjugateIrrep[gaugeGroup[model],reps\[UnderBracket]gauge[[Abs[el]]]]],If[el>0,reps\[UnderBracket]lorentz[[Abs[el]]],Reverse[reps\[UnderBracket]lorentz[[Abs[el]]]]],"DS",If[el>0,0,1],reps\[UnderBracket]copies[[Abs[el]]]},{el,fieldPositions}];
Return[interaction];
]

(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Some auxiliar functions XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)


(* Use it to get {nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,massDimensions,reps\[UnderBracket]discreteCharges,reps\[UnderBracket]grassmannNature} of a model *)
ExtractModelData[model_]:=Module[{nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,massDimensions,reps\[UnderBracket]discreteCharges,data,reps\[UnderBracket]grassmannNature},
If[!ValueQ[fieldsAUX[model]]||(ValueQ[fieldsAUX[model]]&&fieldsAUX[model]==={}),
{nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies}=Transpose[fields[model]];
massDimensions=massDs[model];
reps\[UnderBracket]discreteCharges=discSym[model];
reps\[UnderBracket]grassmannNature=grassmannNature[model];

(* Just converts the simplified input into Dynkin coefficients *)
reps\[UnderBracket]gauge=SimpleRepInputConversion[gaugeGroup[model],#]&/@reps\[UnderBracket]gauge;
reps\[UnderBracket]lorentz=SimpleLorentzInputConversion/@reps\[UnderBracket]lorentz;
,
{nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,massDimensions,reps\[UnderBracket]discreteCharges,reps\[UnderBracket]grassmannNature}=Transpose[fieldsAUX[model]][[{3,4,5,6,7,8,9,10}]];
];

data={nameOfTheFields,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,massDimensions,reps\[UnderBracket]discreteCharges,reps\[UnderBracket]grassmannNature};
Return[data];
]

(* Converts dim -dim to the Dynkin coefficients of the rep, anti-rep with that dimension (and no primes). If the input, for a given factor group, are already the Dynkin coefficients, that is fine too (no processing is done). *)
SimpleRepInputConversion[group_,rep_]:=SimpleRepInputConversion[group,rep]=Module[{dim,reps,result},
If[!IsSimpleGroupQ[group],Return[SimpleRepInputConversion@@@Transpose[{group,rep}]]];
If[group==U1,Return[rep]];
If[Head[rep]===List,Return[rep]];

dim=Abs[rep];
reps=RepsUpToDimNNoConjugates[group,dim];
reps=DeleteCases[reps,x_/;DimR[group,x]!=dim];
reps=Sort[reps,OrderedQ[{Join[{DimR[group,#1],RepresentationIndex[group,#1]},ConjugacyClass[group,#1],-#1],Join[{DimR[group,#2],RepresentationIndex[group,#2]},ConjugacyClass[group,#2],-#2]}]&];

result=If[rep>0,reps[[1]],ConjugateIrrep[group,reps[[1]]]];
Return[result];
]

(* One can give a Lorentz rep by name ("S","L","R","V"), {jL,jR}, or with Dynkin coeffiecients *)
SimpleLorentzInputConversion[rep_]:=SimpleLorentzInputConversion[rep]=Module[{},
If[Head[rep]===String,
Return[rep/.{"R"->{{0},{1}},"L"->{{1},{0}},"S"->{{0},{0}},"V"->{{1},{1}}}];
,
If[Depth[rep]==2,Return[2{{rep[[1]]},{rep[[2]]}}];,Return[rep];];
];
];

(* Auxiliar function to LorentzInvariants. Try nToMod[{2,3,4}] *)
nToMod[bIs_]:=Module[{ns,ii,res},
ns=Range[1,(Times@@bIs)]-1;
res={};
Do[
ii=Mod[ns,el];
ns=(ns-ii)/el;
PrependTo[res,ii];
,{el,Reverse[bIs]}];
Return[Transpose[res]+1];
]

(* Converts a Sn rep given by a partition of n into a formated name [specifically for use with GenerateListOfCouplings] *)
SnRepNameTemp[\[Lambda]_]:=Which[Length[\[Lambda]]==1,Style["S",{Bold,Orange}],Total[\[Lambda]]==Length[\[Lambda]],Style["A",{Bold,Blue}],True,Style[\[Lambda],{Bold,Purple}]]

SnRepNameTemp2[\[Lambda]_]:=Which[Length[\[Lambda]]==1,Style[\[Lambda],{Bold,Orange}],Total[\[Lambda]]==Length[\[Lambda]],Style[\[Lambda],{Bold,Blue}],True,Style[\[Lambda],{Bold,Purple}]]

ConjugateRepSpecial[gaugeGroup_,rep_]:={rep[[1]],ConjugateIrrep[gaugeGroup,rep[[2]]],Reverse[rep[[3]]],Conjugate[rep[[4]]],rep[[5]],rep[[6]]};

ContainsSingletQ[gaugeGroup_,reps_]:=Module[{testLorentz,testGauge,testDiscrete},
testDiscrete=((Times@@reps[[All,4]])==1+0(Times@@reps[[All,4]]));
If[!testDiscrete,Return[False,Module]];
Do[testLorentz=MemberQ[ReduceRepProduct[SU2,Sort[reps[[All,3,i]]]],x_/;x[[1]]0==x[[1]]];If[!testLorentz,Return[False,Module]];,{i,2}];
testGauge=True; (* this line is important for the case when gaugeGroup={} *)
Do[testGauge=MemberQ[ReduceRepProduct[gaugeGroup[[i]],Sort[reps[[All,2,i]]]],x_/;x[[1]]0==x[[1]]];If[!testGauge,Return[False,Module]];,{i,Length[gaugeGroup]}];
Return[testDiscrete&&testLorentz&&testGauge];
];

(* Return the Sn rep which equals the product of \[Lambda] times the totally antisymmetric representation of Sn *)
SnRepTimesTotallyAntisymmetric[\[Lambda]_]:=Module[{n,antiS,result},
n=Total[\[Lambda]];
antiS=ConstantArray[1,n];
result=IntegerPartitions[n][[Position[DecomposeSnProduct[{\[Lambda],antiS}],1][[1,1]]]];
Return[result];
]

(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Some auxiliar functions for handling derivatives and gauge bosons XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)

(* Returns list of elements corresponding to the model fields, with the format {<new field index>,{<original field index>,<n derivatives>},<name>,<gauge q.n.>,<lorentz q.n.>,<"R" or "C">,<# flavors>,<mass dimension>,<discrete charges>,<is Grassmann?>} *)

ModifyModelFieldsToIncludeDerivatives[model_,listOfFieldsAppearing_,useEOMQ_,derIdx_,massDimensions_,reps\[UnderBracket]discreteCharges_]:=Module[{n,counter,result,fullName,gaugeQN,lorentzQN,realQ,nF,nameOfTheFields,reps\[UnderBracket]gaugeIn,reps\[UnderBracket]gauge,reps\[UnderBracket]lorentzIn,reps\[UnderBracket]lorentz,reps\[UnderBracket]RI,reps\[UnderBracket]copies,rootName,rootLorentzQN,FmunuQN,massD,discreteQN,grassmann,derivative},
If[useEOMQ,
n=Length[fields[model]];

{nameOfTheFields,reps\[UnderBracket]gaugeIn,reps\[UnderBracket]lorentzIn,reps\[UnderBracket]RI,reps\[UnderBracket]copies}=Transpose[fields[model]];
(* Just converts the simplified input into Dynkin coefficients *)
reps\[UnderBracket]gauge=SimpleRepInputConversion[gaugeGroup[model],#]&/@reps\[UnderBracket]gaugeIn;

reps\[UnderBracket]lorentz=SimpleLorentzInputConversion/@reps\[UnderBracket]lorentzIn;

FmunuQN=Table[Adjoint[gaugeGroup[model]]UnitVector[Length[gaugeGroup[model]],gI],{gI,Length[gaugeGroup[model]]}];
counter=0;
result=Reap[Do[

rootName=If[el[[1]]<=n,nameOfTheFields[[el[[1]]]],"F"<>ToString[el[[1]]-n-1]];

gaugeQN=If[el[[1]]<=n,reps\[UnderBracket]gauge[[el[[1]]]],FmunuQN[[el[[1]]-n-1]]];
rootLorentzQN=If[el[[1]]<=n,reps\[UnderBracket]lorentz[[el[[1]]]],{{2},{0}}];
realQ=If[el[[1]]<=n,reps\[UnderBracket]RI[[el[[1]]]],"C"];
nF=If[el[[1]]<=n,reps\[UnderBracket]copies[[el[[1]]]],1];
massD=If[el[[1]]<=n,massDimensions[[el[[1]]]],2]+el[[2]];
discreteQN=If[el[[1]]<=n,reps\[UnderBracket]discreteCharges[[el[[1]]]],1];
grassmann=If[el[[1]]<=n,grassmannNature[model][[el[[1]]]],False];

If[el[[2]]===0,
fullName=rootName;
lorentzQN=rootLorentzQN;
,
fullName=\!\(\*
TagBox[
StyleBox["\"\<D\>\"",
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)<>ToString[el[[2]]]<>"("<>rootName<>")";
lorentzQN=rootLorentzQN+el[[2]];
];


If[counter==n,
 (* jump one step to give space for the derivative *)
counter++;
derivative={counter,{derIdx,0},"D",SimpleRepInputConversion[gaugeGroup[model],If[#===U1,0,1]&/@gaugeGroup[model]],SimpleLorentzInputConversion["V"],"R",1,1,1,True};
Sow[derivative];
];
counter++;

Sow[{counter,el,fullName,gaugeQN,lorentzQN,realQ,nF,massD,discreteQN,grassmann}];
,{el,listOfFieldsAppearing}]][[2,1]];
];
Return[result];
]


AddGaugeBosons[model_,initListOfInteractionsIn_,massDimensions_,maxOrder_,derIndex_]:=Module[{initListOfInteractions,startFPos,FmunuPos,aux,aux2,FmunuListOfAdditions,massDimsOfTerms,result},
initListOfInteractions=Prepend[initListOfInteractionsIn,{}];

startFPos=Length[fields[model]]+2;
FmunuPos=Range[startFPos,startFPos-1+Length[gaugeGroup[model]]];
FmunuPos=Sort[Join[-FmunuPos,FmunuPos]];

aux=Table[DeleteDuplicates[Sort/@Tuples[FmunuPos,i]],{i,0,maxOrder/2}];
FmunuListOfAdditions=Table[Flatten[aux[[1;;i]],1],{i,1,Length[aux]}];


massDimsOfTerms=Total/@(massDimensions[[#]]&/@(1+Abs[initListOfInteractions]));
aux=Floor[(maxOrder-massDimsOfTerms)/2];
aux2=FmunuListOfAdditions[[aux+1]]; (* combinations of Fmunu which can be added to each term *)

result=Flatten[Table[Join[initListOfInteractions[[i]],el],{i,Length[initListOfInteractions]},{el,aux2[[i]]}],1];
result=Drop[result,1];
result=DeleteCases[result,x_/;Length[x]==Count[x,derIndex]];
Return[result];
]

ApplyDerivativesToFields[initListOfInteractions_,derIndex_]:=Module[{fieldList1,result,nDers,aux,aux2,aux3,groupPartitions},
fieldList1=DeleteCases[initListOfInteractions,x_/;Length[x]==Count[x,derIndex]];
result=Reap[Do[
nDers=Count[initInteraction,derIndex];
aux=Tally[DeleteCases[initInteraction,derIndex]];
Do[
groupPartitions=Flatten[Permutations/@IntegerPartitionsMOD[nDi,Length[aux]],1];
aux2=Flatten[Table[Tuples[IntegerPartitionsMOD[#1,#2]&@@@Transpose[{gp,aux[[All,2]]}]],{gp,groupPartitions}],1];
aux2=Transpose[aux2];
aux3=Flatten[#,1]&/@Transpose[Table[Thread[List[ConstantArray[aux[[i,1]],aux[[i,2]]],#]]&/@aux2[[i]],{i,Length[aux]}]];
(* If[nDi\[NotEqual]nDers,aux3=Prepend[#,{derIndex,nDers-nDi}]&/@aux3]; *)
If[nDi!=nDers,aux3=Join[ConstantArray[{derIndex,0},nDers-nDi],#]&/@aux3];
Sow[aux3];
,{nDi,0,nDers}];
,{initInteraction,fieldList1}]][[2]];
If[Length[result]>0,result=Flatten[result[[1]],1]];
Return[result];
]

(* Used by ApplyDerivativesToFields *)
IntegerPartitionsMOD[n_,parts_]:=PadRight[#,parts]&/@IntegerPartitions[n,parts]

(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX LorentzInvariants XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)

Options[LorentzInvariants]={DistinguishFields->False};
LorentzInvariants[gaugeG_,gaugeReps_,lorentzReps_,conjs_,OptionsPattern[]]:=Module[{fullGroup,lorentzRepsDI,pConj,fullReps,pConjs,indicesOriginal,indicesConjugated,indicesConjugatedMod,repOriginal,rep,invariants,invariantsAr,rule},
fullGroup=Join[gaugeG,{SU2,SU2}];
lorentzRepsDI=SimpleLorentzInputConversion/@lorentzReps;
pConj=Flatten[Position[conjs,True]];
If[pConj!={},
lorentzRepsDI[[pConj]]=Reverse/@lorentzRepsDI[[pConj]];
];

fullReps=Join@@@Transpose[{gaugeReps,lorentzRepsDI}];

invariants=Invariants[fullGroup,fullReps,Conjugations->conjs,TensorForm->True,DistinguishFields->OptionValue[DistinguishFields]];

(* There is still a problem to be solved at this stage. Conjugation switches the Lorentz L < - > R groups, so we get the result in invariants. But then the component ordereing needs to be fixed. For example (L,R)={{2},{1}} gets treated as having components {11,12,13,21,22,23} in invariants, when in reality we would like the ordering {11,12,21,22,31,32}. This issue is fixed in the following code. *)
pConjs=Flatten[Position[conjs,True]];
If[pConjs!={},
invariantsAr=ArrayRules[invariants[[1]]][[1;;-2]];
];
Do[
rep=fullReps[[pos]];
repOriginal=rep;repOriginal[[-2]]=rep[[-1]];repOriginal[[-1]]=rep[[-2]];

indicesOriginal=nToMod[DimR[fullGroup,repOriginal]];
indicesConjugated=nToMod[DimR[fullGroup,rep]];
indicesConjugatedMod=indicesConjugated;indicesConjugatedMod[[All,-1]]=indicesConjugated[[All,-2]];indicesConjugatedMod[[All,-2]]=indicesConjugated[[All,-1]];

rule=MapThread[Rule,{Range[Length[indicesConjugatedMod]],Flatten[Position[indicesOriginal,#]&/@indicesConjugatedMod]}];
invariantsAr[[All,1,pos+1]]=invariantsAr[[All,1,pos+1]]/.rule;

,{pos,pConjs}];

If[pConjs!={},
invariants[[1]]=SparseArray[invariantsAr,Dimensions[invariants[[1]]]];
];


Return[invariants];
]

(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX SelectLines XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)

(*
SelectLines[sellectedLines_]:=Module[{aux},
aux=Sym2Int\[UnderBracket]GenerateListOfCouplings\[UnderBracket]Table;
aux[[1]]=aux[[1,Join[{1},sellectedLines+1]]];
aux[[2,2,2]]=aux[[2,2,2,1;;(Length[sellectedLines]+2)]];
Print[aux];
]
*)

(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX OTHER THINGS XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)
(* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *)


(* Code to check the conservation of a global U1: CheckConservationOfAGlobalU1 and the auxiliar function ContractAt *)
(* Code to check the conservation of a global U1: CheckConservationOfAGlobalU1 and the auxiliar function ContractAt *)
(* Code to check the conservation of a global U1: CheckConservationOfAGlobalU1 and the auxiliar function ContractAt *)

(* contracts the i-th index of tensor with matrix  *)
ContractAt[tensor_,matrix_,i_]:=Module[{indices,result,n},
n=Length[Dimensions[tensor]];
indices=Range[n];
indices[[{i,n}]]={n,i};
result=Transpose[Transpose[tensor,indices].matrix,indices];
Return[result];
]

(* The output eqs are equations which should equal 0. The position of each equation follows the results of GenerateListOfCouplings *)
CheckConservationOfAGlobalU1[result_,nameOfTheFields_]:=Module[{eqs,symMats,aux,tensors},
eqs={};
Do[
tensors=result[[elI,7,1]];

If[Length[tensors]>0,
symMats=Table[Sign[result[[elI,2,i]]]DiagonalMatrix[Array[nameOfTheFields[[Abs[result[[elI,2,i]]]]],Dimensions[tensors][[i+1]]]],{i,Length[result[[elI,2]]]}];
aux=Table[Sum[ContractAt[tensor,symMats[[j]],j],{j,Length[result[[elI,2]]]}],{tensor,tensors}];
,aux={}];
AppendTo[eqs,aux];
,{elI,Length[result]}];

Return[eqs];
]


(* ::Input::Initialization:: *)
(* Input 'term' is the list of field numbers, as it appears in fieldsAUX[model]. Output is {<operator type>, <total derivatives>,<total dummy derivatives>,<interaction was conjugated (False) or not (True)>}. <operator type> is a list of field numbers, as ordered in fields[model], with all derivatives removed, and possibly conjugated *)
OperatorType[model_,term_]:=Block[{termMod,freeDs,appliedDs,formA,formAConj,operatorType,orderedQ,totalDs},
(* Get rid of pure derivatives *)
termMod=DeleteCases[term,derivativeIndex[model]];

freeDs=Count[term,derivativeIndex[model]];
appliedDs=Total[fieldsAUX[model][[Abs[termMod],2,2]]];
totalDs=freeDs+appliedDs;

formA=Sign[termMod]fieldsAUX[model][[Abs[termMod],2,1]];
formAConj=formA(If[#==="C",-1,1]&/@fieldsAUX[model][[Abs[termMod],6]]);

formA=Sort[formA];
formAConj=Sort[formAConj];

If[OrderedQ[{{Count[formA,x_/;x<0],formA},{Count[formAConj,x_/;x<0],formAConj}}],
operatorType=Sort[formA];
orderedQ=True;
,
operatorType=Sort[formAConj];
orderedQ=False;
];

Return[{operatorType,totalDs,freeDs,orderedQ}];
]


(* ::Input::Initialization:: *)
OperatorsOfEachType[model_,data_,calculateSnDataQ_]:=Module[{selfConjQ,operatorTypes,operatorTypeGroups,aux,numberOfOperators,groupedUpTerms,numberOfFreeDs,noNeedToConjugateQ,selfConjTypeQ},

selfConjQ=Table[Sort[resultI[[1]]]==Sort[If[fieldsAUX[model][[Abs[#],6]]=="R",#,-#]&/@resultI[[1]]],{resultI,data}];
operatorTypeGroups=GatherBy[Table[{i,OperatorType[model,data[[i,1]]]},{i,Length[data]}],#[[2,1;;2]]&];

operatorTypes=operatorTypeGroups[[All,1,2,1;;2]];
groupedUpTerms=operatorTypeGroups[[All,All,1]];
numberOfFreeDs=operatorTypeGroups[[All,All,2,3]];
noNeedToConjugateQ=operatorTypeGroups[[All,All,2,4]];


selfConjTypeQ=Table[Sort[resultI[[1]]]==Sort[If[fieldsAUX[model][[Abs[#],6]]=="R",#,-#]&/@resultI[[1]]],{resultI,operatorTypes}];

If[calculateSnDataQ,
aux=Table[data[[el[[i,1]],2,2,All,2]].data[[el[[i,1]],2,2,All,3]](-1)^el[[i,2,3]]If[selfConjQ[[el[[i,1]]]],1,2],{el,operatorTypeGroups},{i,Length[el]}];
numberOfOperators=Simplify[Total/@aux];
numberOfOperators=numberOfOperators/(selfConjTypeQ/.{True->1,False->2});
,
numberOfOperators=ConstantArray[Null,Length[operatorTypeGroups]];
];

Return[{operatorTypes,groupedUpTerms,numberOfFreeDs,noNeedToConjugateQ,numberOfOperators,selfConjTypeQ,selfConjQ[[#]]&/@groupedUpTerms}];
]


DropElementsFromList[list_,elPos_]:=list[[Sort[Complement[Range[Length[list]],elPos]]]]
DropElementsFromLists[lists_,elPos_]:=lists[[All,Sort[Complement[Range[Length[lists[[1]]]],elPos]]]]

(* 1. Delete derivative symmetry information *)
(* Example: derIndex=59; termInfo={{4,4,-4,-4,59,59},{{{1,2},{3,4},{5,6}},{{{{2},{1,1},{2}},1/4 (-1+nn) nn^2 (1+nn),1},{{{1,1},{2},{2}},1/4 (-1+nn) nn^2 (1+nn),1}}},Null} *)
RemoveDerivativeInfo[model_,termInfo_]:=Module[{aux,aux2,convertIndicesRule,result},
aux=Sort[Flatten[Position[termInfo[[1]],derivativeIndex[model]]]];
If[aux==={},Return[termInfo]];

result=termInfo;
result[[1]]=DropElementsFromList[result[[1]],aux];

(* If is not activated if CalculateSnSymmetries\[Rule]False *)
If[result[[2]]=!=Null,

aux2=Flatten[Position[termInfo[[2,1]],aux]];
result[[2,1]]=DropElementsFromList[result[[2,1]],aux2];

result[[2,2,All,1]]=DropElementsFromLists[result[[2,2,All,1]],aux2];

(* sym info grouping at result[[2,1]] may still be incorrect if derivatives we dropped from the middle of result[[1]]. For example, if the fields were {A,Der,B,B}, the groups were {{1},{2},{3,4}} and at this point they would be converted to {{1},{3,4}}. But with the field in positin 2 dropped, it show read {{1},{2,3}}. Correct it now. *)
convertIndicesRule=DropElementsFromList[Range[Length[termInfo[[1]]]],aux];
convertIndicesRule=Table[convertIndicesRule[[i]]->i,{i,Length[convertIndicesRule]}];
result[[2,1]]=result[[2,1]]/.convertIndicesRule;
];

Return[result];
]

FieldToFieldType[model_,field_]:=Sign[field]fieldsAUX[model][[Abs[field],2,1]]


(* ::Input::Initialization:: *)

(* This function is to be applied right after AnalizeInteraction. It changes the output format, and in particular operators which only differ by where the derivatives are applied, are collapsed into a single 'term type'. *)

(* Ouptut is {<termType>, <selfConjTypeQ>,<symInformation>,<total number of operators (real if selfConjTypeQ=True, complex otherwise)>,{<min terms>,<max terms>},<explicitInvariants>}. If <min terms>=<max terms>, then only a single number is shown. *)

ReorganizeOperatorData[model_,operatorTypeData_,operatorsData_,calculateSnDataQ_]:=Module[{derIndex,maxTerms,maxTerms1,maxTerms2,maxTerm3,minTerms,data,fieldTypes,posGroupFields,symReGrouped,aux,aux2,res,output,termType,symInformation,totalOp,explicitInvariants,fldGroups,fusedSymInformation,selfConjTypeValue,selfConjTermValues,nFlavs,posAlwaysValidSnIrreps,repeatedFields,massDimension,selfConjTypeQ},

derIndex=derivativeIndex[model];

selfConjTypeValue=If[operatorTypeData[[6]],1,2];
selfConjTermValues=If[#,1,2]&/@operatorTypeData[[7]];


termType=Join[ConstantArray[derIndex,operatorTypeData[[1,2]]],operatorTypeData[[1,1]]];
massDimension=Total[fieldsAUX[model][[Abs[termType],8]]];
selfConjTypeQ=operatorTypeData[[6]];

(* Invariants are only calculated when CalculateSnSymmetries\[Rule]False, in which case there is just one element in the list operatorsData *)
explicitInvariants=operatorsData[[1,3]];

(* .................................................................................... *)
(* ..................... 1. Remove information about derivatives .....................  *)
(* .................................................................................... *)

data=RemoveDerivativeInfo[model,#]&/@operatorsData;

(* .................................................................................... *)
(* ........... 2. Conjugate terms if that is necessary to make it canonical ........... *)(* ........... (i.e. with type as given by the function OperatorsOfEachType) .......... *)
(* .................................................................................... *)

Do[
If[!operatorTypeData[[4,i]],
data[[i,1]]=data[[i,1]](If[#==="C",-1,1]&/@fieldsAUX[model][[Abs[data[[i,1]]],6]]);
];
,{i,Length[operatorTypeData[[2]]]}];

If[calculateSnDataQ,
(* .................................................................................... *)
(* ..... 3. Group symmetries according to the fields with the derivatives removed ..... *)
(* .................................................................................... *)

fieldTypes=DeleteDuplicates[operatorTypeData[[1,1]]];

maxTerms1=0;
res=Reap[Do[
aux=Flatten[Position[FieldToFieldType[model,data[[i,1]]],#]]&/@fieldTypes;
aux=Table[{j,Position[aux,data[[i,2,1,j,1]]][[1,1]]},{j,Length[data[[i,2,1]]]}];
posGroupFields=SortBy[GatherBy[aux,#[[2]]&],#[[1,2]]&][[All,All,1]];
(* posGroupFields contains the positions of the groups of fields associated to fieldTypes *)

symReGrouped=Table[el[[grp]],{el,data[[i,2,2,All,1]]},{grp,posGroupFields}];
symReGrouped=Transpose[{symReGrouped,(-1)^operatorTypeData[[3,i]]selfConjTermValues[[i]]/selfConjTypeValue data[[i,2,2,All,3]]}];

If[operatorTypeData[[3,i]]==0,
maxTerms1+=Max[data[[i,2,2,All,3]]]selfConjTermValues[[i]]/selfConjTypeValue;
];
Sow[symReGrouped];
,{i,Length[operatorTypeData[[2]]]}]][[2,1]];

(* .................................................................................... *)
(* ................................. 4. Produce output ................................ *)
(* .................................................................................... *)


fldGroups=operatorTypeData[[1,2]]+Flatten[Position[operatorTypeData[[1,1]],#]]&/@fieldTypes;
(* symInformation={fldGroups,TallyWithMultiplicity[Flatten[res,1]]}; *)
symInformation={fldGroups,Flatten[res,1]};


(* ++++++++++++++    keep only the non-trivial symmetry information, ie the one about truly repeated fields ++++++++++++++ *)
aux2=Flatten[Position[symInformation[[1]],x_/;Length[x]>1,{1}]];
symInformation[[2,All,1]]=symInformation[[2,All,1,aux2]];
symInformation[[1]]=symInformation[[1,aux2]];

repeatedFields=termType[[symInformation[[1,All,1]]]];

totalOp=operatorTypeData[[5]];


(* maxTerms1 counts the number of terms needed ignoring operator redundancies. This is done above. *)
(* maxTerms2 counts the total number of irreps of the permutation group after the procedure described in the paper. *)
(* minTerms counts the total number of operators for 1 flavor. If that is negative, then minTerms is forced to be 0. *)

aux=Table[{#[[1]],el[[2]]#[[2]]}&/@TuplesWithMultiplicity[Table[LittlewoodRichardsonCoefficients[elI],{elI,el[[1]]}]],{el,symInformation[[2]]}];
fusedSymInformation=DeleteCases[TallyWithMultiplicity[Flatten[aux,1]],x_/;x[[2]]==0];

(* LittlewoodRichardsonCoefficients ignored the fact that some Sn irreps are irrelevant due to the limited number of flavors. This can even lead to some irreps having negative multiplicities in fusedSymInformation, which is nonsense. Correct this now. *)
nFlavs=fieldsAUX[model][[Abs[repeatedFields],7]];

fusedSymInformation=DeleteCases[fusedSymInformation,x_/;MemberQ[nFlavs-(Length/@x[[1]]),y_Integer/;y<0]];

maxTerms2=Total[fusedSymInformation[[All,2]]];
minTerms=(Times@@@Map[HookContentFormula[#,1]&,fusedSymInformation[[All,1]],{2}]).fusedSymInformation[[All,2]];
minTerms=Max[0,minTerms];

maxTerm3=If[NumericQ[totalOp],totalOp,\[Infinity]];
maxTerms=Min[maxTerms1,maxTerms2,maxTerm3];

(* If there are 0 derivative ... (I) all #flavs are numeric in which case minTerms=maxTerms1, (II) some #flavs are symbolic, in which case assume that they can be as small as 1 and as large as \[Infinity]. *)
If[operatorTypeData[[1,2]]==0,

(* No need for so many parentesis: there are not 'groups of group of repeated fields' if there are no D's *)
symInformation[[2,All,1]]=Map[Flatten,symInformation[[2,All,1]],{2}];

nFlavs=fieldsAUX[model][[Abs[repeatedFields],7]];
 nFlavs=If[NumberQ[#],#,1]&/@nFlavs;
posAlwaysValidSnIrreps=Flatten[Position[symInformation[[2]],x_/;(nFlavs-(Length/@x[[1]]))==Abs[(nFlavs-(Length/@x[[1]]))],{1},Heads->False]];
minTerms=Max[symInformation[[2,posAlwaysValidSnIrreps,2]]];
minTerms=Max[0,minTerms];

maxTerms=maxTerms1;
];

output={termType/.derIndex->0,massDimension,selfConjTypeQ,totalOp,If[minTerms==maxTerms,maxTerms,{minTerms,maxTerms}],repeatedFields,symInformation,explicitInvariants};
,
(* If CalculateSnSymmetries\[Rule]False ... *)
output={termType/.derIndex->0,massDimension,selfConjTypeQ,Null,Null,Null,Null,explicitInvariants};
];

Return[output];
]


(* ::Input::Initialization:: *)
SConsolas[expr_]:=Style[expr,FontFamily->"Consolas"]

DYT[\[Lambda]_,style_]:=If[style==="NoTableaux",SnRepNameTemp2[\[Lambda]],DrawYoungDiagramRaster[\[Lambda],9]]
BuildTableToPrint[model_,data_,calculateSnDataQ_,style_]:=Module[{tableToPrint,item,row,fieldsInInteraction,repeatedFieldsInInteraction,input,output,aux,divsH,divsV,fieldNames,table},

fieldNames=Join[{"\[ScriptCapitalD]"},fields[model][[All,1]],Table["F"<>ToString[i],{i,0,Length[gaugeGroup[model]]}]];
tableToPrint=Reap[Do[
item=data[[i]];

(* If CalculateSnSymmetries\[Rule]False  (item[[5]] is the number of op and it is not calculated in this case) *)
If[!calculateSnDataQ,
row=item[[1;;4]];

fieldsInInteraction=fieldNames[[1+Abs[#]]]<>If[#>=0,"","*"]&/@item[[2]];
row[[2]]=Row[fieldsInInteraction,"  "];
,

row=ConstantArray[Null,8];

row[[{1,3,4,5,6}]]=item[[{1,3,4,5,6}]];

fieldsInInteraction=fieldNames[[1+Abs[#]]]<>If[#>=0,"","*"]&/@item[[2]];
row[[2]]=Row[fieldsInInteraction,"  "];

repeatedFieldsInInteraction=fieldNames[[1+Abs[#]]]<>If[#>=0,"","*"]&/@item[[7]];
repeatedFieldsInInteraction=If[repeatedFieldsInInteraction==={},Null,repeatedFieldsInInteraction];
If[Length[repeatedFieldsInInteraction]==1,repeatedFieldsInInteraction=repeatedFieldsInInteraction[[1]]];
row[[7]]=repeatedFieldsInInteraction;

(* ++++++++++ symmetry information ++++++++++ *)
input=item[[8,2]];
If[item[[8,1]]==={},
output=Null;,
aux=If[!MemberQ[item[[2]],0],Transpose[{input[[All,2]],Map[DYT[#,style]&,input[[All,1]],{2}]}]
,
Transpose[{input[[All,2]],Map[Row[#,"\[Times]"]&,Map[DYT[#,style]&,input[[All,1]],{3}],{2}]}]];
aux[[2;;-1,1]]=Which[#==1,"+",#==-1,"-",#>0,"+"<>ToString[#],#<0,"-"<>ToString[Abs[#]]]&/@aux[[2;;-1,1]];
aux[[1,1]]=Which[aux[[1,1]]==1,"",aux[[1,1]]==-1,"-",aux[[1,1]]>0,ToString[aux[[1,1]]],aux[[1,1]]<0,ToString[aux[[1,1]]]];

(* If there is only one group of repeated fields, remove one set of parentesis *)
If[Length[item[[8,1]]]==1,aux[[All,2]]=aux[[All,2,1]]];

output=Row[DeleteCases[Flatten[aux,1],"",{1}]," "];
];
row[[8]]=output;
];
Sow[row];
,{i,0+ Length[data]}]][[2,1]];

divsH=Join[{1->{Black,Thick},2->Black,-1->{Black,Thick}},Table[i->Lighter[Gray,0.8],{i,3,Length[data]+1}]];

If[calculateSnDataQ,
divsV={1->{Black,Thick},2->Lighter[Gray,0.8],3->Lighter[Gray,0.8],4->Lighter[Gray,0.8],5->Lighter[Gray,0.8],6->Lighter[Gray,0.8],7->Lighter[Gray,0.8],8->Lighter[Gray,0.8],9->{Black,Thick}};

table=Grid[Prepend[tableToPrint,Style[#,{Bold,Darker[Red],FontFamily->"Consolas"}]&/@{"#",Column[{"Operator","type"},Center],"Dim.",Column[{"Self","conj.?"},Center],Column[{"Number of","operators"},Center],Column[{"Number of","terms"},Center],Column[{"Repeated","fields"},Center],Column[{"Permutation","symmetry"},Center]}],Dividers->{divsV,divsH},Spacings->{1,1},Background->{None,{None,{Lighter[Gray,0.95],None}}}];

,
divsV={1->{Black,Thick},2->Lighter[Gray,0.8],3->Lighter[Gray,0.8],4->Lighter[Gray,0.8],5->{Black,Thick}};

table=Grid[Prepend[tableToPrint,Style[#,{Bold,Darker[Red],FontFamily->"Consolas"}]&/@{"#",Column[{"Operator","type"},Center],"Dim.",Column[{"Self","conj.?"},Center]}],Dividers->{divsV,divsH},Spacings->{1,1},Background->{None,{None,{Lighter[Gray,0.95],None}}}];

];

Return[table];
];


ModelStatistics[data_,maxOrder_]:=Module[{opDimensions,nRealOperators,nRealTermsMin,nRealTermsMax,nRealTypes,aux,nOps,nTs,nTTs,statistics,divsH,divsV,tableWithStatistics,frameColor},

nRealOperators={};
nRealTermsMin={};
nRealTermsMax={};
nRealTypes={};
opDimensions={};
Do[
aux=Cases[data,x_/;x[[3]]==dimOp];
nOps=Simplify[aux[[All,5]].(aux[[All,4]]/.{True->1,False->2})];
nTs=Total[aux[[All,6]](aux[[All,4]]/.{True->1,False->2})];
nTTs=Total[aux[[All,4]]/.{True->1,False->2}];
If[NumericQ[nTs],nTs={nTs,nTs}];

AppendTo[opDimensions,dimOp];
AppendTo[nRealOperators,nOps];
AppendTo[nRealTermsMin,nTs[[1]]];
AppendTo[nRealTermsMax,nTs[[-1]]];
AppendTo[nRealTypes,nTTs];
,{dimOp,Min[data[[All,3]]],maxOrder}];

statistics=Transpose[{opDimensions,nRealOperators,nRealTermsMin,nRealTermsMax,nRealTypes}];

aux=Table[{el[[1]],el[[2]],If[el[[3]]==el[[4]],el[[3]],Row[{el[[3]]," to ",el[[4]]}]],el[[5]]},{el,statistics}];
PrependTo[aux,Style[#,{Bold,Darker[Green],FontFamily->"Consolas"}]&/@{"Dimension","# real operators","# real terms", "# types of real operators"}];

frameColor=Black;
divsH=Join[{1->{frameColor,Thick},2->frameColor,-1->{frameColor,Thick}},Table[i->Lighter[Gray,0.8],{i,3,Length[statistics]+1}]];
divsV={1->{frameColor,Thick},2->Lighter[Gray,0.8],3->Lighter[Gray,0.8],4->Lighter[Gray,0.8],5->{frameColor,Thick}};

tableWithStatistics=Grid[aux,Dividers->{divsV,divsH},Spacings->{1,1},Background->{None,{None,{Lighter[Gray,0.95],None}}}];

Return[{statistics,tableWithStatistics}];
]


(* Prints data obtained with GenerateListOfCouplings in a table. The advantage is that one can select only part of the full data *)
PrintOperatorTable[model_,data_]:=BuildTableToPrint[model,data,True,True]
