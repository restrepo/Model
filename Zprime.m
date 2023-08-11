Off[General::spell]

Model`Name = "Zprime";
Model`NameLaTeX ="Singlet scalar Dark Matter";
Model`Authors = "Diego Restrepo (based on SM model by F.Staub)";
Model`Date = "2015-11-16";

(* 2013-01-24: changed normalization of lambda term to convention of hep-ph/0207271 *)
(* 2013-06-24: using new name conventions (without inital "S" and "F" for scalar and matter fields) *)
(* 2013-09-01: changing to new conventions for FermionFields/MatterFields *)
(* 2013-11-20: Singlet Scalar DM implemented *)
(* 2014-11-06: Changed sign in Lagrangian *)
(* 2015-11-16: changed SPheno.m *)

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Global Symmetries *)

Global[[1]] = {Z[2], Z2};


(* Gauge Groups *)


Gauge[[1]]={B,   U[1], hypercharge, g1,False,1};
Gauge[[2]]={WB, SU[2], left,        g2,True,1};
Gauge[[3]]={G,  SU[3], color,       g3,False,1};


(* Matter Fields *)
EvenSingletScalar = True; (*TODO: move to config*)
NS = 1;

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3,1};  
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1,1};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3,1};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3,1};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1,1};

ScalarFields[[NS]] =  {H, 1, {Hp, H0},     1/2, 2,  1,1};
NS = NS + 1;
ScalarFields[[NS]] =  {S, 1, ss,     0, 1,  1, -1};
NS = NS + 1;
If [EvenSingletScalar,
    ScalarFields[[NS]] =  {bi,1, BiD,    0, 1,  1, 1};
];

RealScalars = {S};



        
(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput] = {
	{LagHC,    {AddHC->True}},
	{LagNoHC,  {AddHC->False}}
};

If [EvenSingletScalar,    
  DEFINITION[GaugeES][LagrangianInput] =
    Join[DEFINITION[GaugeES][LagrangianInput],
        {{LagNoHCbi,{ AddHC->False}}}
    ];
  ];


LagNoHC = -(mH2 conj[H].H + lambda1/2 conj[H].H.conj[H].H + MS2/2 S.S + LamSH S.S.conj[H].H  + LamS/2 S.S.S.S);
LagHC =  -(Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q);
If [EvenSingletScalar,
    LagNoHCbi = -(MuP conj[bi].bi - L2 conj[bi].bi.conj[bi].bi - L3 conj[bi].bi.conj[H].H );
];




			  		  

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
};     
        

If [EvenSingletScalar,
  EWSBMatterSectorList =
  {
     {{phiH,phiB},{hh,ZH}},
     {{sigmaH,sigmaB},{Ah,ZA}},
     {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
     {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
     {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}
    };
];


(* WARNING: Avoid code after here *)



(* ----- VEVs ---- *)

If [EvenSingletScalar,
 DEFINITION[EWSB][VEVs] = 
{    {H0, {vH, 1/Sqrt[2]},  {sigmaH, \[ImaginaryI]/Sqrt[2]},{phiH, 1/Sqrt[2]}},
     {BiD,{vX, 1/Sqrt[2]}, {sigmaB, \[ImaginaryI]/Sqrt[2]},{phiB, 1/Sqrt[2]}}
  };,
DEFINITION[EWSB][VEVs] = 
{    {H0, {vH, 1/Sqrt[2]}, {Bh, \[ImaginaryI]/Sqrt[2]},{Rh, 1/Sqrt[2]}}     };
];



DEFINITION[EWSB][MatterSector] = EWSBMatterSectorList;

(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  DL, conj[DR]},
 Fe ->{  EL, conj[ER]},
 Fu ->{  UL, conj[UR]},
 Fv ->{  vL, 0}};

DEFINITION[EWSB][GaugeES]={
 Fd1 ->{  FdL, 0},
 Fd2 ->{  0, FdR},
 Fu1 ->{  Fu1, 0},
 Fu2 ->{  0, Fu2},
 Fe1 ->{  Fe1, 0},
 Fe2 ->{  0, Fe2}};



