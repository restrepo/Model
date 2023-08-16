Off[General::spell]

Model`Name = "Zprime"; (*Must coincide with directory name*)
Model`NameLaTeX ="Singlet scalar Dark Matter with complex even singlet scalar";
Model`Authors = "Diego Restrepo (based on SM model by F.Staub)";
Model`Date = "2015-11-16";
(*Reads configuration file in both unix and windows*)
Get[FileNameJoin[{$sarahDir,"Models",Model`Name,"config.m"}]];
(* Adjust configuration variables *)
If [GaugeU1,
    EvenSingletScalar = True;
];

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
If [GaugeU1,
    (* name of the charge needs at least 3 characters *)
    Gauge[[4]]={Bp,  U[1], XXX,       g1p, False,1}; (* False as in the official B-L Model *)
(* Get charges from Module *)
  anomaly[d_] :=
  Module[{X, Y, YL, S, NR = -1},
    S = -2 NR;
    X = {Xu, Xu, Xu, d, d, d, -Xq, -Xq, -Xq, -Xq, -Xq, -Xq, -Xl, -Xl, 
      Xe, NR};
    Y = Simplify[
      Solve[{2 d + Xu == 1, -Xq + XH + d == 0, -Xq - XH + Xu == 
        0, -Xl + XH + Xe == 0,
        Total[X] == 0}, {Xq, Xu, Xl, Xe, XH}]];
    {XXq, XXu, XXl, XXe, XXH} = {Xq, Xu, Xl, Xe, XH} /. Y[[1]];
    (*Y=Join[{Xd->d},Y[[1]],{Xbi->S,Xz->NR,nMG->3}];*)
    YL = {Xq -> XXq, Xl -> XXl, Xu -> -XXu, Xd -> -d, Xe -> -XXe, 
      XH -> XXH, Xbi -> -S, Xz -> -NR, nMG -> 3};
    Return[YL]
    ]
  (****
  {Xq, Xl, Xd, Xu, Xe, XH, Xbi} /. anomaly[1/3]
  {1/3, -1, -(1/3), -(1/3), 1, 0, 2}
  *****)      
];

(*{Xq,Xl,Xd,Xu,Xe,XH,Xbi}={0, 0, 0, 0, 0, 0, 2};*)
{Xq,Xl,Xd,Xu,Xe,XH,Xbi}={1/3, -1, -(1/3), -(1/3), 1, 0, -2}

(* Matter Fields *)

(* SM fermion fields *)
{NS,nF} = {1,1};
(* Even fermion fields *)
FermionFields[[nF]] = {q, 3, {uL, dL},     1/6, 2,  3};
If [GaugeU1,
    FermionFields[[nF]] = Join[FermionFields[[nF]],{Xq}];
];
nF = nF + 1;
FermionFields[[nF]] = {l, 3, {vL, eL},    -1/2, 2,  1};
If [GaugeU1,
    FermionFields[[nF]] = Join[FermionFields[[nF]],{Xl}];
];
nF = nF + 1;
FermionFields[[nF]] = {d, 3, conj[dR],     1/3, 1, -3};
If [GaugeU1,
    FermionFields[[nF]] = Join[FermionFields[[nF]],{Xd}];
];
nF = nF + 1;
FermionFields[[nF]] = {u, 3, conj[uR],    -2/3, 1, -3};
If [GaugeU1,
    FermionFields[[nF]] = Join[FermionFields[[nF]],{Xu}];
];
nF = nF + 1;
FermionFields[[nF]] = {e, 3, conj[eR],       1, 1,  1};
If [GaugeU1,
    FermionFields[[nF]] = Join[FermionFields[[nF]],{Xe}];
];
(* Z2 charges *)
Do [ 
  FermionFields[[i]] = Join[FermionFields[[i]],{1}];, {i,1,nF}
   ];
(* Odd fermion fields *)

(*Even scalar fields*)
ScalarFields[[NS]] =  {H, 1, {Hp, H0},     1/2, 2,  1};
If [GaugeU1,
    ScalarFields[[NS]] = Join[ScalarFields[[NS]],{XH}];
];
If [EvenSingletScalar,
    NS = NS + 1;
    ScalarFields[[NS]] =  {bi,1, BiD,    0, 1,  1};
    If [GaugeU1,
        ScalarFields[[NS]] = Join[ScalarFields[[NS]],{Xbi}];
    ];
];
Do [ 
  ScalarFields[[i]] = Join[ScalarFields[[i]],{1}];, {i,1,NS}
   ];

(* Odd neutral scalar fields*)
NS = NS + 1;
NSO = NS;
ScalarFields[[NS]] =  {S, 1, ss,     0, 1,  1};
If [GaugeU1,
    ScalarFields[[NS]] = Join[ScalarFields[[NS]],{0}];
];
Do [ 
  ScalarFields[[i]] = Join[ScalarFields[[i]],{-1}];, {i,NSO,NS}
   ];

RealScalars = {S};

(* XXX charges initialization *)
{Xn,Xp,Xr,Xs,Xt,Xw,Xx,Xy,Xz}={0,0,0,0,0,0,0,0,0};
{xMajorana,yDirac,yMajorana}={False,False,False};

(******* BEGIN: XXX-charged BSM chiral or vector-like fermion fields *******)
If[GaugeU1,
  (** 
    Anomaly solution: Capital letters may include generations
    {D,i,r,a,c,0,...,m,a,j,o,r,A,n,    A}: 
     D->nDG,                   A->nWG, A -> nMG 
     {nDG,nWG,nMG} from config.m
  **)
  {Xn,Xp,Xr,Xs,Xt,Xw,Xx,Xy,Xz}={1/5,-1/5,0,0,0,0,0,0,0};

  (* Multi-generation Dirac Fermions -> Fix PDG numbers in particles.m *)
  If[Xn != 0 && Xp != 0,
    nF=nF+1;
    FermionFields[[nF]] = {n, nDG, nL,	    0, 1,  1, Xn, -1};
    nF=nF+1;
    FermionFields[[nF]] = {p, nDG, conj[pR], 0, 1,  1, Xp, -1};
  ];
];

  (* Single family Dirac Fermion *)
  If[Xr != 0 && Xs != 0,
    nF=nF+1;
    FermionFields[[nF]] = {r, 1, rL,	    0, 1,  1, Xr, -1};
    nF=nF+1;   
    FermionFields[[nF]] = {s, 1, conj[sR],   0, 1,  1, Xs, -1};
  ];
  (* Single family Dirac Fermion *)
  If[Xt != 0 && Xw != 0,
    nF=nF+1;
    FermionFields[[nF]] = {t, 1, tL,	    0, 1,  1, Xt, -1};
    nF=nF+1;   
    FermionFields[[nF]] = {w, 1, conj[wR],   0, 1,  1, Xw, -1};
  ];

  (* (Xx,Xy) pair cases *)
  (************** BEGIN CASES *****************************)
  (* Cases: [1,-1]; [1,1] or [-1,-1] (S=2); [1,3] or [-1,-3] (S=4); [1,3] (S=6); [1,0] (S=2)   *)
  If[Xx != 0 || Xy != 0,
    xMajorana = False;
    yDirac = False;
    yMajorana = False;

    If [
          Xx + Xy == 0, (*vector-like*)
              yDirac = True;,
          If [ (*non vector-like*)
                Abs[Xx + Xy] == Abs[Xbi], (*Massive chiral*)
                If[ (*chiral*)
                    Xx != Xy,
                      yDirac = True;,
                      yMajorana = True;
                    ],
                    (*massles Dirac*)
                    If [
                          Abs[2 Xy] == Abs[Xbi],
                            yMajorana = True;,
                          If [
                                Abs[2 Xx] == Abs[Xbi],
                                  xMajorana = True;
                            ]                           
                      ]             
              ]
        ]
    ]
  (* Cases implementation *)
  If[Xx != 0 && Xy != 0 && yDirac,
    nF=nF+1;
    FermionFields[[nF]] = {x, 1, tL,	    0, 1,  1, Xx, -1};
    nF=nF+1;   
    FermionFields[[nF]] = {y, 1, conj[yR],   0, 1,  1, Xy, -1};
  ];

  (* Multi-generation Weyl Fermion -> Fix PDG numbers in particles.m (Dirac neutrinos) *)
  If[Xx != 0 && (!yDirac && !yMajorana && !xMajorana),
    nF=nF+1;
    FermionFields[[nF]]  = {x, nWG, xL,	    0, 1,  1, Xx, -1};
    ];

  (* Single Majorana Fermion *)
  If[Xy != 0 && xMajorana,
    nF=nF+1;
    FermionFields[[nF]]  = {y, 1, yL,	    0, 1,  1, Xy, -1}; (* try conj[yR] if errors *)
    ];
  (************** END CASES *******************************)
  (*TODO: If SM Dirac Yukawa redefines the mixing nad 4-spinor*)
  If[Xz != 0,
    nF=nF+1;
    FermionFields[[nF]]  = {z, nMG, zL,	    0, 1,  1, Xz, -1}; (* try conj[zR] if errors *)
      ];

(******* END: XXX-charged BSM chiral or vector-like fermion fields *********)

        
(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput] = {
	{LagFer,    {AddHC->True}},
	{LagNoHC,  {AddHC->False}}
};

If [EvenSingletScalar,    
  DEFINITION[GaugeES][LagrangianInput] =
    Join[DEFINITION[GaugeES][LagrangianInput],
        {{LagNoHCbi,{ AddHC->False}}}
    ];
  ];


LagNoHC = -(mH2 conj[H].H + Lambda1/2 conj[H].H.conj[H].H + MS2/2 S.S + LamSH S.S.conj[H].H  + LamS/2 S.S.S.S);
LagFer =  -(Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q);
If [EvenSingletScalar,
    LagNoHCbi = -(MuP conj[bi].bi - L2 conj[bi].bi.conj[bi].bi - L3 conj[bi].bi.conj[H].H );
];

(*** BEGIN: Extends LagFer with chiral XXX-charged fermions ***)
(* (Xn,Xp) cases*)
(*Massive chiral Dirac fermion*)
If[ Xn !=0 && Xn + Xp + Xbi == 0,
    LagFer = LagFer + Ynp n.p.bi;
    ];
(*Vector like Dirac fermion*)
If[ Xn !=0 && Xn + Xp == 0,
    LagFer = LagFer + Mnp n.p;
    ];

If[ Xr !=0 && Xr + Xs - Xbi == 0,
    LagFer = LagFer + Yrs r.s.conj[bi];
    ];
If[ Xt !=0 && Xt + Xw - Xbi == 0,
    LagFer = LagFer + Ytv t.w.conj[bi];
    ];

If[ yDirac && Xx + Xy + Xbi == 0,
 LagFer = LagFer +  Yxy x.y.bi;
];
If[ yMajorana && 2 Xy + Xbi == 0,
 LagFer = LagFer +  Yxy y.y.bi;
];
If[ xMajorana && 2 Xx + Xbi == 0,
 LagFer = LagFer +  Yxy x.x.bi;
];

If[ Xz !=0 && 2*Xz  + Xbi == 0,
 LagFer = LagFer +  Yz z.z.bi;
];

(*** END: Extends LagFer with chiral XXX-charged fermions *****)


(* Gauge Sector *)

If [GaugeU1,
  DEFINITION[EWSB][GaugeSector] = { 
    {{VB,VWB[3],VBp},{VP,VZ,VZp},ZZ},  
    {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
  };,
  DEFINITION[EWSB][GaugeSector] = { 
    {{VB,VWB[3]},{VP,VZ},ZZ},  
    {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
  };
];


(* Last safe place to implement code here. Define lists to be used later on *)
If [EvenSingletScalar,
  EWSBMatterSectorList =
  {
     {{phiH,phiB},{hh,ZH}},
     {{sigmaH,sigmaB},{Ah,ZA}},
     {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
     {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
     {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}
    };,
  EWSBMatterSectorList =
  {
     {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
     {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
     {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}
    };
];

(*** BEGIN: Extends EWSBMatterSectorList with chiral XXX-charged fermions ***)
If[Xn != 0 && nDG > 1,
  EWSBMatterSectorList = Join[EWSBMatterSectorList,
    {
      {{{nL}, {conj[pR]}}, {{NL,Vn}, {NR,Un}}}
    }
   ];
];
If[Xz != 0 && nMG > 1,
   EWSBMatterSectorList = Join[EWSBMatterSectorList,,
      {
	      {{zL},{ZL, Uz}}	 
      }
				       ];
 ];
(*** END: Extends EWSBMatterSectorList with chiral XXX-charged fermions ***)

(*** BEGIN: 4-spinor definitions for XXX-charged chiral fields ***)
If[Xn != 0,
      If[nDG > 1,
      FnList={Fn ->{  NL, conj[NR]}};,
      FnList={Fn ->{  nL, conj[pR]}};
      ]
   ];

If[Xz != 0 && nDM > 1,
      FzList={Fz ->{  ZL, conj[ZL] }};
   ];
(*** END: 4-spinor definitions for XXX-charged chiral fields ***)


(* WARNING: Avoid code after here *)



(* ----- VEVs ---- *)

If [EvenSingletScalar,
 DEFINITION[EWSB][VEVs] = 
{    {H0, {vH, 1/Sqrt[2]},  {sigmaH, \[ImaginaryI]/Sqrt[2]},{phiH, 1/Sqrt[2]}},
     {BiD,{vX, 1/Sqrt[2]}, {sigmaB, \[ImaginaryI]/Sqrt[2]},{phiB, 1/Sqrt[2]}}
  };,
DEFINITION[EWSB][VEVs] = 
{    {H0, {vH, 1/Sqrt[2]}, {Ah, \[ImaginaryI]/Sqrt[2]},{hh, 1/Sqrt[2]}}     };
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

If[Xn != 0,
   DEFINITION[EWSB][DiracSpinors]=Join[
         DEFINITION[EWSB][DiracSpinors],
         FnList   
                                      ];
   ]; 

If[Xz != 0,
   DEFINITION[EWSB][DiracSpinors]=Join[
         DEFINITION[EWSB][DiracSpinors],
         FzList   
                                      ];
   ]; 


DEFINITION[EWSB][GaugeES]={
 Fd1 ->{  FdL, 0},
 Fd2 ->{  0, FdR},
 Fu1 ->{  Fu1, 0},
 Fu2 ->{  0, Fu2},
 Fe1 ->{  Fe1, 0},
 Fe2 ->{  0, Fe2}};



