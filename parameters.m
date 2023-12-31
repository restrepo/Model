(*Reads configuration file in both unix and windows*)
Get[FileNameJoin[{$sarahDir,"Models",Model`Name,"config.m"}]];
If [GaugeU1,
    EvenSingletScalar = True;,
    seesaw = False;
];



If[EvenSingletScalar,
      Lambda1list={Lambda1,  {OutputName -> lam1,
                              LaTeX -> "\\lambda_1",
                              LesHouches -> {BL,1}}};,
      Lambda1list={Lambda1, { Description -> "SM Higgs Selfcouplings",
                              DependenceNum -> Mass[hh]^2/(vH^2)}};
];

If [GaugeU1,
    ZZlist = {ZZ, {Description ->   "Photon-Z-Z' Mixing Matrix"}};
    ,
    ZZlist = {ZZ, {Description -> "Photon-Z Mixing Matrix"}};
];

ParameterDefinitions = { 

{g1,        { Description -> "Hypercharge-Coupling"}},
{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    

{g11p,        {Description -> "Mixed Gauge Coupling 2"}},
{g1p1,        {Description -> "Mixed Gauge Coupling 1"}},
{g1p,       {   Description -> "B-L-Coupling"}},

{MZp,       {   Description -> "Z' mass"}},

{AlphaS,    {Description -> "Alpha Strong"}},	
{e,         { Description -> "electric charge"}}, 

{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},

{Yu,        { Description -> "Up-Yukawa-Coupling",
			 DependenceNum ->  Sqrt[2]/vH* {{Mass[Fu,1],0,0},
             									{0, Mass[Fu,2],0},
             									{0, 0, Mass[Fu,3]}}}}, 
{Yd,        { Description -> "Down-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/vH* {{Mass[Fd,1],0,0},
             									{0, Mass[Fd,2],0},
             									{0, 0, Mass[Fd,3]}}}},
{Ye,        { Description -> "Lepton-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/vH* {{Mass[Fe,1],0,0},
             									{0, Mass[Fe,2],0},
             									{0, 0, Mass[Fe,3]}}}},                                                                                                             
{vH,          { Description -> "EW-VEV",
               DependenceNum -> Sqrt[4*Mass[VWp]^2/(g2^2)],
               DependenceSPheno -> None  }},

{vX,      {  LaTeX -> "v_x",
             Dependence ->  None, 
             OutputName -> vX,
             Real -> True,
             LesHouches -> {BL,43} }},

{ThetaW,    { Description -> "Weinberg-Angle",
              DependenceNum -> ArcSin[Sqrt[1 - Mass[VWp]^2/Mass[VZ]^2]]  }},

{ThetaWp,  { Description -> "Theta'", DependenceNum -> None  }},

ZZlist,

{ZW, {Description -> "W Mixing Matrix",
       Dependence ->   1/Sqrt[2] {{1, 1},
                  {\[ImaginaryI],-\[ImaginaryI]}} }},

{L2, {OutputName -> lam2,
      LaTeX -> "\\lambda_2",
      LesHouches -> {BL,2}}},

{L3, {OutputName -> lam3,
      LaTeX -> "\\lambda_3",
      LesHouches -> {BL,3}}},

{MuP, {OutputName -> MUP,
      LaTeX -> "\\mu'",
      LesHouches -> {BL,10}}},

{Ynp, {OutputName -> Ynp,
      LaTeX -> "Y_{np}",
      LesHouches -> YNP }},

{Yrs, {OutputName -> Yrs,
      LaTeX -> "Y_{rs}",
      LesHouches -> YRS }},
						       
{Ytv, {OutputName -> Ytv,
      LaTeX -> "Y_{tv}",
      LesHouches -> YTV }},

{Yxy, {OutputName -> Yxy,
      LaTeX -> "Y_{xy}",
      LesHouches -> YXY }},

{Yz, {OutputName -> Yz,
      LaTeX -> "Y_{z}",
      LesHouches -> YZ }},

{Yv, {OutputName -> Yv,
      LaTeX -> "Y_\\nu",
      LesHouches -> Yv}},      

{Mnp, {OutputName -> Mnp,
      LaTeX -> "M_{np}",
      LesHouches -> MNP }},

{Vu,        {Description ->"Left-Up-Mixing-Matrix"}},
{Vd,        {Description ->"Left-Down-Mixing-Matrix"}},
{Uu,        {Description ->"Right-Up-Mixing-Matrix"}},
{Ud,        {Description ->"Right-Down-Mixing-Matrix"}}, 
{Ve,        {Description ->"Left-Lepton-Mixing-Matrix"}},
{Ue,        {Description ->"Right-Lepton-Mixing-Matrix"}},



						       
(* Scalar sector *)
{mH2,         { Description -> "SM Mu Parameter",
                LaTeX -> "\\mu^2",
	        OutputName -> mH2}},

{ZH,        { Description->"Scalar-Mixing-Matrix", 
               Dependence -> None,
               DependenceOptional -> None,
               DependenceNum -> None   }},
{ZA,        { Description->"Pseudo-Scalar-Mixing-Matrix", 
              Dependence -> None,
              DependenceOptional -> None,
             DependenceNum -> None   }},

{Vn, {LaTeX -> "V^n_R",
      LesHouches -> VNRMIX,
      OutputName -> ZNR }},

{Un, {LaTeX -> "U^n_L",
      LesHouches -> UNLMIX,
      OutputName -> ZNL }},

{Uz, {LaTeX -> "U^z_L",
      LesHouches -> UZLMIX,
      OutputName -> ZZL }},

{MS2,       {Description -> "Softbreaking Up-Higgs Mass",
             LaTeX -> "M_S^2",
	     OutputName-> MS2}},

{ZM,	    {Description -> "Neutrino-Mixing-Matrix"}},           

Lambda1list,

{LamSH,     {OutputName ->"LSH",
             LesHouches -> {HDM,2}}},

{LamS,     {OutputName ->"LS",
             LesHouches -> {HDM,3}}}


 }; 
 

