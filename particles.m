(*Reads configuration file in both unix and windows*)
Get[FileNameJoin[{$sarahDir,"Models",Model`Name,"config.m"}]];
If [GaugeU1,
    EvenSingletScalar = True;,
    seesaw = False;
];


  If [EvenSingletScalar,
        Ahlist = {Ah   ,  {  Description -> "Pseudo-Scalar Higgs"}};
        hhlist = {hh   ,  { Description -> "Higgs",
                 PDG -> {25,35},
                 Width -> Automatic, 
                 Mass ->LesHouches,
                 FeynArtsNr -> 1,
                 LaTeX -> "h",
                 OutputName -> "h" }};
        VZlist = {VZ,   { Description -> "Z-Boson",
      			 Goldstone -> Ah[{1}] }};,
        Ahlist = {Ah   ,  {  Description -> "Pseudo-Scalar Higgs",
                 FeynArtsNr -> 2,
                 PDG -> {0},
                 PDG.IX ->{0},
                 Mass -> {0},
                 Width -> {0} }};
        hhlist = {hh   ,  {  Description -> "Higgs",
                PDG -> {25},
		PDG.IX -> {101000001},
                Mass -> LesHouches,
                FeynArtsNr -> 1,
                LaTeX -> "h",
                ElectricCharge -> 0,
                LHPC -> {1},
		OutputName -> "h"  }};
        VZlist = {VZ,   { Description -> "Z-Boson",
      			 Goldstone -> Ah }};                
     ];
If [GaugeU1,
   VZplist = {VZp,    { Description -> "Z'-Boson",
 		             Goldstone -> Ah[{2}]}
             };
   gZplist={gZp,    { Description -> "Z'-Ghost" }};
   ,
   VZplist = {VZp,    {LaTeX -> "None" }
             };
   gZplist={gZp,    { LaTeX -> "None" }};                  
];
If [seesaw,
    FvList =       {Fv,   { Description -> "Neutrinos",
 			PDG ->{12,14,16,8810012,8810014,8810016}
                          }
                   };,
   FvList =       {Fv,   { Description -> "Neutrinos"}
                  };
]
FnPDGList = {1000001};
FnPDGIXList = {2000001};
If [nDG >= 2,
   FnPDGList = {1000001,1000002};
   FnPDGIXList = {2000001,2000002};
];
If [nDG >= 3,
   FnPDGList = {1000001,1000002,1000003};
   FnPDGIXList = {2000001,2000002,2000003};
];

FnList = {Fn,   { Description -> "Darkinos",
      PDG -> FnPDGList,
      PDG.IX -> FnPDGIXList,
      Mass -> LesHouches,
         FeynArtsNr -> 11,
         LaTeX -> "{\\chi}^0_n",
         ElectricCharge -> 0,
         LHPC -> {3, "orange"},
         OutputName -> "N"}
 };



ParticleDefinitions[GaugeES] = {
      {H0,  {    PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> "H^0",
                 OutputName -> "H0" }},                         
      
      
      {Hp,  {    PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 2,
                 LaTeX -> "H^+",
                 OutputName -> "Hp" }}, 
                 
               
    
      {VB,   { Description -> "B-Boson"}},                                                   
      {VG,   { Description -> "Gluon"}},          
      {VWB,  { Description -> "W-Bosons"}},          

      {VBp,   { Description -> "Bp-Boson",
                PDG -> {0},
                PDG.IX -> {0},
                Width -> 0, 
                Mass -> 0,
                FeynArtsNr -> 1,
                LaTeX -> "B'",
                OutputName -> "Bp" }},
      
      {gB,   { Description -> "B-Boson Ghost"}},
      {gG,   { Description -> "Gluon Ghost" }},          
      {gWB,  { Description -> "W-Boson Ghost"}},
      {gBp,  { Description -> "Bp-Boson Ghost",
               PDG -> {0},
               PDG.IX -> {0},
               Width -> 0, 
               Mass -> 0,
               FeynArtsNr -> 1,
               LaTeX -> "\\eta^Bp",
               OutputName -> "gBp"}}      
      
      };
      

      
      
  ParticleDefinitions[EWSB] = {
    hhlist, 

    {ss   , {  Description -> "Singlet",
	       PDG -> {6666635},
	       PDG.IX -> {101000002},
               FeynArtsNr -> 10,
               Mass -> LesHouches,               
               LaTeX -> "S",
               ElectricCharge -> 0,
               LHPC -> {"gold"},
       	       OutputName -> "Ss" }},

    
    Ahlist,

    
      
     {Hp,     { Description -> "Charged Higgs", 
                 PDG -> {0},
                 PDG.IX ->{0},
                 Width -> {0}, 
                 Mass -> {0},
                 ElectricCharge->+1,                 
                 LaTeX -> {"H^+","H^-"},
                 OutputName -> {"Hp","Hm"}
                 }},                                                   
      
      {VP,   { Description -> "Photon"}},
      VZlist,
      {VG,   { Description -> "Gluon" }},          
      {VWp,  { Description -> "W+ - Boson",
      			Goldstone -> Hp }},
      VZplist,         
      {gP,   { Description -> "Photon Ghost"}},                                                   
      {gWp,  { Description -> "Positive W+ - Boson Ghost"}}, 
      {gWpC, { Description -> "Negative W+ - Boson Ghost" }}, 
      {gZ,   { Description -> "Z-Boson Ghost" }},
      {gG,   { Description -> "Gluon Ghost" }},          
      gZplist,                            
      {Fd,   { Description -> "Down-Quarks"}},   
      {Fu,   { Description -> "Up-Quarks"}},   
      {Fe,   { Description -> "Leptons" }},
      FvList,
      FnList,
      {Fr,   { Description -> "Darkinos",
	      PDG -> {1000011},
	      PDG.IX -> {2000011},
              FeynArtsNr -> 12,
              ElectricCharge -> 0,
              LaTeX -> "{\\chi}^0_r",
              OutputName -> "R"}
      },
      {Ft,   { Description -> "Darkinos",
	      PDG -> {1000021},
	      PDG.IX -> {2000021},
              FeynArtsNr -> 13,
              ElectricCharge -> 0,              
              LaTeX -> "{\\chi}^0_t",
              OutputName -> "T"}
      },
      {Fx,   { Description -> "Darkinos",
	      PDG -> {1000031},
	      PDG.IX -> {2000031},
              FeynArtsNr -> 14,
              ElectricCharge -> 0,              
              LaTeX -> "{\\chi}^0_x",
              OutputName -> "X"}
      },
      {Fz,   { Description -> "Darkinos",
	      PDG -> {1000051,1000052,1000053},
	      PDG.IX -> {2000051,2000052,2000053},
              FeynArtsNr -> 15,
              ElectricCharge -> 0,
              LaTeX -> "{\\chi}^0_x",
              OutputName -> "X"}
       }
       };    
        
 WeylFermionAndIndermediate = {
     
    {H,      {   PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 LaTeX -> "H",
                 OutputName -> "" }},

    {S,      {   PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 LaTeX -> "S",
                 OutputName -> "" }},

   {dR,     {LaTeX -> "d_R" }},
   {eR,     {LaTeX -> "e_R" }},
   {lep,     {LaTeX -> "l" }},
   {uR,     {LaTeX -> "u_R" }},
   {q,      {LaTeX -> "q" }},
   {eL,     {LaTeX -> "e_L" }},
   {dL,     {LaTeX -> "d_L" }},
   {uL,     {LaTeX -> "u_L" }},
   {vL,     {LaTeX -> "\\nu_L" }},

   {DR,     {LaTeX -> "D_R" }},
   {ER,     {LaTeX -> "E_R" }},
   {UR,     {LaTeX -> "U_R" }},
   {EL,     {LaTeX -> "E_L" }},
   {DL,     {LaTeX -> "D_L" }},
   {UL,     {LaTeX -> "U_L" }},

   {nL,     {LaTeX -> "n_L" }},
   {pR,     {LaTeX -> "p_R" }},
   {zL,     {LaTeX -> "z_L" }},
   
   {NR,     {LaTeX -> "N_R" }},
   {NL,     {LaTeX -> "N_R" }},
   {ZL,     {LaTeX -> "Z_R" }}
   
        };       


