(*Reads configuration file in both unix and windows*)
Get[FileNameJoin[{$sarahDir,"Models",Model`Name,"config.m"}]];
OnlyLowEnergySPheno = True;


MINPAR={{1,Lambda1INPUT},
        {2,LamSHINPUT},
        {3,LamSINPUT},
        {4,MSinput}
        };
If[EvenSingletScalar,
  MINPAR = Join[MINPAR,
              {{5,Lambda2INPUT},
               {6,Lambda3INPUT},
               {20,vXINPUT}
               }
             ];
  ];

ParametersToSolveTadpoles = {mH2};
If[EvenSingletScalar,
  ParametersToSolveTadpoles = Join[ParametersToSolveTadpoles,
                                   {MuP}];
  ];
  

DEFINITION[MatchingConditions]= {
 {Ye, YeSM},
 {Yd, YdSM},
 {Yu, YuSM},
 {g1, g1SM},
 {g2, g2SM},
 {g3, g3SM},
 {vH, vSM}
 };

BoundaryLowScaleInput={
  {Lambda1,Lambda1INPUT},
  {LamSH,LamSHINPUT},
  {LamS,LamSINPUT},
  {MS2, MSInput}
};
If [EvenSingletScalar,
  BoundaryLowScaleInput = Join[BoundaryLowScaleInput,
                              {  {L2, Lambda2INPUT},
                                 {L3, Lambda3INPUT},
                                 {vX,vXINPUT}
                               }
                              ];
  ];


ListDecayParticles = {Fu,Fe,Fd};
If[EvenSingletScalar,
  ListDecayParticles = Join[ListDecayParticles,
                            {Ah,hh}];,
  ListDecayParticles = Join[ListDecayParticles,
                            {hh}];
  ];

ListDecayParticles3B = {{Fu,"Fu.f90"},{Fe,"Fe.f90"},{Fd,"Fd.f90"}};



DefaultInputValues ={Lambda1INPUT -> 0.28};

AddTreeLevelUnitarityLimits=True;
FlagLoopContributions = True;

RenConditionsDecays={
{dCosTW, 1/2*Cos[ThetaW] * (PiVWp/(MVWp^2) - PiVZ/(mVZ^2)) },
{dSinTW, -dCosTW/Tan[ThetaW]},
{dg2, 1/2*g2*(derPiVPheavy0 + PiVPlightMZ/MVZ^2 - (-(PiVWp/MVWp^2) + PiVZ/MVZ^2)/Tan[ThetaW]^2 + (2*PiVZVP*Tan[ThetaW])/MVZ^2)  },
{dg1, dg2*Tan[ThetaW]+g2*dSinTW/Cos[ThetaW]- dCosTW*g2*Tan[ThetaW]/Cos[ThetaW]}
};
