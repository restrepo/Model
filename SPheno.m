OnlyLowEnergySPheno = True;


MINPAR={{1,lambda1INPUT},
        {2,Lambda2INPUT},
        {3,Lambda3INPUT},
        {4,LamSHINPUT},
        {5,LamSINPUT},
        {4,MSinput},
        {20, vXINPUT}
        };


ParametersToSolveTadpoles = {MuP,mH2};

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
  {lambda1,lambda1INPUT},
  {L2, Lambda2INPUT},
  {L3, Lambda3INPUT},
  {LamSH,LamSHINPUT},
  {LamS,LamSINPUT},
  {MS2, MSInput},
  {vX,vXINPUT}
};


ListDecayParticles = {Fu,Fe,Fd,Ah,hh};
ListDecayParticles3B = {{Fu,"Fu.f90"},{Fe,"Fe.f90"},{Fd,"Fd.f90"}};



DefaultInputValues ={lambda1INPUT -> 0.28};

AddTreeLevelUnitarityLimits=True;
FlagLoopContributions = True;

RenConditionsDecays={
{dCosTW, 1/2*Cos[ThetaW] * (PiVWp/(MVWp^2) - PiVZ/(mVZ^2)) },
{dSinTW, -dCosTW/Tan[ThetaW]},
{dg2, 1/2*g2*(derPiVPheavy0 + PiVPlightMZ/MVZ^2 - (-(PiVWp/MVWp^2) + PiVZ/MVZ^2)/Tan[ThetaW]^2 + (2*PiVZVP*Tan[ThetaW])/MVZ^2)  },
{dg1, dg2*Tan[ThetaW]+g2*dSinTW/Cos[ThetaW]- dCosTW*g2*Tan[ThetaW]/Cos[ThetaW]}
};
