(** Configuration file for shared variables between SARAH model files **)

EvenSingletScalar=True; (* SSDM + Complex even singlet scalar *)
GaugeU1 = True; (* U(1) gauge symmetry *)

(** 
    Anomaly solution: Capital letters may include generations
    {D,i,r,a,c,0,...,m,a,j,o,r,A,n,    A}: 
     D->nDG,                   A->nWG, A -> nMG 
**)
{nDG, nWG, nMG} = {1, 1, 3}; (* number of DG, WG, MG *)
seesaw = True; (* Type-I seesaw mechanism *)
SSDM=False; (* Singlet scalar dark matter *)