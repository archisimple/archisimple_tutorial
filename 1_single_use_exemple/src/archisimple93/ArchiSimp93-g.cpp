/**************************************************************************************/
/*  ARCHI SIMP9.3, d?cembre 2019 (issu de la version 9.1)                             */
/*  R?alis? pour Christophe L afin de prendre en compte les probas de ramif           */
/**************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
//#include <windows.h>   // pour la version windows
#include <sys/time.h>  // pour la version linux
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include "external/tinyxml2/tinyxml2.h"
/*#include <vector>
#include <chrono>
#include <Wincrypt.h>
#include <cstdint>*/


#define NBPASMAX 151 /* Nombre maximal de pas de temps */
#define NBHORMAX 60 /* Nombre maximal d'horizons de sol */
#define MAXLINE 250  /* Longueur maxi de la ligne dans fichiers texte */
#define NBCASEMAX 301  /* Nombre maximal de cases de sol en X, Y et Z */
#define CONVERROR 1000  /* Mximal number of trials allowed when truing to converge in functions */

const std::string filename = "parameter.xml";
const int deltaT=1;          /* Time steps, in days */
const double epsilon=1.0e-8; /* Small value, close to 0 */
const double pi=3.141592653589793238462;  /* close value of PI (close enough) */
const double epaissHor=50.0;  /* Epaisseur horizons de sol (en mm) */
const float longSegNorm=5.0;  /* Longueur habituelle des segments form?s (mm) */
const float longSegMin=2.0;  /* Longueur minimale des segments form?s, quand faible croissance (mm) */
const float dureeSansCreation=1.9; /* Dur?e maximale sans cr?ation de noeud, quand faible croissance (jour) */
const double mailleMin=6.0;  /* Valeur de la maille minimale de sol (mm) */
const double d1=2.0;   /* Premi?re valeur de distance (mm) */
const double d2=4.0;  /* Deuxi?me valeur de distance (mm) */
const double d3=8.0;   /* Troisi?me valeur de distance (mm) */
const double d4=16.0;  /* Quatri?me valeur de distance (mm) */
const double d5=32.0;   /* Cinqui?me valeur de distance (mm) */

typedef float r2[2];  /* Tableau 2D */
typedef float r3[3];  /* Tableau 3D */

typedef struct SysRac *pTSysRac; /* Pour le syst?me racinaire entier */
typedef struct Axe *pTAxe;   /* Pour chacun des axes */
typedef struct Pointe *pTPointe; /* Pour les pointes, parties apicales des racines */
typedef struct Seg *pTSeg; /* Pour les segments, qui sont des portions d'axes */

struct SysRac /* Ensemble d'axes */
  {
  long int nbAxeForm;  /* Nombre d'axes form?s */
  long int nbAxeSup;   /* Nombre d'axes supprim?s */
  long int nbSegForm;    /* Nombre de segments form?s */
  long int nbSeg; /* Nombre de segments tels qu'ils sont compt?s aux 3 dates */
  int nbSem;  /* Nombre de s?minales ?mises */
  int nbAdv;  /* Nombre de racines adventives ?mises */
  float angDep;        /* Orientation */
  r3 origine;           /* Position de l'origine */
  pTAxe premAxe;       /* Premier axe du syst?me (acc?s ?? la liste) */
  pTAxe dernAxe;       /* Dernier axe produit */
  float biomMax[NBPASMAX+1];  /* Biomasse racinaire maximale disponible pendant chaque pas de temps */
  float biomUtilisee[NBPASMAX+1];  /* Biomasse racinaire r?alis?e pendant chaque pas de temps */
  float tSatis[NBPASMAX+1];  /* Taux de satisfaction de la demande ?? chaque pas de temps */
  float longueur; /* Longueur totale de racines */
  float profMax, profMoy; /* Profondeurs maximale et moyenne */
  float distMax, distMoy; /* Distances maximale et moyenne ?? l'axe du syst?me */
  float diamMax; /* Diam?tre maximal, du plus gros segment */
  float xbinf,ybinf,zbinf,xbsup,ybsup,zbsup; /* Bornes en x, y et z */
  float volProd,volPrim,volTot,volTotPrec; /* Volumes racinaires : produit, primaire, total et total au pas pr?c?dent */
  float secPointe; /* Section totale des pointes matures et non s?niles */
  float tSatisMoy; /* Taux de satisfaction moyen */
  float volSolD1, volSolD2, volSolD3, volSolD4, volSolD5; /* Volumes de sol ?? distance d1 ? d5 */
  } ;

struct Pointe /* M?rist?me apical, ou pointe de chaque racine */
  {
  float distPrimInit;  /* Distance de l'apex au dernier primordium initi? */
  float longueur;  /* Longueur non encore exprim?e en allongement de l'axe */
  int dateDerniereCreation; /* Date ?? laquelle il y a eu cr?ation d'un noeud */
  r3 coord;           /* Coordonn?es de la pointe */
  r3 dirCroiss;       /* Direction de croissance */
  r3 dirInit;         /* Direction initiale */
  float age;          /* Age du m?rist?me */
  float diametre;     /* Diam?tre de la pointe */
  unsigned char stop;           /* Stopp?e ?, ou encore en croissance ... */
  unsigned char senile;         /* S?nile ?, ou encore actif ... */
  unsigned char mature;         /* Mature ?, ou encore au stade primordium ... */
  } ;

struct Axe /* Ensemble constitu? d'un m?rist?me et liste de noeuds */
  {
  long int num;      /* Num?ro de l'axe */
  int nbSeg;       /* Nombre de noeuds */
  pTPointe pointe; /* M?rist?me apical */
  pTAxe suivant;     /* Suivant de la liste */
  pTAxe precedent;   /* Pr?c?dent de la liste */
  pTAxe pere;        /* Axe p?re, sur lequel celui-ci est branch? */
  pTSeg premSeg; /* Premier segment de l'axe, sa base */
  pTSeg dernSeg; /* Dernier segment de l'axe, apical */
  } ;

struct Seg
  {
  long int num;      /* Num?ro d'ordre de cr?ation */
  int jourForm;      /* Date de formation (en jours) */
  unsigned char complet; /* Complet (1), c'est-??-dire avec ses deux points, ou non (0) */
  float diametre;    /* Diametre */
  r3 posO;            /* Position dans l'espace de son origine */
  r3 posE;            /* Position dans l'espace de son extr??mit? */
  pTSeg suiv;      /* Suivant dans le prolongement (NULL sinon) */
  pTSeg prec;     /* Pr?c?dent, sur le m?me axe quand non base, sur axe p?re sinon */
  pTAxe axe;        /* Axe auquel appartient le segment */
  unsigned char necrose;       /* Necrose ? 0 : non; 1 : oui */
  } Seg ;

struct Horizon  /* Horizon de sol */
  {
  float croiss;  /* Coefficient de croissance, compris entre 0 et 1 */
  float ramif;   /* Coefficient multiplicateur de distance inter-ramif  */
  float iCMeca;  /* Intensit? de la contrainte m?canique */
  int oCMeca;    /* Orientation de la contrainte m?canique (O iso, ou 1 vert) */
  } ;

typedef Horizon TSol[NBHORMAX];  /* Sol pour la croissance, tableau d'horizons */

/* Fichiers */
FILE *FSeg;      /* Fichier contenant la structure sous forme de segments */
FILE *FPar;      /* Param?tres */
FILE *FSol;      /* Informations sur le sol, par horizons */
FILE *FBiom;      /* Informations sur la biomasse racinaire possible, ?? chaque pas de temps */
// FILE *FAudit;  /* Audit sur le d?roulement de la simulation */
FILE *FSynth;  /* Fichier contenant des variables de synth?se */
// FILE *FVox;    /* Informations sur les voxels colonis?s */

/* Param?tres, lus dans fichier param?tres */
int P_duree=70; /* Dur?e de la simulation, en jours */


// Export
int P_exportType = 1; /* 1 = txt, 2 = RSML */
char P_exportName [200] = "myroot"; /* Name of the exported file */
char P_exportName2 [250] = "myroot"; /* Name of the exported file */
char P_exportExt [50] = ".txt"; /* Name of the extension */

// Caract?risation de l'?mission des racines primaires
float P_vitEmissionSem=0.5; /* Vitesse d'?mission des primaires (en jour-1) */
int P_nbMaxSem=1; /* Nombre maximal de racines primaires */
float P_propDiamSem=1.0; /* Proportion du diam?tre des s?minales (par rapport au diam?tre max) */
float P_angInitMoyVertSem=0.7854; /* Angle d'insertion moyen par rapport ?? la verticale pour les primaires */
float P_angInitETVertSem=0.35;  /* ?cart-type de l'angle d'insertion des primaires */

// Caract?risation de l'?mission des adventives
float P_ageEmissionAdv=12.0; /* ??ge de commencement de l'?mission des racines adventives */
float P_vitEmissionAdv=2.0; /* Vitesse d'?mission des adventives (en jour-1) */
float P_dBaseMaxAdv=30.0; /* Distance ?? la base maximale pour les adventives (mm) */
float P_propDiamAdv=1.0; /* Proportion du diam?tre des adventives (par rapport aux diam?tre max) */
int P_nbMaxAdv=40; /* Nombre maximal de racines adventives */

float P_angInitMoyVertAdv=1.4; /* Angle d'insertion moyen par rapport ?? la verticale pour les adventives */
float P_angInitETVertAdv=0.7;  /* ?cart-type de l'angle d'insertion des adventives */

// Croissance radiale
float P_coeffCroissRad=0.6; // coefficient de croissance radiale

// Allongement (croissance axiale)
float P_diamMin=0.10;  /* Diam?tre minimal en deça duquel il n'y a pas de croissance (mm) */
float P_diamMax=1.1;   /* Diam?tre maximal donn? aux racines primaires (mm) */
float P_penteVitDiam=12.0; /* pente de la relation entre vitesse de croissance et diam?tre (mm.mm.jour-1) */
int P_tendanceDirTropisme=2;  /* Type de tropisme (0: plagio; -1: geo-; +1: geo+; 2: exo */
float P_intensiteTropisme=0.2; /* Coefficient multipli? par le diam?tre pour chaque racine */
float P_penteDureeCroissDiam2=3000.0; /* pente de la relation dur?e de croissance versus diam?tre^2 */

// Ramification
float P_ageMaturitePointe=4.5;  /* ?ge de maturit? des m?rist?mes (jours) */
float P_distRamif=4.0; /* distance inter-ramification (mm) */
float P_probEmergeDmax=0.8;  /* probabilit? d'?mergence d'une lat?rale sur un axe de diam Dmax */
float P_probEmergeDmin=0.0;  /* probabilit? d'?mergence d'une lat?rale sur un axe de diam Dmin */
float P_propDiamRamif=0.2; /* proportion de diam?tre des filles par rapport à leur m?re */
float P_coeffVarDiamRamif=0.30; /* coefficient de variation du diam?tre des ramifs */
float P_angLat=1.3; /* angle d'insertion des racines lat?rales */

// Mortalit?
float P_TMD=0.2; /* Tissue mass density, ou masse volumique */
float P_penteDureeVieDiamTMD=2000.0; /* pente de la relation dur?e de vie versus diam?tre et TMD */

/* Variables globales diverses */
int temps=0;  /* Le temps, en jours */
r3 orig;      /* Position d'origine du syst?me racinaire */

unsigned char vox[NBCASEMAX+1][NBCASEMAX+1][NBCASEMAX+1];  /* tableau sol-voxel dont les cases vont contenir des entiers entre 1 et 6 */

float maille=mailleMin; /* Valeur initialis?e de la maille de sol */
double volElemSol;  /* Volume ?l?mentaire de sol associ? à la maille (mm3) */

pTSysRac sR;  /* Le syst?me racinaire */
TSol sol;     /* Le sol */



/****************************************************************************/
/****************************************************************************/
double dRandUnif(void)
  /* Cette fonction tire un aléatoire uniforme réel entre 0 et 1 */
{
  double tirage;
  tirage=(double) rand()/(double) RAND_MAX;
  if (tirage<epsilon) { tirage=epsilon; }
  return(tirage);
}


/****************************************************************************/
/****************************************************************************/
double drandUnifEntre(double min, double max)
{
  return ((double) rand()/(double) RAND_MAX) * (max-min) + min;
}


/****************************************************************************/
/****************************************************************************/
void norme(r3 u, r3 un)
/* Cette fonction norme le vecteur u de l'espace de dimension 3.
  Le vecteur norme de retour est appele un. */
{
double norU;
  norU=sqrt((u[0]*u[0])+(u[1]*u[1])+(u[2]*u[2]));
  if (norU<epsilon)
  {
  printf("WARNING, null vector. Norm is : %f \n",norU);
  exit(1);
  }
  else
  {
   un[0]=u[0]/norU;
   un[1]=u[1]/norU;
   un[2]=u[2]/norU;
  }
}  /* Fonction Norme */


/****************************************************************************/
/****************************************************************************/
double prodScal(r3 u,r3 v)
/* Cette fonction retourne le produit scalaire de 2 vecteurs u et v de
  l'espace a 3 dimensions. */
{
double prodScal;
  prodScal=(u[0]*v[0])+(u[1]*v[1])+(u[2]*v[2]);
  return(prodScal);
}  /* Fonction prodScal */


/****************************************************************************/
/****************************************************************************/
void prodVect(r3 u, r3 v, r3 u_vect_v)
/* Cette fonction calcule le produit vectoriel de deux vecteurs u et v
  de l'espace de dimension 3. Le vecteur resultant est u_vect_v. */
{
  u_vect_v[0]=(u[1]*v[2])-(v[1]*u[2]);
  u_vect_v[1]=(u[2]*v[0])-(v[2]*u[0]);
  u_vect_v[2]=(u[0]*v[1])-(v[0]*u[1]);
}   /* Fonction prodVect */


/****************************************************************************/
/****************************************************************************/
void rotVect(double omega, r3 u, r3 x, r3 rot_x)

/* Cette fonction calcule le vecteur rot_x dans l'espace de dimension 3,
  issu de la rotation du vecteur x autour d'un axe dont u est un vecteur
  unitaire. La rotation se fait d'un angle omega radians. Elle appelle
  PRODSCAL, PRODVECT. */
{
double uscalx;   /* produit scalaire u.x  */
r3    uvectx;   /* produit vectoriel u^x */

  uscalx=prodScal(u,x);
  prodVect(u,x,uvectx);

  rot_x[0]=((1-cos(omega))*uscalx*u[0])
      +(cos(omega)*x[0])+(sin(omega)*uvectx[0]);
  rot_x[1]=((1-cos(omega))*uscalx*u[1])
      +(cos(omega)*x[1])+(sin(omega)*uvectx[1]);
  rot_x[2]=((1-cos(omega))*uscalx*u[2])
      +(cos(omega)*x[2])+(sin(omega)*uvectx[2]);

}  /* Fonction rotVect */


/****************************************************************************/
/****************************************************************************/
void rotZ(r3 u, r3 v, double teta)
/* Cette fonction fait tourner "u" d'un angle "teta" autour de l'axe (Oz);
  le vecteur calcule est "v" */
{
  v[0]=(u[0]*cos(teta))-(u[1]*sin(teta));
  v[1]=(u[0]*sin(teta))+(u[1]*cos(teta));
  v[2]=u[2];
}


/****************************************************************************/
/****************************************************************************/
int iRandUnif(int imax)

/* Cette fonction tire un al?atoire uniforme entier entre 0 et imax */
{
  int tirage;

  tirage=imax+1;
  while (tirage>imax) tirage=rand();
  return tirage;
}


/****************************************************************************/
/****************************************************************************/
void ouvreFichiers(void)
/* Cette fonction ouvre les fichiers, en lecture et ?criture */
{

  if(P_exportType == 2){
    sprintf(P_exportExt, ".rsml");
  }
  sprintf(P_exportName2, "%s%s", P_exportName, P_exportExt);

  FSeg=fopen(P_exportName2,"w");    // Fichier contenant les segments

//  FAudit=fopen("audit.txt","w");
  FPar=fopen("paramarch93.txt","rt");   // Param?tres de simulation

  FSol=fopen("sol.txt","rt");
  FBiom=fopen("biomrac.txt","rt");   // Biomasse racinaire autoris?e ? chaque pas
//  FSynth=fopen("synth.txt","w");
//  FVox=fopen("vox.txt","w");   //  Fichier des voxels
} /* Fonction ouvreFichiers */


/****************************************************************************/
/****************************************************************************/
void litSol(void)
/* Fonction de lecture des caract?ristiques du sol, une ligne par horizon */
{
int hor;              /* Compteur des horizons */
char bid[MAXLINE];    /* Chaîne qui accueille les caract?res suppl?mentaires */

fgets(bid,MAXLINE-1,FSol);          /* Ligne entête */
for (hor=0; hor<NBHORMAX; hor++)
{
  fscanf(FSol,"%f %f %f %i",&sol[hor].croiss,&sol[hor].ramif,&sol[hor].iCMeca,&sol[hor].oCMeca);
//  fscanf(FSol,"%f",&sol[hor].croiss); // Favorable à la croissance
//  fscanf(FSol,"%f",&sol[hor].ramif);  // Favorable à la ramification
//  fscanf(FSol,"%f",&sol[hor].iCMeca); // Intensit? de la contrainte
//  fscanf(FSol,"%d",&sol[hor].oCMeca); // Orientation 0: iso, 1: verticale
  fgets(bid,MAXLINE-1,FSol);
}

} /* Fonction litSol */


/****************************************************************************/
/****************************************************************************/
void litBiomasseMaxSR(pTSysRac sR)
/* Fonction de lecture des biomasses maximales à chaque pas de temps */
{
int pas;              /* Compteur des pas de temps */
char bid[MAXLINE];    /* Chaîne qui accueille les caract?res suppl?mentaires */

fgets(bid,MAXLINE-1,FBiom);          /* Ligne entête */
for (pas=1; pas<NBPASMAX; pas++)
{
  fscanf(FBiom,"%f",&(sR->biomMax[pas]));
  fgets(bid,MAXLINE-1,FBiom);
}
  sR->biomMax[0]=sR->biomMax[1];

} /* Fonction litBiomasseMaxSR */


/****************************************************************************/
/****************************************************************************/
double croissSol(TSol sol, double profondeur)
/* Renvoie le coefficient de croissance du sol à la Profondeur donn?e */
{
int hor;

  hor=(int) floor(profondeur/epaissHor);
  if (hor>=NBHORMAX) hor=NBHORMAX-1;
  if (hor<0) hor=0;

  return(sol[hor].croiss);
} /* Fonction croissSol */


/****************************************************************************/
/****************************************************************************/
double ramifSol(TSol sol, double profondeur)
/* Renvoie le coefficient de ramification du sol à la profondeur donn?e */
{
int hor;

  hor=(int) floor(profondeur/epaissHor);
  if (hor>=NBHORMAX) hor=NBHORMAX-1;
  if (hor<0) hor=0;

  return(sol[hor].ramif);
} /* Fonction ramifSol */


/****************************************************************************/
/****************************************************************************/
double iCMecaSol(TSol sol, double profondeur)
/* Renvoie l'intensit? de la contrainte m?ca du sol à la Profondeur donn?e */
{
int hor;

  hor=(int) floor(profondeur/epaissHor);
  if (hor>=NBHORMAX) hor=NBHORMAX-1;
  if (hor<0) hor=0;

  return(sol[hor].iCMeca);
} /* Fonction iCMecaSol */


/****************************************************************************/
/****************************************************************************/
int oCMecaSol(TSol sol, double profondeur)
/* Renvoie l'indice de la direction de contrainte : 0 pour iso, 1 pour verti */
{
int hor;

  hor=(int) floor(profondeur/epaissHor);
  if (hor>=NBHORMAX) hor=NBHORMAX-1;
  if (hor<0) hor=0;

  return(sol[hor].oCMeca);
} /* Fonction oCMecaSol */


/****************************************************************************/
/****************************************************************************/
double tireGaussien(float moy, float et)
{  /* R?alise un tirage gaussien dans une distribution de moyenne moy et ?cart-type et */
  double tireGaussien,tire1,tire2;

  tire1=dRandUnif();
  tire2=dRandUnif();
  tireGaussien=moy+(et*sqrt(-log(tire1))*cos(pi*tire2)*1.414);
  return(tireGaussien);
} /* Fonction tireGaussien */


/****************************************************************************/
/****************************************************************************/
double tireAngRad(void)
{   /* Tire l'angle radial dans l'intervalle 0 - 2*Pi */

return (2.0*pi*dRandUnif());
} /* Fonction TireAngRad */


/****************************************************************************/
/****************************************************************************/
void increNbSegSR(pTSysRac sR)
/* Incr?mente le nombre de noeuds qui a ?t? form? dans ce syst?me sR */
{
  sR->nbSegForm++;
} /* Fonction increNbSegSR */


/****************************************************************************/
/****************************************************************************/
pTSeg creeSeg(void)
/* Cette fonction retourne une nouvelle variable de type pTSeg,
  c'est-à-dire un pointeur sur le type Seg */
{
pTSeg seg;
  seg=(pTSeg) malloc(sizeof(Seg));
  if (seg==NULL)
  { printf("Memory issue allocation dans creeSeg \n"); exit(1); }

return seg;
} /* Fonction creeSeg */


/****************************************************************************/
/****************************************************************************/
pTSeg initialiseSeg(long int num, r3 posOrig, r3 posExtrem, double diam, pTAxe axeSeg, unsigned char comp, pTSeg precedent)
/* Cette fonction retourne une nouvelle variable de type pTSeg,
  dont une partie des valeurs est initialis?e */
{
pTSeg seg;

  seg=creeSeg();

  seg->num=num;
  seg->jourForm=temps;
  seg->necrose=0;
  seg->complet=comp;

  seg->diametre=diam;
  seg->axe=axeSeg;

  seg->posO[0]=posOrig[0];
  seg->posO[1]=posOrig[1];
  seg->posO[2]=posOrig[2];

  seg->posE[0]=posExtrem[0];
  seg->posE[1]=posExtrem[1];
  seg->posE[2]=posExtrem[2];

  seg->suiv=NULL;  // pour l'instant
  seg->prec=precedent;

return seg;
} /* Fonction initialiseSeg */


/****************************************************************************/
/****************************************************************************/
void detruitSeg(pTSeg segADetruire)
/* Supprime un noeud en m?moire */
{
  free(segADetruire);
} /* Fonction detruitSeg */


/****************************************************************************/
/****************************************************************************/
pTPointe creePointe(void)
/* Cette fonction retourne une nouvelle variable de type pTPointe,
  c'est-à-dire un pointeur sur le type TPointe */
{
pTPointe pointe;
  pointe=(pTPointe) malloc(sizeof(Pointe));
  if (pointe==NULL)
  { printf("Memory issue allocation dans creePointe \n"); exit(1); }

return pointe;
} /* Fonction creePointe */


/****************************************************************************/
/****************************************************************************/
pTPointe initialisePointe(float diam, r3 position, r3 direction)
/* Cette fonction retourne une nouvelle variable de type pTPointe,
  dont les valeurs sont en partie initialis?es */
{
pTPointe pointe;

  pointe=creePointe();

  pointe->distPrimInit=0.0;
  pointe->longueur=0.0;
  pointe->age=0.0;
  pointe->diametre=diam;
  pointe->stop=0;
  pointe->senile=0;
  pointe->mature=0;

  pointe->coord[0]=position[0];
  pointe->coord[1]=position[1];
  pointe->coord[2]=position[2];

  pointe->dirCroiss[0]=direction[0];
  pointe->dirCroiss[1]=direction[1];
  pointe->dirCroiss[2]=direction[2];

  pointe->dirInit[0]=direction[0];
  pointe->dirInit[1]=direction[1];
  pointe->dirInit[2]=direction[2];


return pointe;
} /* Fonction initialisePointe */


/****************************************************************************/
/****************************************************************************/
void deflecMecaPointe(pTPointe pointe, r3 dirApresMeca, double elong)
{
const double teta=15.0; /* Angle autour de G, en degres */

r3 vTire,vTireN,dirInt;
double profondeur, cont;

  profondeur=pointe->coord[2];
  cont=iCMecaSol(sol,profondeur);
  if (oCMecaSol(sol,profondeur)==1)  /* Contrainte anisotrope verticale */
  {
  /* Tirage vecteur dans l'angle Teta autour de G */
  do
  {
  vTire[0]=(2.0*dRandUnif()-1.0)*sin(pi*teta/180.0);
  vTire[1]=(2.0*dRandUnif()-1.0)*sin(pi*teta/180.0);
  do { vTire[2]=dRandUnif(); } while (vTire[2]>cos(pi*teta/180.0));
  norme(vTire,vTireN);
  }  while (vTireN[2]>cos(pi*teta/180.0));
  dirInt[0]=pointe->dirCroiss[0]+(elong*vTireN[0]*cont);
  dirInt[1]=pointe->dirCroiss[1]+(elong*vTireN[1]*cont);
  dirInt[2]=pointe->dirCroiss[2]+(elong*vTireN[2]*cont);
  }
  else    /* Contrainte isotrope [oCMecaSol(Profondeur)==0] */
    {
    vTire[0]=2.0*dRandUnif()-1.0;
    vTire[1]=2.0*dRandUnif()-1.0;
    vTire[2]=2.0*dRandUnif()-1.0;
    norme(vTire,vTireN);
  if (prodScal(vTireN,pointe->dirCroiss)<0.0)
    {
    vTireN[0]=-vTireN[0];
    vTireN[1]=-vTireN[1];
    vTireN[2]=-vTireN[2];
  }
  dirInt[0]=pointe->dirCroiss[0]+(elong*vTireN[0]*cont);
  dirInt[1]=pointe->dirCroiss[1]+(elong*vTireN[1]*cont);
  dirInt[2]=pointe->dirCroiss[2]+(elong*vTireN[2]*cont);
  }
  norme(dirInt,dirApresMeca);

} /* Fonction deflecMecaPointe */


/****************************************************************************/
/****************************************************************************/
void deflecGeoPointe(pTPointe pointe, r3 dirApresMeca, r3 dirApresGeo, double elong)
/* Version avec plagiotropisme */
{
r3 dirInt,vGeoInt,vGeo;

  switch (P_tendanceDirTropisme) {
    case -1 : vGeo[0]=0.0;                  /* Gravitropisme n?gatif */
              vGeo[1]=0.0;
              vGeo[2]=-1.0;
              break;
    case 0 : vGeoInt[0]=pointe->dirInit[0]; /* Plagiotropisme */
             vGeoInt[1]=pointe->dirInit[1];
             vGeoInt[2]=0.0;
             norme(vGeoInt,vGeo);
             break;
    case 1 : vGeo[0]=0.0;                  /* Gravitropisme positif */
             vGeo[1]=0.0;
             vGeo[2]=1.0;
              break;
    case 2 : vGeoInt[0]=pointe->dirInit[0]; /* Exotropisme */
             vGeoInt[1]=pointe->dirInit[1];
             vGeoInt[2]=pointe->dirInit[2];
             norme(vGeoInt,vGeo);
             break;
    default : vGeo[0]=0.0;                 /* Gravitropisme positif */
              vGeo[1]=0.0;
              vGeo[2]=1.0;
              break;
  }

  dirInt[0]=dirApresMeca[0]+(vGeo[0]*P_intensiteTropisme*elong*pointe->diametre);
  dirInt[1]=dirApresMeca[1]+(vGeo[1]*P_intensiteTropisme*elong*pointe->diametre);
  dirInt[2]=dirApresMeca[2]+(vGeo[2]*P_intensiteTropisme*elong*pointe->diametre);

  norme(dirInt,dirApresGeo);
} /* Fonction deflecGeoPointe */


/****************************************************************************/
/****************************************************************************/
void deflecSurfPointe(pTPointe Pointe, r3 dirApresGeo, r3 dirApresSurf)
{
const double profLim=50.0*dRandUnif();
r3 dirInt;
  dirInt[0]=dirApresGeo[0];
  dirInt[1]=dirApresGeo[1];
  dirInt[2]=dirApresGeo[2];

  if ((dirInt[2]<0.0) && ((Pointe->coord[2])<profLim)) dirInt[2]=dirInt[2]/10.0;
  norme(dirInt,dirApresSurf);
} /* Fonction deflecSurfPointe */


/****************************************************************************/
/****************************************************************************/
void reorientePointe(pTPointe pointe, double elong)
{
r3 dirInt1, dirInt2, nouvDir;

  deflecMecaPointe(pointe,dirInt1,elong);
  deflecGeoPointe(pointe,dirInt1,dirInt2,elong);
  deflecSurfPointe(pointe,dirInt2,nouvDir);

  pointe->dirCroiss[0]=nouvDir[0];
  pointe->dirCroiss[1]=nouvDir[1];
  pointe->dirCroiss[2]=nouvDir[2];


} /* Fonction reorientePointe */


/****************************************************************************/
/****************************************************************************/
double calcElongationPointe(pTPointe pointe, TSol sol)
/* Calcul de l'?longation potentielle affect?e par le sol */
{
  if ((pointe->mature) && (!pointe->stop) && (!pointe->senile) && (pointe->diametre>P_diamMin))
//    return (pointe->diametre-P_diamMin)*deltaT*P_penteVitDiam*croissSol(sol,pointe->coord[2]);
    return (pointe->diametre)*deltaT*P_penteVitDiam*croissSol(sol,pointe->coord[2]);
    
  else return 0.0;

} /* Fonction calcElongationPointe */


/****************************************************************************/
/****************************************************************************/
void developpePointe(pTPointe pointe)
{ /* Assure l'?volution de la pointe, en la faisant vieillir et en changeant ses variables d'?tat au cours de sa vie */
 
  pointe->age=pointe->age+deltaT; /* Incr?mente l'?ge du m?rist?me selon le pas de temps */ 
  
  if ((!pointe->mature)&&(pointe->age>P_ageMaturitePointe))
  {
    pointe->mature=1;  /* Le primordium devient m?rist?me vrai */
    pointe->age=0.0;   /* Son ?ge est r?initialis? à 0 en tant que pointe mature */
  }
  
  if ((pointe->mature)&&(!pointe->stop)&&(pointe->age>(P_penteDureeCroissDiam2*pointe->diametre*pointe->diametre)))
  {
    pointe->stop=1;  /* La pointe stoppe sa croissance */
  }
  
  if ((pointe->mature)&&(pointe->stop)&&(!pointe->senile)&&(pointe->age>(P_penteDureeVieDiamTMD*pointe->diametre*P_TMD)))
  {
    pointe->senile=1;  /* La pointe devient s?nile */
  }
  
} /* Fonction developpePointe */


/****************************************************************************/
/****************************************************************************/
void deplacePointe(pTPointe pointe, double elong)
{ /* Assure le d?placement du m?rist?me suite à croissance axiale */

  /* Sa position est modifi?e */
  pointe->coord[0]=pointe->coord[0]+(elong*pointe->dirCroiss[0]);
  pointe->coord[1]=pointe->coord[1]+(elong*pointe->dirCroiss[1]);
  pointe->coord[2]=pointe->coord[2]+(elong*pointe->dirCroiss[2]);

  /* Son attribut distPrimInit est modifi? */
  pointe->distPrimInit+=elong;

} /* Fonction deplacePointe */


/****************************************************************************/
/****************************************************************************/
double distInterRamifPointe(pTPointe pointe, TSol sol)
{ /* Renvoie la valeur locale de la distance inter-ramification de la pointe */

  return (P_distRamif*ramifSol(sol,pointe->coord[2]));

} /* Fonction distInterRamifPointe */


/****************************************************************************/
/****************************************************************************/
void detruitPointe(pTPointe pointeADetruire)
/* Supprime une pointe */
{
  free(pointeADetruire);
} /* Fonction detruitPointe */


/****************************************************************************/
/****************************************************************************/
pTAxe creeAxe(void)
/* Cette fonction retourne une nouvelle variable de type pTAxe,
  c'est-à-dire un pointeur sur le type Axe */
{
pTAxe axe;
  axe=(pTAxe) malloc(sizeof(Axe));
  if (axe==NULL)
  { printf("Memory issue allocation dans creeAxe \n"); exit(1); }

return axe;
} /* Fonction creeAxe */


/****************************************************************************/
/****************************************************************************/
pTAxe initialiseAxe(long int numAxe, float diamPointe, r3 origine, r3 dirInit, pTAxe axePere, pTSeg segPorteur)
/* Cette fonction retourne une nouvelle variable de type pTAxe,
  c'est-à-dire un pointeur sur le type Axe */
{
  pTAxe nouvAxe;
  pTSeg premierSeg;

  nouvAxe=creeAxe();
  premierSeg=initialiseSeg(sR->nbSegForm+1,origine,origine,diamPointe,nouvAxe,0,segPorteur);
  nouvAxe->pointe=initialisePointe(diamPointe,origine,dirInit);
  nouvAxe->premSeg=premierSeg;
  nouvAxe->dernSeg=premierSeg;
  nouvAxe->nbSeg=1;
  nouvAxe->num=numAxe;
  nouvAxe->pere=axePere;

  nouvAxe->suivant=NULL;
  nouvAxe->precedent=NULL;

  return nouvAxe;
} /* Fonction initialiseAxe */


/****************************************************************************/
/****************************************************************************/
void ajouteSegProlongeAxe(pTAxe axe, pTSeg segAAjouter)
/* Cette fonction ajoute un segment de prolongement en position apicale
à l'axe concern?, et incr?mente son compteur de segments */
{
pTSeg ancienSegTerm;

  ancienSegTerm=axe->dernSeg;

  // Si ce dernier segment est complet, il faut prolonger la liste
  if (ancienSegTerm->complet) {
    ancienSegTerm->suiv=segAAjouter;
    segAAjouter->prec=ancienSegTerm;
    axe->dernSeg=segAAjouter;
    axe->nbSeg++;
  }
  // Sinon, il faut juste compl?ter le dernier segment
  else {
    // rien à faire
    // les mises à jour sont faites dans developpeAxeSR
    // on ne doit pas passer ici
  }

} /* Fonction ajouteSegProlongeAxe */


/****************************************************************************/
/****************************************************************************/
void ajouteAxeSR(pTSysRac sR, pTAxe axeAAjouter)
/* Cette fonction ins?re un axe dans la chaîne des axes du syst?me racinaire,
elle incr?mente en même temps le compteur d'axes et de segments */
{
  if (sR->premAxe==NULL)  /* Le syst?me racinaire est vide */
  {
    axeAAjouter->suivant=NULL;
    axeAAjouter->precedent=NULL;
    sR->premAxe=axeAAjouter;
    sR->dernAxe=axeAAjouter;
  }
  else /* Le syst?me contient d?jà des axes, assure le chaînage double des axes */
  {
    axeAAjouter->suivant=NULL;
    axeAAjouter->precedent=sR->dernAxe;
    sR->dernAxe->suivant=axeAAjouter;
    sR->dernAxe=axeAAjouter;
  }
  sR->nbAxeForm++;
  sR->nbSegForm++;   // à chaque axe, un segment

} /* Fonction ajouteAxeSR */


/****************************************************************************/
/****************************************************************************/
int latEmerge(float diamAxePere)
{   /* Renvoie 1 ou 0 suivant que la racine lat?rale ?merge ou pas */
    /* Cela d?pend du diam?tre de l'axe p?re, et des deux param?tres P_probEmergeDmin et P_probEmergeDmax */
    
  /* On commence par calculer la probabilit? d'?mergence pour un axe qui a pour diam?tre diamAxePere */
  float probaEmerge = P_probEmergeDmin + ((diamAxePere-P_diamMin)*(P_probEmergeDmax-P_probEmergeDmin)/(P_diamMax-P_diamMin));
  
  /* On compare un tirage al?atoire uniforme entre 0 et 1 ? la valeur de cette probabilit? et on retourne la comparaison */
  return(probaEmerge > dRandUnif());

} /* Fonction latEmerge */


/****************************************************************************/
/****************************************************************************/
int axeRamifiable(pTAxe axe)
{   /* Renvoie 1 ou 0 suivant que l'axe est ramifiable ou non */
    // 3 conditions doivent ?tre r?unies (diam?tre, position, ?ge)

  return((axe->pointe->diametre > 1.1*P_diamMin) && (axe->pointe->distPrimInit > P_distRamif)  && (((P_penteDureeCroissDiam2*axe->pointe->diametre*axe->pointe->diametre) - axe->pointe->age ) > P_ageMaturitePointe));

} /* Fonction axeRamifiable */


/****************************************************************************/
/****************************************************************************/
float tireDiamPointeFille(pTAxe axePere)
{   /* Tire le diam?tre d'un m?rist?me de ramification suivant celui du p?re
       pour la ramification s?quentielle */

	float moy=(axePere->pointe->diametre*P_propDiamRamif) + (P_diamMin*(1.0-P_propDiamRamif));
	float et=moy*P_coeffVarDiamRamif;
	float diamPFille=10000.0;  // initialisation à une forte valeur pour boucle de tirage
	int count = 1;
  while (diamPFille>(0.95*axePere->pointe->diametre)){
    count = count + 1;
    if(count >= CONVERROR){
      printf("---------- \n \n");
      printf("Error, P_propDiamRamif is too large, diameter cannot converge after %3i tries \n", CONVERROR);
      printf("-> Change the value of P_propDiamRamif in the input file  \n");
      printf("---------- \n \n");
      exit(1);
    }
    diamPFille=tireGaussien(moy,et);
  }

	return diamPFille;  /* La racine fille a un diam?tre inf?rieur ? celui de sa m?re */

} /* Fonction tireDiamPointeFille */


/****************************************************************************/
/****************************************************************************/
float tireDiamSem(void)
{   /* Tire le diam?tre d'une pointe de s?minale, avec une certaine variation */
// Pas utilis? dans toutes les versions
  double dmin=0.90*P_diamMax;
  double dmax=1.0*P_diamMax;
  return drandUnifEntre(dmin,dmax);
} /* Fonction tireDiamSem */
/****************************************************************************/
/****************************************************************************/
void origineAdv(pTSeg segPere, r3 origineFils)
{   /* Calcule la position du point d'origine d'une tardive sur le seg p?re */

  double rel=dRandUnif();  /* definira la position relative sur le segment */
  origineFils[0]=(rel*segPere->posO[0]) + ((1.0-rel)*segPere->posE[0]);
  origineFils[1]=(rel*segPere->posO[1]) + ((1.0-rel)*segPere->posE[1]);
  origineFils[2]=(rel*segPere->posO[2]) + ((1.0-rel)*segPere->posE[2]);

} /* Fonction origineAdv */
/****************************************************************************/
/****************************************************************************/
void origineRamif(pTAxe axePere, r3 origineFils)
{   /* Calcule la position du point d'origine d'une ramification */
origineFils[0]=axePere->pointe->coord[0]-
                  (axePere->pointe->distPrimInit*axePere->pointe->dirCroiss[0]);
origineFils[1]=axePere->pointe->coord[1]-
                  (axePere->pointe->distPrimInit*axePere->pointe->dirCroiss[1]);
origineFils[2]=axePere->pointe->coord[2]-
                  (axePere->pointe->distPrimInit*axePere->pointe->dirCroiss[2]);
} /* Fonction origineRamif */
/****************************************************************************/
/****************************************************************************/
void orienteRamif(pTAxe axePere, r3 dirFils)
{   /* Calcule la direction d'un axe fils issu de ramification */
r3 vAxeRot,rotDirCroiss;
double norVProjHor,angRot;

/* Calcul de la norme de la projection direction sur plan horizontal */
norVProjHor=sqrt((axePere->pointe->dirCroiss[0]*axePere->pointe->dirCroiss[0])+
                 (axePere->pointe->dirCroiss[1]*axePere->pointe->dirCroiss[1]));
if (norVProjHor<epsilon)
{
  vAxeRot[0]=1.0; /* Vecteur initial vertical */
  vAxeRot[1]=0.0;
  vAxeRot[2]=0.0; /* Vecteur (1,0,0) choisi pour axe de rotation */
}
else
{
  vAxeRot[0]=axePere->pointe->dirCroiss[1]/norVProjHor;
  vAxeRot[1]=-axePere->pointe->dirCroiss[0]/norVProjHor;
  vAxeRot[2]=0.0;
}
/* On fait tourner dirCroiss autour de vAxeRot d'un angle d'insertion */
angRot=P_angLat;
rotVect(angRot,vAxeRot,axePere->pointe->dirCroiss,rotDirCroiss);

/* On fait tourner rotDirCroiss autour de dirCroiss d'un angle radial */
angRot=tireAngRad();
rotVect(angRot,axePere->pointe->dirCroiss,rotDirCroiss,dirFils);
} /* Fonction orienteRamif */
/****************************************************************************/
/****************************************************************************/
void ramifieAxe(pTAxe axePere)
{
pTAxe nouvAxe;
float diamRamif;
r3 origRamif, dirRamif;

  /* D?cr?mente la distance au dernier primordium initi? */
  axePere->pointe->distPrimInit-=distInterRamifPointe(axePere->pointe,sol);

  /* Calcul des attributs d'une ramification */
  diamRamif=tireDiamPointeFille(axePere);    /* Tire le diam?tre de sa pointe */

  if ((diamRamif > P_diamMin) && latEmerge(axePere->pointe->diametre))  /* Si la racine est assez grosse et si elle a la chance d'?merger (? modifier) */
  {
    origineRamif(axePere,origRamif);         /* Calcule sa position */
    orienteRamif(axePere,dirRamif);          /* Calcule sa direction */

    nouvAxe=initialiseAxe(sR->nbAxeForm+1,diamRamif,origRamif,dirRamif,axePere,axePere->dernSeg);

    ajouteAxeSR(sR,nouvAxe);

  }

} /* Fonction ramifieAxe */
/****************************************************************************/
/****************************************************************************/
void developpeAxe(pTAxe axe,float taux)
/* Assure le d?veloppement de l'axe, avec diff?rentes composantes */
{
double elongation;
pTSeg nouvSeg;

  elongation=taux*calcElongationPointe(axe->pointe,sol);

  axe->pointe->longueur+=elongation;

  while (axe->pointe->longueur > longSegNorm) { // on fait un segment "normal"

    axe->pointe->dateDerniereCreation=temps;

    axe->pointe->longueur-=longSegNorm;

    /* Calcule et affecte la nouvelle direction de croissance du m?rist?me */
    reorientePointe(axe->pointe,longSegNorm);

    /* Le m?rist?me se d?place */
    deplacePointe(axe->pointe,longSegNorm);

    if (axe->dernSeg->complet) {
      /* Il g?n?re un nouveau segment sur cet axe à sa nouvelle position */
      increNbSegSR(sR);
      nouvSeg=initialiseSeg(sR->nbSegForm,axe->dernSeg->posE,axe->pointe->coord,axe->pointe->diametre,axe,1,axe->dernSeg);
      ajouteSegProlongeAxe(axe,nouvSeg);
    }
    else { // le premier segment est incomplet, on le modifie
      axe->dernSeg->complet=1;
      axe->dernSeg->posE[0]=axe->pointe->coord[0];
      axe->dernSeg->posE[1]=axe->pointe->coord[1];
      axe->dernSeg->posE[2]=axe->pointe->coord[2];
      axe->dernSeg->jourForm=temps;
    }
//    printf(" Dans developpeAxe\n");

    while (axeRamifiable(axe)) ramifieAxe(axe); // on ramifie ?ventuellement

  } // fin du while  (axe->pointe->longueur > longSegNorm)

  if (((temps - axe->pointe->dateDerniereCreation) > dureeSansCreation)&&(axe->pointe->longueur > longSegMin)) { /* production segment court  */

    axe->pointe->dateDerniereCreation=temps;


    /* Calcule et affecte la nouvelle direction de croissance de la pointe */
    reorientePointe(axe->pointe,axe->pointe->longueur);

    /* La pointe se d?place */
    deplacePointe(axe->pointe,axe->pointe->longueur);

    /* Elle g?n?re un nouveau segment sur cet axe à sa nouvelle position */
    if (axe->dernSeg->complet) {
      /* Il g?n?re un nouveau segment sur cet axe à sa nouvelle position */
      increNbSegSR(sR);
      nouvSeg=initialiseSeg(sR->nbSegForm,axe->dernSeg->posE,axe->pointe->coord,axe->pointe->diametre,axe,1,axe->dernSeg);
      ajouteSegProlongeAxe(axe,nouvSeg);
    }
    else { // le premier segment est incomplet, on le modifie
      axe->dernSeg->complet=1;
      axe->dernSeg->posE[0]=axe->pointe->coord[0];
      axe->dernSeg->posE[1]=axe->pointe->coord[1];
      axe->dernSeg->posE[2]=axe->pointe->coord[2];
      axe->dernSeg->jourForm=temps;
    }

    axe->pointe->longueur=0.0; // remet la longueur en attente du m?rist?me à 0

    while (axeRamifiable(axe)) ramifieAxe(axe); // on ramifie ?ventuellement

  } // fin du if (production d'un segment court)


} /* Fonction developpeAxe */
/****************************************************************************/
/****************************************************************************/
void calcTSatisMoySR(pTSysRac sR)
{
/* Calcul du taux de satisfaction moyen sur la p?riode ?coul?e (de 1 à temps) */

  double tSatisCum=0.0;

  for (int date=1; date<=temps; date++) /* Boucle sur la p?riode ?coul?e */
  {
    tSatisCum+=sR->tSatis[date];
  }
  sR->tSatisMoy=tSatisCum/temps;

}  /* Fonction calcTSatisMoySR */
/****************************************************************************/
/****************************************************************************/
float calcTauxSatis(pTSysRac sR)
{
  float taux1, taux2;
  float biomU0,biomU1,biomU2;  // biomasses utilis?es aux diff?rents pas (-2, -1, actuel)
	float biomD1,biomD2;  // biomasses disponibles aux pas pr?c?dent et actuel

  calcTSatisMoySR(sR);
  
  if (sR->biomUtilisee[temps-1] < 0.0002) return 1.0;   // en d?but de croissance notamment, on ne g?re pas la limitation ?ventuelle

  else {
  	biomU0=0.0002;
  	if (temps>1) biomU0=sR->biomUtilisee[temps-2];   // biomasse utilis?e au pas de temps - 2, s'il existe
  	biomU1=sR->biomUtilisee[temps-1];    // biomasse utilis?e au pas de temps - 1
  	biomU2=biomU1*(1+((biomU1-biomU0)/(0.5*(biomU0+biomU1))));   // biomasse au pas de temps t (non connue, c'est juste une estimation) 
  	
  	biomD1=sR->biomMax[temps-1]; // biomasse disponible au pas de temps -1
  	biomD2=sR->biomMax[temps];   // biomasse disponible au pas de temps actuel
  	
  	taux1=biomD1/biomU1;
    if (taux1>1.0) taux1=1.0;
  	
  	taux2=biomD2/biomU2;
    if (taux2>1.0) taux2=1.0;
  }
 return pow(taux1,4.0); 
//  if (taux1<taux2) { return taux1; }
//	            else { return taux2; }  // on renvoie le plus petit des 2 taux calcul?s

} /* Fonction calcTauxSatis */
/****************************************************************************/
/****************************************************************************/
double calcDemandeVolume(pTAxe axe)
{
/* Calcule la demande en volume correspondant à la croissance en longueur
   pour un axe donn? */

  return pi*(axe->pointe->diametre)*(axe->pointe->diametre)*calcElongationPointe(axe->pointe,sol)/4.0;

} /* Fonction calcDemandeVolume */
/****************************************************************************/
/****************************************************************************/
void detruitAxe(pTAxe axeADetruire)
/* Supprime un axe en supprimant ses segments, puis l'axe lui-même */
{
pTSeg segCour, segAEnlever;

  /* Lib?rer tous les segments de cet axe */
  segCour=axeADetruire->premSeg;
  while (segCour->suiv!=NULL)
  {
    segAEnlever=segCour;
    segCour=segCour->suiv;
//    if (ndCour->suivSPere!=NULL) { printf("Probl?me : Axe ramifi? à enlever\n"); exit(1); }
    detruitSeg(segAEnlever);
  }
  detruitSeg(segCour); /* Enl?ve le segment apical */

  detruitPointe(axeADetruire->pointe);

  /* Enlever l'axe en m?moire */
  free(axeADetruire);

} /* Fonction detruitAxe */
/****************************************************************************/
/****************************************************************************/
int axeToutNecrose(pTAxe axe)
/* Cette fonction retourne la valeur 1 si l'axe a tous ses segments n?cros?s et 0 sinon */
{
pTSeg segCour;
int resu=1; // on initialise la valeur r?sultat à vrai (1)

  if (!axe->pointe->senile) resu=0; // non tout n?cros? si pointe non s?nile
  segCour=axe->premSeg;
  while (segCour!=NULL) {
    if (!segCour->necrose) resu=0;  // non tout n?cros? si un segment non n?cros?
    segCour=segCour->suiv;
  }
  return resu;

} /* Fonction axeToutNecrose */
/****************************************************************************/
/****************************************************************************/
void affecValNecroseAxe(pTAxe axe, int valNecrose)
/* Cette fonction affecte à chacun des segments de l'axe
   la valeur de necrose (0 ou 1) */
{
pTSeg segCour;

  segCour=axe->premSeg;
  while (segCour!=NULL)
  {
    segCour->necrose=valNecrose;
    segCour=segCour->suiv;
  }  // fin du while

} /* Fonction affecValNecroseAxe */
/****************************************************************************/
/****************************************************************************/
void affecValNecroseAmont(pTAxe axe, int valNecrose)
/* Cette fonction affecte a chacun des segments en amont de l'axe
la valeur de necrose (0 ou 1) */
{
pTSeg segCour;

  segCour=axe->dernSeg;
  while (segCour!=NULL)
  {
    segCour->necrose=valNecrose;
    segCour=segCour->prec;
  }

} /* Fonction affecValNecroseAmont */
/****************************************************************************/
/****************************************************************************/
void affecValDiamAxe(pTAxe axe, float diam)
/* Cette fonction affecte à chacun des segments de l'axe
   la valeur de diam?tre diam */
{
pTSeg segCour;

  segCour=axe->premSeg;
  while (segCour!=NULL)
  {
    segCour->diametre=diam;
    segCour=segCour->suiv;
  }  // fin du while

} /* Fonction affecValDiamAxe */
/****************************************************************************/
/****************************************************************************/
void increValDiamAmont(pTAxe axe, double diam, double coeff)
/* Cette fonction incremente le diametre de chacun des noeuds en amont de l'axe */
{
pTSeg segCour;
double section,diamInit;

  segCour=axe->premSeg->prec; // segment duquel l'axe est segment lat?ral
  while (segCour!=NULL)
  {
    diamInit=segCour->diametre;
    section=(pi*diamInit*diamInit/4.0)+(pi*coeff*diam*diam/4.0);
    segCour->diametre=sqrt(4.0*section/pi);
    segCour=segCour->prec;
  } // fin du while

} /* Fonction increValDiamAmont */
/****************************************************************************/
/****************************************************************************/
pTSysRac creeSR(void)
/* Cette fonction retourne une nouvelle variable de type pTSysRac,
  c'est-à-dire un pointeur sur le type SysRac */
{
pTSysRac sR;
  sR=(pTSysRac) malloc(sizeof(SysRac));
  if (sR==NULL)
  { printf("Memory issue allocation dans CreeSR \n"); exit(1); }

return sR;
} /* Fonction creeSR */
/****************************************************************************/
/****************************************************************************/
void enleveAxeSR(pTSysRac sR, pTAxe axeAEnlever)
/* Cette fonction enl?ve un axe dans la chaîne des axes du syst?me racinaire */
{
unsigned char axeDestructible=0;

  if (sR->premAxe==NULL)  /* Le syst?me racinaire est vide */
  {
    printf("ATTENTION, probleme dans enleveAxeSR, sR vide \n");
    exit(1);
  }
  else
  {
    if ((axeAEnlever->precedent!=NULL)&&(axeAEnlever->suivant!=NULL)) {
      // On pourra le supprimer, on refait le chaînage
      axeAEnlever->precedent->suivant=axeAEnlever->suivant;
      axeAEnlever->suivant->precedent=axeAEnlever->precedent;
      axeDestructible=1;
    } // fin du if !=NULL && !=NULL

    if ((axeAEnlever->precedent==NULL)&&(axeAEnlever->suivant!=NULL)) {
      // On pourra le supprimer, on refait le chaînage
      axeAEnlever->suivant->precedent=NULL;
      sR->premAxe=axeAEnlever->suivant;
      axeDestructible=1;
    } // fin du if ==NULL && !=NULL

    if ((axeAEnlever->precedent!=NULL)&&(axeAEnlever->suivant==NULL)) {
      // On pourra le supprimer, on refait le chaînage
      axeAEnlever->precedent->suivant=NULL;
      sR->dernAxe=axeAEnlever->precedent;
      axeDestructible=1;
    } // fin du if !=NULL && ==NULL

    if ((axeAEnlever->precedent==NULL)&&(axeAEnlever->suivant==NULL)) {
      // On ne pourra pas le supprimer, car il est seul
      axeDestructible=0;
    } // fin du if ==NULL && ==NULL

    if (axeDestructible) {
      sR->nbAxeSup++;
      detruitAxe(axeAEnlever); // D?truit ses segments, sa pointe, et lui-même
    }
  }
} /* Fonction enleveAxeSR */
/****************************************************************************/
/****************************************************************************/
pTSysRac initialiseSR(r3 origine)
{
/* Initialisation du syst?me racinaire */

pTSysRac sR;

  sR=creeSR();  /* Cr?ation d'un syst?me racinaire */

  sR->nbAxeForm=0;  /* Initialisation des variables */
  sR->nbAxeSup=0;
  sR->nbSegForm=0;
  sR->nbSem=0;  /* Nombre de racines s?minales ?mises */
  sR->nbAdv=0;  /* Nombre de racines adventives ?mises*/
  sR->longueur=0.0;
  sR->premAxe=NULL;
  sR->dernAxe=NULL;
  sR->tSatisMoy=1;
  sR->volTotPrec=0.0;
  sR->volTot=0.0;

  sR->origine[0]=origine[0];  /* Coordonn?es de l'origine du syst?me racinaire */
  sR->origine[1]=origine[1];
  sR->origine[2]=origine[2];

  sR->angDep=2.0*pi*dRandUnif();  /* Orientation */

  for (int i=0; i<NBPASMAX; i++) sR->tSatis[i]=1.0;

  return(sR);
}  /* Fonction initialiseSR */
/****************************************************************************/
/****************************************************************************/
float longSeg(pTSeg seg)
/* Calcule la longueur d'un segment */
{
  return sqrt(((seg->posE[0]-seg->posO[0])*(seg->posE[0]-seg->posO[0]))+
              ((seg->posE[1]-seg->posO[1])*(seg->posE[1]-seg->posO[1]))+
              ((seg->posE[2]-seg->posO[2])*(seg->posE[2]-seg->posO[2])));

}  /* Fonction longSeg */
/****************************************************************************/
/****************************************************************************/
int calcNouvNbSem(void)
{
/* Calcul du nouveau nombre de racines s?minales */

  int nouvNbSem;

  nouvNbSem=int (P_vitEmissionSem*temps);

  if (nouvNbSem>=P_nbMaxSem) nouvNbSem=P_nbMaxSem;

  return nouvNbSem;

}  /* Fonction calcNouvNbSem */
/****************************************************************************/
/****************************************************************************/
int calcNouvNbAdv(void)
{
/* Calcul du nouveau nombre de racines adventives */

  int nouvNbAdv;

  nouvNbAdv=int (P_vitEmissionAdv*(temps-P_ageEmissionAdv));

  if (nouvNbAdv>P_nbMaxAdv) nouvNbAdv=P_nbMaxAdv;

  return nouvNbAdv;

}  /* Fonction calcNouvNbAdv */
/****************************************************************************/
/****************************************************************************/
void emissionSemSR(pTSysRac sR)
{
/* Emission de nouveaux axes s?minaux sur le syst?me racinaire */

  pTAxe nouvAxe;
  int numSem, nbSemAEmettre;
  r3 vInit, dirInit;
  double angRot,angI;
//  float diamSem=0.0f;

  nbSemAEmettre=calcNouvNbSem() - sR->nbSem; /* Nombre de s?minales à ?mettre */
  for ((numSem=1); (numSem<=nbSemAEmettre); (numSem++)) /* Pour les nouvelles s?minales à ?mettre */
  {
//    printf("Je suis dans emissionSemSR %3i \n",sR->nbPrim);
    /* Calcul de la direction initiale de l'axe */
    if (sR->nbSem==0) angI=tireGaussien(0.0,0.06); // ?mission de la radicule proche verticale (gravitropisme initial fort)
      else angI=tireGaussien(P_angInitMoyVertSem,P_angInitETVertSem); // angle par rapport à la verticale
    vInit[0]=sin(angI);
    vInit[1]=0.0;
    vInit[2]=cos(angI);
    angRot=sR->angDep+tireAngRad();
    rotZ(vInit,dirInit,angRot);

    /* G?n?ration de l'axe et int?gration dans le syst?me racinaire */
//    diamSem=tireDiamSem();    
	  nouvAxe=initialiseAxe(sR->nbAxeForm+1,P_propDiamSem*P_diamMax,sR->origine,dirInit,NULL,NULL);
    ajouteAxeSR(sR,nouvAxe);
    sR->nbSem++;
  }

  }  /* Fonction emissionSemSR */
/****************************************************************************/
/****************************************************************************/
void emissionAdvSR(pTSysRac sR)
{
/* Emission de nouveaux axes adventifs sur le syst?me racinaire, 
   elles apparaissent sur la radicule (premi?re s?minale) */

  pTAxe nouvAxe;
  pTSeg segPere;  /* Segment sur lequel la racine adventive sera ?mise */
  int numAdv, nbAdvAEmettre;
  r3 vInit, dirInit, posInit;
  double angRot,angI,dBaseAdv,dBaseCour;

  nbAdvAEmettre=calcNouvNbAdv() - sR->nbAdv; /* Nombre de racines adventives à ?mettre */
  for ((numAdv=1); (numAdv<=nbAdvAEmettre); (numAdv++)) /* Pour les nouvelles adventives à ?mettre */
  {
//    printf("Je suis dans emissionAdvSR %3i \n",sR->nbTard);
	/* Calcul de la position initiale de l'axe */
	  /* Tirage de la distance à la base de cette adventive */
	  dBaseAdv=dRandUnif()*P_dBaseMaxAdv;

	  /* D?termination du segment p?re, sur le premier axe */
	  segPere=sR->premAxe->premSeg;
	  dBaseCour=longSeg(segPere);
	  while ((dBaseCour < dBaseAdv) && (segPere->suiv!=NULL)) {
		  segPere=segPere->suiv;
		  dBaseCour+=longSeg(segPere);
	  }

	  /* Position sur ce segment */
	  origineAdv(segPere,posInit);
    /* Calcul de la direction initiale de l'axe */
    angI=tireGaussien(P_angInitMoyVertAdv,P_angInitETVertAdv); // angle par rapport à la verticale
    vInit[0]=sin(angI);
    vInit[1]=0.0;
    vInit[2]=cos(angI);
    angRot=tireAngRad();
    rotZ(vInit,dirInit,angRot);

    /* G?n?ration de l'axe et int?gration dans le syst?me racinaire */
    nouvAxe=initialiseAxe(sR->nbAxeForm+1,P_propDiamAdv*P_diamMax,posInit,dirInit,sR->premAxe,segPere);
    ajouteAxeSR(sR,nouvAxe);
    sR->nbAdv++;
  }

  }  /* Fonction emissionAdvSR */
/****************************************************************************/
/****************************************************************************/
void calcVolProdSR(pTSysRac sR)
{
/* Calcul du volume racinaire produit sur la p?riode ?coul?e */
  int date;

  sR->volProd=0.0;
  for ((date=1); (date<=temps); (date++)) /* Boucle sur la p?riode ?coul?e */
  {
    //  sR->volProd+=sR->volDem[date]*sR->tSatis[date];     faux apr?s changement des calculs
  }

}  /* Fonction calcVolProdSR */



/****************************************************************************/
/****************************************************************************/
void readParamXML(void)

/* Fucntion to read the xml input file for the model */
{

  printf("Reading XML parameter file \n");
  tinyxml2::XMLDocument xmlDoc;
  xmlDoc.LoadFile(filename.c_str());

  // Get the root
  tinyxml2::XMLNode * pRoot = xmlDoc.FirstChildElement( "Plant" );

  P_duree = pRoot->FirstChildElement( "simtime" ) -> IntText();
  P_vitEmissionSem = pRoot->FirstChildElement( "erSem" )->FloatText();
  P_propDiamSem = pRoot->FirstChildElement( "dSem" )->FloatText();
  P_nbMaxSem = pRoot->FirstChildElement( "maxSem" )->IntText();
  P_ageEmissionAdv = pRoot->FirstChildElement( "ageAdv" )->FloatText();
  P_dBaseMaxAdv = pRoot->FirstChildElement( "distAdv" )->FloatText();
  P_vitEmissionAdv = pRoot->FirstChildElement( "erAdv" )->FloatText();
  P_propDiamAdv = pRoot->FirstChildElement( "dAdv" )->FloatText();
  P_nbMaxAdv = pRoot->FirstChildElement( "maxAdv" )->IntText();
  P_diamMin = pRoot->FirstChildElement( "dmin" )->FloatText();
  P_diamMax = pRoot->FirstChildElement( "dmax" )->FloatText();
  P_penteVitDiam = pRoot->FirstChildElement( "EL" )->FloatText();
  P_tendanceDirTropisme = pRoot->FirstChildElement( "TrT" )->IntText();
  P_intensiteTropisme = pRoot->FirstChildElement( "TrInt" )->FloatText();
  P_ageMaturitePointe = pRoot->FirstChildElement( "PDT" )->FloatText();
  P_distRamif = pRoot->FirstChildElement( "IPD" )->FloatText();
  P_probEmergeDmax = pRoot->FirstChildElement( "pdmax" )->FloatText();
  P_probEmergeDmin = pRoot->FirstChildElement( "pdmin" )->FloatText();
  P_propDiamRamif = pRoot->FirstChildElement( "RDM" )->FloatText();
  P_coeffVarDiamRamif = pRoot->FirstChildElement( "CVDD" )->FloatText();
  P_TMD = pRoot->FirstChildElement( "TMD" )->FloatText();
  P_penteDureeCroissDiam2 = pRoot->FirstChildElement( "GDs" )->FloatText();
  P_penteDureeVieDiamTMD = pRoot->FirstChildElement( "LDC" )->FloatText();
  P_coeffCroissRad = pRoot->FirstChildElement( "SGC" )->FloatText();

  P_exportType = pRoot->FirstChildElement( "exportType" )->IntText();
  sprintf(P_exportName, "%s", pRoot->FirstChildElement( "exportName" )->GetText());


}
/****************************************************************************/
/****************************************************************************/
void litParam(void)

/* Fonction de lecture des parametres de la simulation */
{
char bid[MAXLINE];

  // Dur?e de simulation
  fscanf(FPar,"%i",&P_duree);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_vitEmissionSem);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_propDiamSem);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%i",&P_nbMaxSem);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_ageEmissionAdv);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_dBaseMaxAdv);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_vitEmissionAdv);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_propDiamAdv);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%i",&P_nbMaxAdv);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_diamMin);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_diamMax);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_penteVitDiam);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%i",&P_tendanceDirTropisme);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_intensiteTropisme);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_ageMaturitePointe);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_distRamif);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_probEmergeDmax);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_probEmergeDmin);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_propDiamRamif);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_coeffVarDiamRamif);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_TMD);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_penteDureeCroissDiam2);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_penteDureeVieDiamTMD);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

  fscanf(FPar,"%f",&P_coeffCroissRad);
  fgets(bid,MAXLINE-1,FPar); // reste de la ligne

} /* Fonction litParam */
/****************************************************************************/
/****************************************************************************/
void origineEmission(pTAxe nouvAxe)
{
nouvAxe->pointe->coord[0]=sR->origine[0];
nouvAxe->pointe->coord[1]=sR->origine[1];
nouvAxe->pointe->coord[2]=sR->origine[2];
} /* Fonction origineEmission */
/****************************************************************************/
/****************************************************************************/
void orienteEmission(pTAxe nouvAxe, int num)
{
double angRot,angI;
r3 vInit;

angI=tireGaussien(P_angInitMoyVertSem,P_angInitETVertSem);
vInit[0]=sin(angI);
vInit[1]=0.0;
vInit[2]=cos(angI);

angRot=sR->angDep+(2*pi*num/P_nbMaxSem);
rotZ(vInit,nouvAxe->pointe->dirCroiss,angRot);
} /* Fonction orienteEmission */
/****************************************************************************/
/****************************************************************************/
float volPrimSeg(pTSeg seg)
/* Calcule le volume primaire du segment */
{
  return 0.25*pi*seg->axe->pointe->diametre*seg->axe->pointe->diametre*
  sqrt(((seg->posE[0]-seg->posO[0])*(seg->posE[0]-seg->posO[0]))+
       ((seg->posE[1]-seg->posO[1])*(seg->posE[1]-seg->posO[1]))+
       ((seg->posE[2]-seg->posO[2])*(seg->posE[2]-seg->posO[2])));

}  /* Fonction volPrimSeg */
/****************************************************************************/
/****************************************************************************/
float volTotalSeg(pTSeg seg)
/* Calcule le volume total du segment */
{
  return 0.25*pi*seg->diametre*seg->diametre*sqrt(((seg->posE[0]-seg->posO[0])*(seg->posE[0]-seg->posO[0]))+
                                                  ((seg->posE[1]-seg->posO[1])*(seg->posE[1]-seg->posO[1]))+
                                                  ((seg->posE[2]-seg->posO[2])*(seg->posE[2]-seg->posO[2])));

}  /* Fonction volTotalSeg */
/****************************************************************************/
/****************************************************************************/
float distHorSeg(pTSeg seg)
/* Calcule la distance horizontale d'un segment */
{
  return sqrt(((seg->posE[0]+seg->posO[0])*(seg->posE[0]+seg->posO[0])/4)+
              ((seg->posE[1]+seg->posO[1])*(seg->posE[1]+seg->posO[1])/4));

}  /* Fonction distHorSeg */
/****************************************************************************/
/****************************************************************************/
void calcLimitesSR(pTSysRac sR)
{
/* Calcul des limites du syst?me racinaire et de quelques autres variables */
  pTAxe axeCour;
  pTSeg segCour;
  float distHor,distHorLong,profLong,longS,profS,amplMax;

  // Initialisation des variables
  sR->volPrim=0.0;  // volume des structures primaires
  sR->secPointe=0.0;  // section totale des pointes actives
  sR->volTot=0.0;   // volume total
  sR->longueur=0.0;  // longueur
  sR->diamMax=-1.0e10; // diam?tre maximal, du plus gros segment
  sR->distMax=-1.0e10;  // extension maximale
  sR->profMax=-1.0e10;  // profondeur maximale

  sR->xbinf=+1.0e10; sR->ybinf=+1.0e10; sR->zbinf=+1.0e10; // initialisation des valeurs
  sR->xbsup=-1.0e10; sR->ybsup=-1.0e10; sR->zbsup=-1.0e10;

  distHorLong=0.0;
  profLong=0.0;

  axeCour=sR->premAxe;
  while (axeCour!=NULL)  // Calcul du volume "demand?"
  {
    sR->volTot+=pi*axeCour->pointe->diametre*axeCour->pointe->diametre*axeCour->pointe->longueur/4;
    sR->volPrim+=pi*axeCour->pointe->diametre*axeCour->pointe->diametre*axeCour->pointe->longueur/4;
    if ((axeCour->pointe->mature)&&(!axeCour->pointe->senile))
      sR->secPointe+=pi*axeCour->pointe->diametre*axeCour->pointe->diametre/4;
    segCour=axeCour->premSeg;
    while (segCour!=NULL) { // Tant que ce segment existe
      // Calculs sur le segment courant segCour
      if (segCour->posO[0] < sR->xbinf) { sR->xbinf=segCour->posO[0]; }
      if (segCour->posE[0] < sR->xbinf) { sR->xbinf=segCour->posE[0]; }
      if (segCour->posO[0] > sR->xbsup) { sR->xbsup=segCour->posO[0]; }
      if (segCour->posE[0] > sR->xbsup) { sR->xbsup=segCour->posE[0]; }

      if (segCour->posO[1] < sR->ybinf) { sR->ybinf=segCour->posO[1]; }
      if (segCour->posE[1] < sR->ybinf) { sR->ybinf=segCour->posE[1]; }
      if (segCour->posO[1] > sR->ybsup) { sR->ybsup=segCour->posO[1]; }
      if (segCour->posE[1] > sR->ybsup) { sR->ybsup=segCour->posE[1]; }

      if (segCour->posO[2] < sR->zbinf) { sR->zbinf=segCour->posO[2]; }
      if (segCour->posE[2] < sR->zbinf) { sR->zbinf=segCour->posE[2]; }
      if (segCour->posO[2] > sR->zbsup) { sR->zbsup=segCour->posO[2]; }
      if (segCour->posE[2] > sR->zbsup) { sR->zbsup=segCour->posE[2]; }

      if (segCour->diametre > sR->diamMax) { sR->diamMax=segCour->diametre; }

      distHor=distHorSeg(segCour);
      if (distHor > sR->distMax) { sR->distMax=distHor; }

      if (segCour->posO[2]>segCour->posE[2]) profS=segCour->posO[2]; else profS=segCour->posE[2];
      if (profS > sR->profMax) { sR->profMax=profS; }

      sR->volTot+=volTotalSeg(segCour);
      sR->volPrim+=volPrimSeg(segCour);
      longS=longSeg(segCour);
      sR->longueur+=longS;
      distHorLong+=distHor*longS;
      profLong+=profS*longS;

      segCour=segCour->suiv;
    }  // fin du while segCour
    axeCour=axeCour->suivant;
  }  // fin du while axeCour

  sR->xbinf=sR->xbinf-d2; sR->xbsup=sR->xbsup+d2;
  sR->ybinf=sR->ybinf-d2; sR->ybsup=sR->ybsup+d2;
  sR->zbinf=sR->zbinf-d2; sR->zbsup=sR->zbsup+d2;

  // Calcul de la maille de sol, en fonction de l'amplitude à balayer

  amplMax=0.0;
  if ((sR->xbsup-sR->xbinf) > amplMax) amplMax=sR->xbsup-sR->xbinf;
  if ((sR->ybsup-sR->ybinf) > amplMax) amplMax=sR->ybsup-sR->ybinf;
  if ((sR->zbsup-sR->zbinf) > amplMax) amplMax=sR->zbsup-sR->zbinf;

  maille=amplMax/(NBCASEMAX-1);
  if (maille<mailleMin) maille=mailleMin;
  volElemSol=maille*maille*maille;

//  maille=5.0; volElemSol=maille*maille*maille;

  sR->distMoy=distHorLong/sR->longueur;
  sR->profMoy=profLong/sR->longueur;

//  printf(" xbinf :%7.2f",sR->xbinf); printf(" xbsup :%7.2f\n",sR->xbsup);
//  printf(" ybinf :%7.2f",sR->ybinf); printf(" ybsup :%7.2f\n",sR->ybsup);
//  printf(" zbinf :%7.2f",sR->zbinf); printf(" zbsup :%7.2f\n",sR->zbsup);
//  printf(" amplMax :%7.2f",amplMax); printf(" maille :%7.2f\n",maille);


}  /* Fonction calcLimitesSR */
/****************************************************************************/
/****************************************************************************/
void calcVolumesBiomassesSR(pTSysRac sR)
{
/* Calcul des diff?rents volumes pour le syst?me racinaire avant le d?velopement sur le pas de temps */
  pTAxe axeCour;
  pTSeg segCour;

  // Initialisation des variables
  sR->volTotPrec=sR->volTot;  // volume qui avait ?t? calcul? au pas de temps pr?c?dent
  
  sR->volPrim=0.0;  // volume des structures primaires
  sR->volTot=0.0;   // volume total ? calculer sur l'ensemble des segments et des pointes de tous les axes

  axeCour=sR->premAxe;
  while (axeCour!=NULL)  // Sur l'ensemble des axes
  {
    sR->volTot+=pi*axeCour->pointe->diametre*axeCour->pointe->diametre*axeCour->pointe->longueur/4;  // volume de la pointe
    sR->volPrim+=pi*axeCour->pointe->diametre*axeCour->pointe->diametre*axeCour->pointe->longueur/4;

    segCour=axeCour->premSeg;
    while (segCour!=NULL) { // Tant qu'il y a des segments sur l'axe
      // Cumul des volumes
      sR->volTot+=volTotalSeg(segCour);
      sR->volPrim+=volPrimSeg(segCour);
      segCour=segCour->suiv;
    }  // fin du while segCour
    axeCour=axeCour->suivant;
  }  // fin du while axeCour

  sR->biomUtilisee[temps-1]=P_TMD*(sR->volTot - sR->volTotPrec)/1000;  // Volume calcul? comme la diff?rence entre volume courant et le volume au jour pr?c?dent (NB peut ?tre n?gatif)
  if (sR->biomUtilisee[temps-1] < 0.0001)  sR->biomUtilisee[temps-1]=0.0001;  // Correction ?ventuelle, pour ?tre s?r de ne pas avoir une masse n?gative ou trop faible
  
}  /* Fonction calcVolumesBiomassesSR */
/****************************************************************************/
/*************************************************************************/
void translateSR(pTSysRac sR)
/* Translate le syst?me racinaire de façon à ce que tout se passe en territoire
 positif et d?marre de 0*/

{
pTAxe axeCour;
pTSeg segCour;

  axeCour=sR->premAxe;
  while (axeCour!=NULL)  // Calcul du volume "demand?"
  {
    // Translation de la pointe de l'axe
    axeCour->pointe->coord[0] -= sR->xbinf;
    axeCour->pointe->coord[1] -= sR->ybinf;
    axeCour->pointe->coord[2] -= sR->zbinf;

    segCour=axeCour->premSeg;
    while (segCour!=NULL) { // Tant qu'il y a des segments sur l'axe
      // Translation du segment segCour
      segCour->posO[0] -= sR->xbinf;
      segCour->posO[1] -= sR->ybinf;
      segCour->posO[2] -= sR->zbinf;

      segCour->posE[0] -= sR->xbinf;
      segCour->posE[1] -= sR->ybinf;
      segCour->posE[2] -= sR->zbinf;

      segCour=segCour->suiv;
    }
    axeCour=axeCour->suivant;
  }

  sR->xbsup-=sR->xbinf; sR->ybsup-=sR->ybinf; sR->zbsup-=sR->zbinf;
  sR->xbinf=0.0; sR->ybinf=0.0; sR->zbinf=0.0;

} /* Fonction translateSR */
/*************************************************************************/
/*************************************************************************/
void initialiseTabSol(void)
/* Initialise le tableau sol. La valeur 6 signifie que le point est à
   distance sup?rieure à d1 et d2 et d3 et d4 et d5 */
{
  for (int i=0; i<=NBCASEMAX; i++) { // pour le sol, on initialise à 6
    for (int j=0; j<=NBCASEMAX; j++) {
      for (int k=0; k<=NBCASEMAX; k++) { vox[i][j][k]=6; }
    }
  }

} /* Fonction initialiseTabSol */
/*************************************************************************/
/*************************************************************************/
int rangCase(float coord)
{ // renvoie le rang de la case du tableau des voxels dans laquelle est cette coordonn?e

  int rang=int (coord/maille);
  if (rang<0) return 0; else if (rang>NBCASEMAX) return NBCASEMAX;
  return rang;

} /* Fonction rangCase */
/*************************************************************************/
/*************************************************************************/
float coordPointCase(int rangCase)
{ // renvoie les coordonn?es d'un point dans la case

  return (((rangCase+0.5)*maille)+(0.5*maille*dRandUnif()));

} /* Fonction coordPointCase */
/*************************************************************************/
/*************************************************************************/
float coordCentreCase(int rangCase)
{ // renvoie les coordonn?es du centre de la case

  return ((rangCase+0.5)*maille);

} /* Fonction coordCentreCase */
/*************************************************************************/
/*************************************************************************/
void calcDistancesSR(pTSysRac sR)
/* Calcule les distances entre mailles du sol et segments racinaires
   et calcule les volumes colonis?s */
{
  pTAxe axeCour;
  pTSeg segCour;

  float xp1,yp1,zp1, xp2,yp2,zp2, dx,dy,dz, dM, distCour,
        xMin, xMax, yMin, yMax, zMin, zMax, xS, yS, zS, xProj, yProj, zProj;
//  float dist1,dist2;   // si besoin

  axeCour=sR->premAxe;
  while (axeCour!=NULL) // Tant qu'il y a des axes dans le syst?me racinaire
  {
    segCour=axeCour->premSeg;
    if (segCour->complet) {
      while (segCour!=NULL)  // Tant qu'il y a des segments sur l'axe
      {
        xp1=segCour->posO[0]; yp1=segCour->posO[1]; zp1=segCour->posO[2];
        xp2=segCour->posE[0]; yp2=segCour->posE[1]; zp2=segCour->posE[2];

        // Calcul des limites du domaine à explorer pour ce segment
        if (xp1<xp2) { xMin=xp1 - d5; xMax=xp2 + d5; }
        else { xMin=xp2 - d5; xMax=xp1 + d5; }
        if (yp1<yp2) { yMin=yp1 - d5; yMax=yp2 + d5; }
        else { yMin=yp2 - d5; yMax=yp1 + d5; }
        if (zp1<zp2) { zMin=zp1 - d5; zMax=zp2 + d5; }
        else { zMin=zp2 - d5; zMax=zp1 + d5; }

          // balayage de ce domaine pertinent et calcul des distances
        for (int caseX=rangCase(xMin); caseX<=rangCase(xMax); caseX++) {
          for (int caseY=rangCase(yMin); caseY<=rangCase(yMax); caseY++) {
            for (int caseZ=rangCase(zMin); caseZ<=rangCase(zMax); caseZ++) {

              xS=coordCentreCase(caseX); yS=coordCentreCase(caseY); zS=coordCentreCase(caseZ);

              // On calcule le projet? du point sol sur la droite contenant p1 et p2
              dx=xp2-xp1; dy=yp2-yp1; dz=zp2-zp1;
              dM=(dx*xS)+(dy*yS)+(dz*zS);

              xProj=((dM*dx)+(xp1*((dy*dy)+(dz*dz)))-(zp1*dx*dz)-(yp1*dx*dy))/((dx*dx)+(dy*dy)+(dz*dz));

              if ((xProj<=xp1 && xProj>=xp2)|(xProj>=xp1 && xProj<=xp2))
              { // Le projet? est entre les deux points du segment
                yProj=((dM*dy)+(yp1*((dx*dx)+(dz*dz)))-(xp1*dx*dy)-(zp1*dz*dy))/((dx*dx)+(dy*dy)+(dz*dz));
                zProj=((dM*dz)+(zp1*((dx*dx)+(dy*dy)))-(yp1*dz*dy)-(xp1*dz*dx))/((dx*dx)+(dy*dy)+(dz*dz));
                distCour=sqrt(((xProj-xS)*(xProj-xS))+((yProj-yS)*(yProj-yS))+((zProj-zS)*(zProj-zS)));
              }  // fin du if
              else
              {  // Le projet? est à l'ext?rieur du segment
              /*
                dist1=sqrt(((xp1-xS)*(xp1-xS))+((yp1-yS)*(yp1-yS))+
                           ((zp1-zS)*(zp1-zS)));
                dist2=sqrt(((xp2-xS)*(xp2-xS))+((yp2-yS)*(yp2-yS))+
                           ((zp2-zS)*(zp2-zS)));
                if (dist1<dist2) { distCour=dist1; } else { distCour=dist2; }
              */
              distCour=2000.0;
              }  // fin du else

              if ((distCour<=d5)&&(vox[caseX][caseY][caseZ]>5)) { vox[caseX][caseY][caseZ]=5; }
              if ((distCour<=d4)&&(vox[caseX][caseY][caseZ]>4)) { vox[caseX][caseY][caseZ]=4; }
              if ((distCour<=d3)&&(vox[caseX][caseY][caseZ]>3)) { vox[caseX][caseY][caseZ]=3; }
              if ((distCour<=d2)&&(vox[caseX][caseY][caseZ]>2)) { vox[caseX][caseY][caseZ]=2; }
              if ((distCour<=d1)&&(vox[caseX][caseY][caseZ]>1)) { vox[caseX][caseY][caseZ]=1; }

            } // for du for caseZ
          } // for du for caseY
        } // for du for caseX
        segCour=segCour->suiv;
      }  // fin du while (segCour!=NULL)
    } // fin du if (segCour->complet)
    axeCour=axeCour->suivant;
  } // fin du while (axeCour!=NULL)

  sR->volSolD1=0.0;  // initialisation avant cumul
  sR->volSolD2=0.0;  // initialisation avant cumul
  sR->volSolD3=0.0;  // initialisation avant cumul
  sR->volSolD4=0.0;  // initialisation avant cumul
  sR->volSolD5=0.0;  // initialisation avant cumul
  
  for (int i=0; i<=NBCASEMAX; i++) {
    for (int j=0; j<=NBCASEMAX; j++) {
      for (int k=0; k<=NBCASEMAX; k++) {
        if (vox[i][j][k]==5) {
          sR->volSolD5+=volElemSol;
        } // fin du if
        if (vox[i][j][k]==4) {
          sR->volSolD5+=volElemSol; sR->volSolD4+=volElemSol;
        } // fin du if
        if (vox[i][j][k]==3) {
          sR->volSolD5+=volElemSol; sR->volSolD4+=volElemSol; sR->volSolD3+=volElemSol;
        } // fin du if
        if (vox[i][j][k]==2) {
          sR->volSolD5+=volElemSol; sR->volSolD4+=volElemSol; sR->volSolD3+=volElemSol; sR->volSolD2+=volElemSol;
        } // fin du if
        if (vox[i][j][k]==1) {
          sR->volSolD5+=volElemSol; sR->volSolD4+=volElemSol; sR->volSolD3+=volElemSol; sR->volSolD2+=volElemSol; sR->volSolD1+=volElemSol;
        } // fin du if
      }  // fin du for sur k
    }  // fin du for sur j
  }  // fin du for sur i


}  /* Fonction calcDistancesSR  */
/****************************************************************************/
/****************************************************************************/
void developpeSR(pTSysRac sR)
{
/* D?veloppement : croissance et ramification de chaque axe du syst?me */
pTAxe axeCour;
double tS;

  calcVolumesBiomassesSR(sR);  // On calcule les volumes et les biomasses du syst?me racinaire pour la demande et le taux de satisfaction
  
  sR->tSatis[temps]=calcTauxSatis(sR);  // On calcule le taux de satifaction en utilisant l'info sur la biomasse utilis?e au jour pr?c?dent

  axeCour=sR->premAxe;
  while (axeCour!=NULL)  // D?veloppement
  {
    developpeAxe(axeCour,sR->tSatis[temps]); // d?veloppe l'axe (croissance, ramif)
    developpePointe(axeCour->pointe);    // modifie les attributs de la pointe
    axeCour=axeCour->suivant;
  }

  // printf(" Volume demand? : %16.5f \n",volumeDem);
  // printf(" NbRac : %6i \n",sR->nbAxeForm);

}  /* Fonction developpeSR */
/****************************************************************************/
/****************************************************************************/
void mortaliteSR(pTSysRac sR)
{
pTAxe axeCour, axeAEnlever;

  // Premier passage : calcul de la s?nilit? et affectation n?crose sur l'ensemble des axes */
  axeCour=sR->premAxe; // Dans le sens de premiers vers les derniers
  while (axeCour!=NULL)
  {
    if (axeCour->pointe->senile)
    { /* L'axe est n?cros? */
      affecValNecroseAxe(axeCour, 1);
    }
    else
    {  /* L'axe n'est pas n?cros? */
      /* Tous les noeuds en amont de la pointe ne sont pas necros?s non plus */
      affecValNecroseAmont(axeCour, 0);
    }

    axeCour=axeCour->suivant;
  }

  // Calcul de l'?lagage, enl?vement des axes tout n?cros?s
  axeCour=sR->dernAxe; // Dans le sens de derniers vers les premiers
  while (axeCour!=NULL)
  {
    if (axeToutNecrose(axeCour))
    {
      axeAEnlever=axeCour;
      axeCour=axeCour->precedent;
      if (axeAEnlever->pere!=NULL) enleveAxeSR(sR,axeAEnlever);
    }
    else axeCour=axeCour->precedent;
  }

}  /* Fonction mortaliteSR */
/****************************************************************************/
/****************************************************************************/
void croissanceRadialeSR(pTSysRac sR, float coeffCroiss)
{
pTAxe axeCour;
float diam;

  /* Premier passage, initialisation aux diametres primaires */
  axeCour=sR->dernAxe;
  while (axeCour!=NULL)
  {
    diam=axeCour->pointe->diametre;
    affecValDiamAxe(axeCour, diam);
    axeCour=axeCour->precedent;
  }

  /* Deuxi?me passage, avec incr?ment des diametres */
  axeCour=sR->dernAxe;
  while (axeCour!=NULL)
  {
    /* les noeuds en amont sont incr?ment?s si axe en croissance (pointe mature et non senile) */
    if ((axeCour->pointe->mature)&&(!axeCour->pointe->senile))
    {
      diam=axeCour->pointe->diametre;
      increValDiamAmont(axeCour, diam, coeffCroiss);
    }
    axeCour=axeCour->precedent;
  }

}  /* Fonction croissanceRadialeSR */
/****************************************************************************/
/****************************************************************************/
void imprimeSeg(pTSeg seg)
/* Imprime un segment sur le fichier des segments */
{
/*
  long int suivant,precedent;

  if (seg->prec==NULL) precedent=0;
  else precedent=seg->prec->num;

  if (seg->suiv==NULL) suivant=0;
  else suivant=seg->suiv->num;

  fprintf(FSeg,"%5li %5i %5li %5li %5li %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
    seg->num,seg->jourForm,seg->axe->num,suivant,precedent,seg->diametre,
    seg->posO[0],seg->posO[1],seg->posO[2],seg->posE[0],seg->posE[1],seg->posE[2]);
*/
  fprintf(FSeg,"%5li %5i %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
    seg->axe->num,seg->jourForm,seg->diametre,
    seg->posO[0],seg->posO[1],seg->posO[2],seg->posE[0],seg->posE[1],seg->posE[2]);

}  /* Fonction imprimeSeg */
/****************************************************************************/
/****************************************************************************/
void imprimeSRGlobal(pTSysRac sR)
{
  /* Imprime un ensemble de variables pour qualifier globalement le syst?me racinaire */
//  fprintf(FSynth,"        longueur           volTot           volPrim          volProd          tSatisMoy          profMax          profMoy          distMax          distMoy         volSolD1         volSolD2\n");

	fprintf(FSynth,"%16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f \n",
          sR->longueur,sR->volTot,sR->tSatisMoy,sR->profMax,sR->distMax,sR->volSolD1/1000.0,sR->volSolD2/1000.0,sR->volSolD3/1000.0,sR->volSolD4/1000.0,sR->volSolD5/1000.0);

//  fprintf(FSynth,"        longueur           volTot           volPrim          volProd          secPointe          diamMax        tSatisMoy          profMax          profMoy          distMax          distMoy         volSolD1         volSolD2\n");
/*
	fprintf(FSynth,"%16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f %16.2f\n",
          sR->longueur,sR->volTot,sR->volPrim,sR->volProd,sR->tSatisMoy,sR->profMax,
          sR->profMoy,sR->distMax,sR->distMoy,sR->volSolD1/1000.0,sR->volSolD2/1000.0);
*/
}  /* Fonction imprimeSRGlobal */
/****************************************************************************/
/*************************************************************************/
void imprimeSolColonise(int distance)
/* Imprime les cellules de sol colonis? */
{
/*
  // Impression de l'entête
  fprintf(FVox,"caseX caseY caseZ\n");

  for (int i=0; i<NBCASEMAX; i++) {
    for (int j=0; j<NBCASEMAX; j++) {
      for (int k=0; k<NBCASEMAX; k++) {
        if (vox[i][j][k]<=distance) fprintf(FVox,"%5i %5i %5i\n",i,j,k);
      } // fin du for k
    }   // fin du for j
  }     // fin du for i
*/

} /* Fonction imprimeSolColonise */
/*************************************************************************/
/****************************************************************************/
void imprimeAudit(void)
{/*

  fprintf(FAudit,"voldispo                voldem             tsatis\n");
  for (int pas=1; pas<NBPASMAX; pas++)
  {
    fprintf(FAudit,"%16.2f %16.2f %8.5f\n",sR->volMax[pas],sR->volDem[pas],sR->tSatis[pas]);
  }
*/
}  /* Fonction imprimeAudit */
/****************************************************************************/
/****************************************************************************/



void imprimeSRSegmentsEnteteRSML(void)
{   /* Imprime l'entête du fichier contenant les noeuds du système racinaire */
  FSeg=fopen(P_exportName2,"w"); //MODIFBEN passage à outputname2, incluant le pas de temps

 // fprintf(FSeg,"NumSeg Jour NumAxe Suiv Prec Diam     X1       Y1       Z1      X2       Y2       Z2\n");
  //fprintf(FSeg,"NumAxe Jour Diam     X1       Y1       Z1      X2       Y2       Z2\n");
fprintf(FSeg,"<?xml version='1.0' encoding='UTF-8'?>\n");
fprintf(FSeg,"<rsml xmlns:po='http://www.plantontology.org/xml-dtd/po.dtd'>\n");
fprintf(FSeg,"  <metadata>\n");
fprintf(FSeg,"    <version>1</version>\n");
fprintf(FSeg,"    <last-modified>today</last-modified>\n");
fprintf(FSeg,"    <software>archisimple</software>\n");
fprintf(FSeg,"    <species>Arabidopsis thaliana</species>\n");
fprintf(FSeg,"    <medium>null</medium>\n");
fprintf(FSeg,"    <treatment>null</treatment>\n");
fprintf(FSeg,"    <age_of_plant>");
fprintf(FSeg,"%i",temps);
fprintf(FSeg,"</age_of_plant>\n");

// Print the parameters of the simulation
fprintf(FSeg,"    <parameters>\n");  

fprintf(FSeg,"      <parameter name='P_duree'>");
fprintf(FSeg,"%i",P_duree);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_vitEmissionSem'>");
fprintf(FSeg,"%f",P_vitEmissionSem);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_nbMaxSem'>");
fprintf(FSeg,"%i",P_nbMaxSem);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_propDiamSem'>");
fprintf(FSeg,"%f",P_propDiamSem);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_ageEmissionAdv'>");
fprintf(FSeg,"%f",P_ageEmissionAdv);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_dBaseMaxAdv'>");
fprintf(FSeg,"%f",P_dBaseMaxAdv);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_vitEmissionAdv'>");
fprintf(FSeg,"%f",P_vitEmissionAdv);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_propDiamAdv'>");
fprintf(FSeg,"%f",P_propDiamAdv);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_nbMaxAdv'>");
fprintf(FSeg,"%i",P_nbMaxAdv);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_diamMin'>");
fprintf(FSeg,"%f",P_diamMin);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_diamMax'>");
fprintf(FSeg,"%f",P_diamMax);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_penteVitDiam'>");
fprintf(FSeg,"%f",P_penteVitDiam);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_tendanceDirTropisme'>");
fprintf(FSeg,"%i",P_tendanceDirTropisme);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_intensiteTropisme'>");
fprintf(FSeg,"%f",P_intensiteTropisme);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_ageMaturitePointe'>");
fprintf(FSeg,"%f",P_ageMaturitePointe);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_distRamif'>");
fprintf(FSeg,"%f",P_distRamif);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_probEmergeDmax'>");
fprintf(FSeg,"%f",P_probEmergeDmax);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_probEmergeDmin'>");
fprintf(FSeg,"%f",P_probEmergeDmin);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_propDiamRamif'>");
fprintf(FSeg,"%f",P_propDiamRamif);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_coeffVarDiamRamif'>");
fprintf(FSeg,"%f",P_coeffVarDiamRamif);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_TMD'>");
fprintf(FSeg,"%f",P_TMD);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_penteDureeCroissDiam2'>");
fprintf(FSeg,"%f",P_penteDureeCroissDiam2);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_penteDureeVieDiamTMD'>");
fprintf(FSeg,"%f",P_penteDureeVieDiamTMD);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"      <parameter name='P_coeffCroissRad'>");
fprintf(FSeg,"%f",P_coeffCroissRad);
fprintf(FSeg,"</parameter>\n");

fprintf(FSeg,"    </parameters>\n");  
fprintf(FSeg,"    <user>globet</user>\n");
fprintf(FSeg,"    <file-key>myimage</file-key>\n");
fprintf(FSeg,"    <property-definitions>\n");
fprintf(FSeg,"      <property-definition>\n");
fprintf(FSeg,"          <label>diameter</label>\n");
fprintf(FSeg,"            <type>float</type>\n");
fprintf(FSeg,"            <unit>cm</unit>\n");
fprintf(FSeg,"      </property-definition>\n");
fprintf(FSeg,"      <property-definition>\n");
fprintf(FSeg,"          <label>age</label>\n");
fprintf(FSeg,"            <type>int</type>\n");
fprintf(FSeg,"            <unit>day</unit>\n");
fprintf(FSeg,"      </property-definition>\n");
fprintf(FSeg,"    </property-definitions>\n");
fprintf(FSeg,"  </metadata>\n");
fprintf(FSeg,"  <scene>\n");
fprintf(FSeg,"    <plant>  \n");

}  /* Fonction imprimeSRSegmentsEntete */

 /****************************************************************************/
/****************************************************************************/

void imprimeSegRSML(pTSeg seg, bool last)
/* Imprime un segment sur le fichier des segments */
{

     if(!last){
          fprintf(FSeg,"             <point x='");
          fprintf(FSeg,"%f",seg->posO[0]);
          fprintf(FSeg,"' y='");
          fprintf(FSeg,"%f",seg->posO[2]);
          fprintf(FSeg,"' z='");
          fprintf(FSeg,"%f",seg->posO[1]);
          fprintf(FSeg,"'/>\n");
      }
       if(last){
          fprintf(FSeg,"             <point x='");
          fprintf(FSeg,"%f",seg->posE[0]);
          fprintf(FSeg,"' y='");
          fprintf(FSeg,"%f",seg->posE[2]);
          fprintf(FSeg,"' z='");
          fprintf(FSeg,"%f",seg->posE[1]);
          fprintf(FSeg,"'/>\n");
      }


}  /* Fonction imprimeSegRSML */


/****************************************************************************/
/****************************************************************************/
void printAgeRSML(pTSeg seg)
/* Imprime un segment sur le fichier des segments */
{

    fprintf(FSeg,"             <sample>");
    fprintf(FSeg,"%i", seg->jourForm);
    fprintf(FSeg,"</sample>\n");

}  /* Fonction imprimeSeg */

 /****************************************************************************/
/****************************************************************************/
void printDiamRSML(pTSeg seg)
/* Imprime un segment sur le fichier des segments */
{

    fprintf(FSeg,"             <sample>");
    fprintf(FSeg,"%f", seg->diametre);
    fprintf(FSeg,"</sample>\n");

}  /* Fonction imprimeSeg */

/****************************************************************************/
/****************************************************************************/
void imprimeAxeSegmentsRSML(pTAxe axe)
{   /* Imprime les segments de l'axe */

  pTSeg segCour;
      fprintf(FSeg,"        <geometry>  \n");
      fprintf(FSeg,"          <polyline>  \n");
      segCour=axe->premSeg;
      while (segCour!=NULL) {
        if (segCour->complet) imprimeSegRSML(segCour, false);
        segCour=segCour->suiv;
      }
      segCour=axe->dernSeg;
      if (segCour->complet) imprimeSegRSML(segCour, true);

      fprintf(FSeg,"          </polyline>  \n");
      fprintf(FSeg,"        </geometry>  \n");

      // Print the diameters
      fprintf(FSeg,"        <functions>  \n");
      fprintf(FSeg,"          <function name='diameter' domain='polyline'>  \n");
      segCour=axe->premSeg;
      while (segCour!=NULL) {
        if (segCour->complet) printDiamRSML(segCour);
        segCour=segCour->suiv;
      }
      segCour=axe->dernSeg;
      if (segCour->complet) printDiamRSML(segCour);
      fprintf(FSeg,"          </function>  \n");

      // Print the age
      fprintf(FSeg,"          <function name='age' domain='polyline'>  \n");
      segCour=axe->premSeg;
      while (segCour!=NULL) {
        if (segCour->complet) printAgeRSML(segCour);
        segCour=segCour->suiv;
      }
      segCour=axe->dernSeg;
      if (segCour->complet) printAgeRSML(segCour);
      fprintf(FSeg,"          </function>  \n");
      fprintf(FSeg,"        </functions>  \n");
}  /* Fonction imprimeAxeSegments */






/****************************************************************************/
/****************************************************************************/
void imprimeSRSegmentsRSML(pTSysRac sR)
{  /* Imprime l'ensemble des segments du système racinaire */

  pTAxe axeCour;
  pTAxe axeCour1;
  pTAxe axeCour2;
  pTAxe axeCour3;
  pTAxe axeCour4;

  imprimeSRSegmentsEnteteRSML();

  axeCour=sR->premAxe;
  while (axeCour!=NULL)
  {

    // Print the primary axes
    if(axeCour->pere==NULL){

      // Count the number of segments
      pTSeg segCour;
      int count = 0;
      segCour=axeCour->premSeg;
      while (segCour!=NULL) {
        if (segCour->complet) count = count+1;
        segCour=segCour->suiv;
      }
      if(count > 1){
        fprintf(FSeg,"      <root ID='");
        fprintf(FSeg,"%li",axeCour->num);
        fprintf(FSeg,"' label='root' po:accession='PO:0009005'>  \n");
        imprimeAxeSegmentsRSML(axeCour);

        // Print the secondary axes
        axeCour1=sR->premAxe;
        while (axeCour1!=NULL)
        {
          if(axeCour1->pere==axeCour){
            int count = 0;
            segCour=axeCour1->premSeg;
            while (segCour!=NULL) {
              if (segCour->complet) count = count+1;
              segCour=segCour->suiv;
            }
            if(count > 1){
              fprintf(FSeg,"      <root ID='");
              fprintf(FSeg,"%li",axeCour1->num);
              fprintf(FSeg,"' label='root' po:accession='PO:0009005'>  \n");
              imprimeAxeSegmentsRSML(axeCour1);

                  // Print the tertiary axes
                  axeCour2=sR->premAxe;
                  while (axeCour2!=NULL)
                  {
                    if(axeCour2->pere==axeCour1){
                      int count2 = 0;
                      segCour=axeCour2->premSeg;
                      while (segCour!=NULL) {
                        if (segCour->complet) count2 = count2+1;
                        segCour=segCour->suiv;
                      }
                      if(count2 > 1){
                        fprintf(FSeg,"      <root ID='");
                        fprintf(FSeg,"%li",axeCour1->num);
                        fprintf(FSeg,"' label='root' po:accession='PO:0009005'>  \n");
                        imprimeAxeSegmentsRSML(axeCour2);

                        // Print the 4 order  axes
                          axeCour3=sR->premAxe;
                          while (axeCour3!=NULL)
                          {
                            if(axeCour3->pere == axeCour2){
                              int count3 = 0;
                              segCour=axeCour3->premSeg;
                              while (segCour!=NULL) {
                                if (segCour->complet) count3 = count3+1;
                                segCour=segCour->suiv;
                              }
                              if(count3 > 1){
                                fprintf(FSeg,"      <root ID='");
                                fprintf(FSeg,"%li",axeCour1->num);
                                fprintf(FSeg,"' label='root' po:accession='PO:0009005'>  \n");
                                imprimeAxeSegmentsRSML(axeCour3);

                                  // Print the 5 order  axes
                                  axeCour4=sR->premAxe;
                                  while (axeCour4!=NULL)
                                  {
                                    if(axeCour4->pere == axeCour3){
                                      int count4 = 0;
                                      segCour=axeCour4->premSeg;
                                      while (segCour!=NULL) {
                                        if (segCour->complet) count4 = count4+1;
                                        segCour=segCour->suiv;
                                      }
                                      if(count4 > 1){
                                        fprintf(FSeg,"      <root ID='");
                                        fprintf(FSeg,"%li",axeCour1->num);
                                        fprintf(FSeg,"' label='root' po:accession='PO:0009005'>  \n");
                                        imprimeAxeSegmentsRSML(axeCour4);
                                        fprintf(FSeg,"      </root>  \n");
                                      }
                                    }
                                    axeCour4=axeCour4->suivant;
                                  }



                                fprintf(FSeg,"      </root>  \n");
                              }
                            }
                            axeCour3=axeCour3->suivant;
                          }
                        fprintf(FSeg,"      </root>  \n");
                      }
                    }
                    axeCour2=axeCour2->suivant;
                  }
              
              fprintf(FSeg,"      </root>  \n");
            }
          }
          axeCour1=axeCour1->suivant;
        }

        fprintf(FSeg,"      </root>  \n");
      }
    }
    axeCour=axeCour->suivant;
  }
  fprintf(FSeg,"    </plant>  \n");
  fprintf(FSeg,"  </scene>\n");
  fprintf(FSeg,"</rsml>\n");
}  /* Fonction imprimeSRSegments */
/****************************************************************************/
/****************************************************************************/





void imprimeSRSegmentsEntete(void)
{   /* Imprime l'entête du fichier contenant les noeuds du syst?me racinaire */

 // fprintf(FSeg,"NumSeg Jour NumAxe Suiv Prec Diam     X1       Y1       Z1      X2       Y2       Z2\n");
  fprintf(FSeg,"NumAxe Jour Diam     X1       Y1       Z1      X2       Y2       Z2\n");

}  /* Fonction imprimeSRSegmentsEntete */
/****************************************************************************/
/****************************************************************************/
void imprimeAxeSegments(pTAxe axe)
{   /* Imprime les segments de l'axe */

pTSeg segCour;

  segCour=axe->premSeg;
  while (segCour!=NULL) {
    if (segCour->complet) imprimeSeg(segCour);
    segCour=segCour->suiv;
  }

}  /* Fonction imprimeAxeSegments */
/****************************************************************************/
/****************************************************************************/
void imprimeSRSegments(pTSysRac sR)
{  /* Imprime l'ensemble des segments du syst?me racinaire */

pTAxe axeCour;

  imprimeSRSegmentsEntete();

  axeCour=sR->premAxe;
  while (axeCour!=NULL)
  {
    imprimeAxeSegments(axeCour);
    axeCour=axeCour->suivant;
  }
}  /* Fonction imprimeSRSegments */
/****************************************************************************/
/****************************************************************************/
void calcResumeSR(pTSysRac sR)
{  /* Calcule les diff?rentes variables r?sum?es et ?criture sur fichier */

//  calcTSatisMoySR(sR);
  calcLimitesSR(sR);
//  calcVolProdSR(sR);
  translateSR(sR);
  initialiseTabSol();
  calcDistancesSR(sR);
  imprimeSRGlobal(sR);

}  /* Fonction calcResumeSR */
/****************************************************************************/
/****************************************************************************/
void fermeFichiers(void)
{
  fclose(FSeg);
  fclose(FPar);
  fclose(FSol);
  fclose(FBiom);
//  fclose(FAudit);
//  fclose(FSynth);
//  fclose(FVox);
}  /* Fonction fermeFichiers */
/****************************************************************************/

int countSegments(pTSysRac sR)
/*MODIFBEN : rend une racine (d?finie par le param?tre condemnedRoot) s?nile et arr?t?e ? un temps d?fini par le param?tre P_sacrificeTime*/
{
  pTAxe axeCour;
  int countSeg = 0;
  axeCour=sR->premAxe;
  while (axeCour!=NULL)
  {
      // Count the number of segments
      pTSeg segCour;
      segCour=axeCour->premSeg;
      while (segCour!=NULL) {
        if (segCour->complet) countSeg = countSeg+1;
        segCour=segCour->suiv;
      }
      axeCour=axeCour->suivant;
  }
  return countSeg;
}


int main(int argc, char *argv[])
{
   struct timeval tv;   // version linux

if (argc>1)
     { orig[0]=atof(argv[1]); orig[1]=atof(argv[2]); orig[2]=atof(argv[3]);}
else { orig[0]=0.0; orig[1]=0.0; orig[2]=20.0; } // semence l?g?rement enterr?e

//   gettimeofday(&tv,NULL);
//  printf("%d\n",  tv.tv_usec);
//   srand(tv.tv_usec); /* Initialisation du g?n?rateur al?atoire, version linux */

  srand( (unsigned) time(NULL) ); /* Initialisation du g?n?rateur al?atoire, version windows */

//  printf("Je passe 1  \n");

readParamXML();
ouvreFichiers();



//printf("Je passe 2  \n");


//litParam();



//printf("Je passe 3  \n");

/*
  printf(" Dur?e : %4i \n",P_duree);
  printf(" Age Maturit? Pointe : %16.5f \n",P_ageMaturitePointe);
  printf(" Vitesse d'emission des seminales : %16.5f \n",P_vitEmissionSem); // Vitesse d'?mission des primaires
  printf(" Nombre maximal de  primaires : %5i \n",P_nbMaxSem);  // Nombre maximal de primaires
  printf(" Diam?tre minimal : %16.5f \n",P_diamMin);  //  Diam?tre minimal (mm)
  printf(" Diam?tre maximal : %16.5f \n",P_diamMax);  // Diam?tre maximal (mm)
  printf(" Pente vitesse de croissance diam?tre : %16.5f \n",P_penteVitDiam);  //  Pente vitesse de croissance diam?tre
  printf(" Type de tropisme : %2i \n",P_tendanceDirTropisme); //   Type de tropisme (gravi positif, exo, plagio)
  printf(" Intensit? du tropisme : %16.5f \n",P_intensiteTropisme); //  Intensit? du tropisme
  printf(" Distance inter-ramifications : %16.5f \n",P_distRamif); //  Distance inter-ramifications (mm)
  printf(" Pente de la relation entre diam?tre des lat?rales : %16.5f \n",P_propDiamRamif); // Pente de la relation entre diam?tre des lat?rales et de la porteuse
  printf(" Coefficient de variabilit? des lat?rales : %16.5f \n",P_coeffVarDiamRamif); // Coefficient de variabilit? des lat?rales
  printf(" Masse volumique des racines : %16.5f \n",P_TMD); // Masse volumique des racines
  printf(" Pente dur?e de croissance : %16.5f \n",P_penteDureeCroissDiam2); // Pente dur?e de croissance
  printf(" Pente dur?e de vie : %16.5f \n",P_penteDureeVieDiamTMD); // Pente dur?e de vie
  printf(" Coefficient de croissance radiale : %16.5f \n",P_coeffCroissRad); // Coefficient de croissance radiale
*/

litSol();

//printf("Je passe 1  \n");

sR=initialiseSR(orig);

//printf("Je passe 12  \n");

litBiomasseMaxSR(sR);

printf("Je passe 13  \n");

// printf("Taille de TSeg: %4li \n",sizeof(TSeg));
// printf("Taille de TAxe: %4li \n",sizeof(TAxe));
// printf("Taille de TPointe: %4li \n",sizeof(TPointe));
// Sleep(5001);
int segs = 0;
while (temps<P_duree)
{
  temps=temps+deltaT;
//   printf("Temps : %3i \n",temps);

  segs = countSegments(sR);
  printf("Segments : %3i \n",segs);

  /* Emission des racines s?minales */
  emissionSemSR(sR);

  /* Emission des racines adentives */
  emissionAdvSR(sR);

  /* D?veloppement du syst?me racinaire */
  developpeSR(sR);

  /* Croissance radiale du syst?me racinaire */
  croissanceRadialeSR(sR, P_coeffCroissRad);

  /* Mortalit? du syst?me racinaire */
  mortaliteSR(sR);

//  if (temps==40) calcResumeSR(sR);       // si on veut calculer les volumes colonis?s ? certaines dates
//  if (temps==50) calcResumeSR(sR);
//  if (temps==60) calcResumeSR(sR);

}


if(P_exportType == 1){
   imprimeSRSegments(sR);   // pour imprimer le syst?me racinaire sous forme d'un ensemble de segments 
}else if(P_exportType == 2){
  imprimeSRSegmentsRSML(sR);
}

//  imprimeAudit();
delete sR;
fermeFichiers();
//Sleep(1001);  // pas possible sur linux
return 0;
}

