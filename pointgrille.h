//terme vecteur
class TermVec
{
 public:
  double cont; //continuité
  double elanX; //élan en direction X
  double elanY;
  double elanZ;
  double nrg; //énergie

  TermVec();
  TermVec(double c, double x, double y, double z, double e);
  ~TermVec();
};

// chaque point sur le quadrillage consiste en un de ces classes
//elle contient 5 objets: U, E, F, G, et J
//U - solution
//E,F,G - termes flux
//J - terme source

class PointGrille
{
 public:
  TermVec U;
  TermVec E;
  TermVec F;
  TermVec G;
  TermVec J;

  PointGrille();
  ~PointGrille();
  void initICs();
  //autre constructeur ici pour accéder à TermVec::TermVec(c,x,y,z,e) ?
};
