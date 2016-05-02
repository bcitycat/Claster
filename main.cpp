/***********************************************************/
/* Rendezetlen kvantum spinrendszerek renormalizacioja     */
/*                                                         */
/* Alapeset: RTFIC lanc, vagy letra                        */
/*                                                         */
/* Kovacs Istvan                                           */
/***********************************************************/
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <set>
#include <iostream>
#include <time.h>
#include <vector>
#include "binheap.h" 
#include "MersenneTwister.h"
#include <map>
#include <iomanip>
/****************************************************************
*Itt kezdodik a program						 *
****************************************************************/

using namespace std;                          
/*miért nyitja ki a standard névteret?
 * ha utána  std:: ot használ?
 * mely típusok definiálásakor van ez kihasználva,
 * törölhető e ez a  rész?
 */

typedef set<int, std::greater<int> > halmaz; /*Itt mi történik?*/
typedef set<int, std::less<int> > felhalmaz; 
/*Halmaz felhalmaz sehol mashol nem fordul elo*/

typedef std::map<int, int> intmap;			/*dinamikus asszociativ tipus*/
typedef std::map<int, double> doublemap; 	/*kulcs int ertek int OR double*/

typedef std::map<int, intmap> intmatrix;
typedef std::map<int, doublemap> doublematrix;


doublematrix c,crg;
/* Mit csinál c, crg változó?*/

struct edge_type 
	{int                        start, end, h_id;
	//double                     weight;
	};
/*létrehozzuk az élek struktúráját
 * Eleje, vége, és valami h_id indexe van, 
 * lehetne súlya is, de az most nins
 * 
 * mia h_id?*/
 
	
typedef std::map<int, edge_type> ellist;
ellist elek;
/*létrehozzuk az élek listáját
 * a kulcs int az érték edge_type (vagyis az él maga)
 * */
int elszam, Nodes, NumPts,NumPts2,mode;
/*Mik ezek a változók?*
 * elszam: itt egyedül van jelen az egész kódban!
 * 
 * mi a szerepük? */

typedef BinaryHeap<double, int, std::less_equal<double> >    min_heap; 
/*létrehozunk egy bináris fát,
 * std::less_equal<double> ez micsoda?
 * */
min_heap                   m_h;
min_heap::pair_type        m_he;
/*::pair_type mi ez? */
//intmap ha,hb;
intmatrix hid;

int joa,job;
double joe, norm;
doublemap hist,sorn,sor,osz;

double D=0,Sa=0, D0=0,Sa0=0;
#define INF 1E9
intmap label;
intmap a;
int ido;


double delta1(int a,int b)
	{//ha csak a sorokat vonjuk ossze, az oszlopokat nem
	double r=0,s=0,s2=0;
	/* ere a fordító azt írja, hogy nem használt változó
	 * Szerepe?*/
	
	//v: foatlot adja majd
	doublemap uj;
	//a)a kisebbik meretu soron iteralunk vegig, mindig megnezve az adott elem megvan-e a masikban is
	//b)letre kell hozni az osszeg vektort es az entropiakat kell kiszamolni
	//az osszeghez siman vegigiteralunk mindketton es egy kulon map sorba eltesszuk a kapott ertekeket
	// a vegen ennek az entropiajat is ki tudjuk szamolni
	//-S(a)
	for(doublemap::iterator i=crg[a].begin();i!=crg[a].end();i++)
		{uj[i->first]=i->second;
		if (i->second>0)
			{r+=i->second*log(i->second);
			//s+=i->second;
			}
		}
    if (sor[a]>0)
        {r-=sor[a]*log(sor[a]);
        }
	//-S(b)	
	for(doublemap::iterator i=crg[b].begin();i!=crg[b].end();i++)
		{if (uj.find(i->first)!=uj.end())
			{uj[i->first]+=i->second;
			}
		else
			{uj[i->first]=i->second;
			}
		if (i->second>0)
			{r+=i->second*log(i->second);
			//s2+=i->second;
			}
		}
    if (sor[b]>0)
        {r-=sor[b]*log(sor[b]);
        }
	
	for(doublemap::iterator i=uj.begin();i!=uj.end();i++)
		{//uj[i->first]=i->second;
		if (i->second>0)
			{r-=i->second*log(i->second);
			}
		//s+=i->second;
		}
    if ((sor[a]+sor[b])>0)
        {r+=(sor[a]+sor[b])*log(sor[a]+sor[b]);
        }
	
	return (r);	
	};

void fuse1()
	{//csak a sorokat vonjuk ossze, az oszlopokat nem
	
	for(doublemap::iterator i=crg[joa].begin();i!=crg[joa].end();i++)
		{crg[Nodes][i->first]=i->second;
		}
	for(doublemap::iterator i=crg[job].begin();i!=crg[job].end();i++)
		{if (crg[Nodes].find(i->first)!=crg[Nodes].end())
			{crg[Nodes][i->first]+=i->second;
			}
		else
			{crg[Nodes][i->first]=i->second;
			}
		}	
	
	crg.erase(joa);
	crg.erase(job);
	
	//hist-et se hasznaljuk:
	hist[Nodes]=joe;
	sor[Nodes]=sor[joa]+sor[job];
	//ez igy nem biztos, hogy korrekt, de jelenleg nem is hasznaljuk
	sorn[Nodes]=sor[joa]+sor[job];
	//S(a) update:
	Sa+=sor[joa]*log(sor[joa])/norm;
	Sa+=sor[job]*log(sor[job])/norm;
	Sa-=sor[Nodes]*log(sor[Nodes])/norm;
    //eredeti D update:
  //  Di+=2*joe/norm;
    //    Di+=sor[joa]*log(sor[joa])/norm;
      //  Di+=sor[job]*log(sor[job])/norm;
      //  Di-=sor[Nodes]*log(sor[Nodes])/norm;
	//D update
	D+=joe/norm;
	label[joa]=Nodes;
	label[job]=Nodes;
	label[Nodes]=-1;
	Nodes++;
	};

void kiszed()
	{m_he=m_h.pop();
	//first v second?
	joe = m_he.first;
	joa=m_he.second/(2*NumPts);
	job=m_he.second%(2*NumPts);
	//mar nem kellenek, kitorolhetjuk
	//Ki kell szedni az osszes kapcsolodot:
	for(doublematrix::iterator i=crg.begin(); i!=crg.end();i++ )
		{if (i->first!=joa)
			{if (i->first!=job)
				{if (i->first<joa)
					{m_h.pop(hid[i->first][joa]);
					hid[i->first].erase(joa);
					}
				else
					{m_h.pop(hid[joa][i->first]);
					}
				if (i->first<job)
					{m_h.pop(hid[i->first][job]);
					hid[i->first].erase(job);
					}
				else
					{m_h.pop(hid[job][i->first]);
					}
				}					
			}
		}
	hid.erase(joa);
	hid.erase(job);
	}
	
void update_in()
	{//kezdeti szamolas: heap feltoltese
	//de be kell hozni a heap-et az optimum kivalasztasahoz :(	
	//joe=INF;
	double jobb;
	for(doublematrix::iterator i=crg.begin(); i!=crg.end();i++ )
		{//tmp_el.start=i;
           // cout << i->first << "\n";
		doublematrix::iterator j=i;
		j++;	
		//joe=INF;
		for(; j!=crg.end();j++ )
			{//if ((i->first)<(j->first))
			//jobb=delta(i->first,j->first)+hist[i->first]+hist[j->first];	
			//hid[i->first][j->first]=m_h.push(delta(i->first,j->first)+hist[i->first]+hist[j->first],i->first*2*NumPts+j->first);
			jobb=delta1(i->first,j->first);
               // cout << jobb << "\n";
			hid[i->first][j->first]=m_h.push(jobb,i->first*2*NumPts+j->first);
			}
		}
	kiszed();
	};

void update2(int be)
	{//csak az uj el ertekeit kell updatelni
	//be adja meg, miaz uj index
	//bar ezt igazabol Nodes rogziti!, valoszinuleg Nodes-1 az erteke	
	//de be kell hozni a heap-et az optimum kivalasztasahoz :(	
	//joe=INF;
        
	double jobb;
	for(doublematrix::iterator i=crg.begin(); i!=crg.end();i++ )
		{jobb=delta1(be,i->first);	
		if (be<i->first)
			{hid[be][i->first]=m_h.push(jobb,be*2*NumPts+i->first);
			}
		else
			{if (be!=i->first)
				{hid[i->first][be]=m_h.push(jobb,i->first*2*NumPts+be);
				}
			}	
		}
	kiszed();
	};

int find(int x)
{//negativ: meretet adja, nemnegativ: azt hogy hova olvadt
    if (label[x]>-1)
    {//x=nodes[x];
        //if (nodes[nodes[x]]<0) {return nodes[x];}
        int y=label[x];
        int z;
        while (y>-1)
        {z=y;
            y = label[z];
        }
        //z mutatja a jo root pontidt
        //u az eredeti pontot
        while (label[x]!=z)
        {y = label[x];
            label[x]=z;
            x=y;
        }
        return z;
    }
    else
    {return x;
    }
}

void clu_out( )
{intmap szin;
    int ii;
    int ij=1;
    char llag[255];
    sprintf(llag,"l_%d_%d_%d_%lf.clu",ido, NumPts,  2*NumPts-Nodes+1,D);
    ofstream ll(llag, ios::out);
    if (!ll)
    {cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
        exit(-1);
    };
    //cout <<D << "\n";
    for(int i=1;i<=NumPts;i++)
    {ii=find(i);
        if (szin.find(ii)==szin.end())
        {szin[ii]=ij;
            ij++;
        }
        ll << szin[ii]<< "\n";
    }
}

void crg_out( )
{char llag[255];
    sprintf(llag,"crg_%d_%d_%d_%lf.dat",ido, NumPts,  2*NumPts-Nodes+1,D);
    ofstream ll(llag, ios::out);
    if (!ll)
    {cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
        exit(-1);
    };
    for(doublematrix::iterator id=crg.begin();id!=crg.end();id++)
    {for(int j=1; j<=NumPts2;j++ )
    {if (crg.find(j)!=crg.end())
				{if (id->second.find(j)!=id->second.end())
                {ll << id->second[j]/norm <<"\11";
                }
                else
                {ll << 0 <<"\11";
                }
                }
    }
        ll << "\n";
    }	
    
}


void RG(ofstream *flo,ofstream *flm, ofstream *fll  )
	{//int elid;
	//double x;
	//double sum=0;	
	crg=c;
	Nodes=NumPts+1;
	for(int i=1;i<Nodes;i++)
		{hist[i]=0;
	label[i]=-1;
		}
		*flm <<D << "\11"<<Sa<<"\11"<< D/D0-1+(D-D0+Sa)/(Sa0-D0) << "\11"<< D+D-D0+Sa<< "\n";
       // cout <<D << "\11"<<Sa<<"\11"<< D/D0-1+(D-D0+Sa)/(Sa0-D0) << "\11"<< D+D-D0+Sa<< "\n";
	//init();
	//cout << "Init\n";
	//while (!m_h.isEmpty()) 
	while (Nodes<NumPts*2)
		{//m_he = m_h.pop(hid[2][3]);
		//update(Nodes-NumPts);
		
		if (Nodes==NumPts+1)
			{update_in();
			}
		else	
			{update2(Nodes-1);
			}
            
		//itt kell kiirni azt is, az epp megszuntetett sor eseten a naiv onatfedes hanyadresze az igazi, amiatt hogy az eredeti elek nem feltetlenul fednek at
cout <<joe/norm << "\11"<<Nodes<<"\11"<<joa<<"\11"<<job<<"\11"<<"\n";
	*flo <<joe/norm << "\11"<<Nodes<<"\11"<<joa<<"\11"<<1.0-(1-mode)*sorn[joa]/(sor[joa]*sor[joa])<<"\11"<<crg[joa].size()<<"\11";
	
		//sum+=joe;
		//kiirjuk a ket erintett sort crg-bol
		//*flo << joa<<"\11";
		//*flo << joe<<"\11";
		for(doublemap::iterator id=crg[joa].begin();id!=crg[joa].end();id++)
			{//kiirjuk melyik elem letezik es mekkora
			//if(id->second>0)
				{*flo << id->first<<"\11"<< id->second/norm<<"\11";
				}
			}
		//*flo <<"\n";
		*flo << "\n" <<joe/norm << "\11"<<Nodes<<"\11"<<job<<"\11"<<1.0-(1-mode)*sorn[job]/(sor[job]*sor[job])<<"\11"<<crg[job].size()<<"\11";
		//*flo << job<<"\11";
		//*flo << joe<<"\11";
		for(doublemap::iterator id=crg[job].begin();id!=crg[job].end();id++)
			{//kiirjuk melyik elem letezik es mekkora
			//if(id->second>0)
				{*flo << id->first<<"\11"<< id->second/norm<<"\11";
				}
			}
		*flo <<"\n";
		//osszeolvasztjuk
		//cout << Nodes<< "\11";
		fuse1();
	
            //clu_out();
            
            if (2*NumPts-Nodes<60)
            {clu_out();
            crg_out();
            }
		*flm <<D << "\11"<<Sa<<"\11"<< D/D0-1+(D-D0+Sa)/(Sa0-D0) << "\11"<< D+D-D0+Sa<<"\n";
		}
	//a legvegen kiirjuk az utolso pontot is:
	*flo <<Nodes-1<<"\11"<<1.0-(1-mode)*sorn[Nodes-1]/(sor[Nodes-1]*sor[Nodes-1])<<"\n";
	};

	
void qinit_v1()
{//ha a rat-ok 1-ek, akkor siman a mutual info a kezdeti ertek
	int i,j;
	double q=0;
	//q kiszamitasa:
	cout << "Norm:" << norm  << "\n";
	for (i = 1; i <= NumPts;i++)
		{if (sor[i]>0)
			{Sa0-=sor[i]/norm*log(sor[i]/norm);
		q-=sor[i]/norm*log(sor[i]/norm);
		for (j = 1; j <= NumPts2;j++)
			{//ha letezik c-ben:
			//if (j!=i)
			{if (c[i].find(j)!=c[i].end())
				{if (c[i][j]>0)
					{q+=c[i][j]/norm*log(c[i][j]/norm);
					}
				}
			}
			/*else
				{if (c[i][i]>0)
					{q+=c[i][i]/norm*log(c[i][i]/norm);
					}
				}
             */
			}
			}
		}
    for (j = 1; j <= NumPts2;j++)
    {q-=osz[j]/norm*log(osz[j]/norm);
    }
	D0=q;
	//D=D0;	
	Sa=Sa0;
    //Di=0;
    D=0;
	cout << "Q:" <<q<< "\11";	
	//return (q);
	}

int main(int argc, char** argv)
/*argv tehát egy char pointerre mutato pointer*/
{	
	//beolvasas
	ido=(int) time(NULL);
	int nov=5;
	int type =atoi(argv[1]);  //beolvasasi mod
		/*atoi ez függvény?
		* hol van definiálva? mit csinál, miért egy tömb a paramétere?
		* */
	mode =atoi(argv[2]);  //van-e foatlo: 0 ha nincs
	NumPts = atoi(argv[3]);   //sorok szama
	NumPts2 = atoi(argv[4]);   //oszlopok szama	
	char *fg;
	fg=argv[5];
	FILE *lg;
	 
	if(!(lg=fopen(fg,"r"))) 
	{	
		fprintf(stderr,"Bemeneti h fajl megnyitasa sikertelen!\n");
		exit(-1);
	} 
	
	int seed=0; 
	seed = atoi(argv[6]);
	MTRand mtrand2(seed);
	double s;
	int j=0; 
	int e1,e2;
	norm=0;
	
	int volte;	
	if (type==1)
	{//ellistakent jon be c
		for(int i=1; i<=NumPts; i++)
        {	//c[i][i]=0;
            sor[i]=0;
            sorn[i]=0;
            osz[i]=0;
        }
        while (!feof(lg))
		{
			j++;
			fscanf(lg,"%d",&e1);
			fscanf(lg,"%d",&e2);
			fscanf(lg,"%lf",&s);
			if (s>0)
			{
				//ha meg nem volt ezt az el:
				volte=0;	
				if (c.find(e1)!=c.end())	
				{
					if (c[e1].find(e2)!=c[e1].end())
					{
						volte=1;
					}
				}
				if (volte==0)	
				{
					c[e1][e2]=s;
					c[e2][e1]=s;
					sor[e1]+=s;
					osz[e2]+=s;
					if (e1!=e2)
					{
						sor[e2]+=s;
						osz[e1]+=s;
						//a foatlot is feltoltjuk:
						if (mode==1)
						{
							sor[e1]+=s;
							sor[e2]+=s;
							c[e1][e1]+=s;
							c[e2][e2]+=s;
							norm+=2*s;
						}
						norm+=2*s;
					}
					else
					{
						//sc[e1]+=2*s;
						//sc[e2]+=2*s;
						norm+=s;
					}
				}
			}
		fscanf(lg,"\n");
		}
		fclose(lg);		
	}
	else
	{
		if (type==0)
		{
			//matrixkent jon be c
			for(int i=1; i<=NumPts; i++)
			{	//c[i][i]=0;
				sor[i]=0;
				sorn[i]=0;
			}
        for(int i=1; i<=NumPts2; i++)
        {//c[i][i]=0;
            osz[i]=0;
            //sorn[i]=0;
        }
        j=0;
	while (!feof(lg))
		{j++;
		for(int i=1; i<=NumPts2; i++)
			{fscanf(lg,"%lf",&s);
			if (s>0)
				{//c[j][i]=s;
				c[j][i]=s;
				sor[j]+=s;
                osz[i]+=s;
				//sor[e2]+=s;
				//sc[i]+=s;
				norm+=s;
				}
			}
		fscanf(lg,"\n");
		}
	fclose(lg);	
	}
        else
        {//ellista, de nem trukkozunk a foatloval, es aszimmetrikus
            //ellistakent jon be c
            for(int i=1; i<=NumPts; i++)
            {//c[i][i]=0;
                sor[i]=0;
                //sorn[i]=0;
            }
            for(int i=1; i<=NumPts2; i++)
            {//c[i][i]=0;
             //   sor[i]=0;
                osz[i]=0;
            }
            j=0;
            while (!feof(lg))
            {j++;
                fscanf(lg,"%d",&e1);
                fscanf(lg,"%d",&e2);
                fscanf(lg,"%lf",&s);
                if (s>0)
                {//ha meg nem volt ezt az el:
                    //volte=0;
                    {c[j][e1]=s;
                    c[j][e2]=s;
                        sor[j]+=2*s;
                        //sor[e2]+=s;
                        osz[e1]+=s;
                        osz[e2]+=s;
                        //a foatlot is feltoltjuk:
                        //sc[e1]+=2*s;
                        //sc[e2]+=2*s;
                        norm+=2*s;
                    }
                }
                fscanf(lg,"\n");
            }
            fclose(lg);
            
            
        }
    }
	
	int fend=0;
	char random_seed; 
	if(argc>nov)
		{fend = atoi(argv[nov]);
	//cout << "OUT";
		}
	else
		{FILE *random_file = fopen("/dev/urandom","r");
	//NEEDS ERROR CHECKING
	if (!random_file)
		{cout << "Zaj fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
	random_seed = getc(random_file);
	fend=random_seed;
		}
	
	MTRand mtrand(fend*ido);
	//MTRand mtrand(fend);
	
	char mag[255];	
	sprintf(mag,"mag_%d_%d_%d.dat",fend,ido, NumPts);
	ofstream lmag(mag, ios::out);
	if (!lmag)
		{cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
		exit(-1);
		};	
	char mmag[255];	
	sprintf(mmag,"mm_%d_%d_%d.dat",fend,ido, NumPts);
	ofstream lmm(mmag, ios::out);
	if (!lmm)
		{cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
		exit(-1);
		};	
	char lag[255];	
	sprintf(lag,"l_%d_%d_%d.dat",fend,ido, NumPts);
	ofstream ll(lag, ios::out);
	if (!ll)
		{cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
		exit(-1);
		};	
	
	//if (mode==1)
		{
	for(int i=1; i<=NumPts; i++)
		{sorn[i]=sor[i]*sor[i];
		}
		}
	char ma[255];	
	sprintf(ma,"init_c_%d.dat", NumPts);
	ofstream lma(ma, ios::out);
	if (!lma)
		{cout << "Kimeneti c fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
	char ca[255];	
	sprintf(ca,"init_a_%d.dat", NumPts);
	ofstream cma(ca, ios::out);
	if (!cma)
		{cout << "Kimeneti a fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
	
	for(int i=1; i<=NumPts;i++ )
		{for(int j=1; j<=NumPts2;j++ )
			{if (c[i].find(j)!=c[i].end())	
				{cma << c[i][j]/norm <<"\11";
				}
			else
				{cma << 0 <<"\11";
				}
			}
		//lmag << "\n";	
		cma << "\n";
		}
	
	for(int i=1; i<=NumPts;i++ )
		{lma << c[i].size() <<"\11";
		for(doublemap::iterator j=c[i].begin(); j!=c[i].end();j++ )
			{lma << j->first<< "\11" <<j->second/norm <<"\11";
			}
		lma << "\n";
		}
	
	lma.close();
	cout << "Qinit:\n";
	qinit_v1();
	cout << "RG\n";
	RG(&lmag, &lmm, &ll);
	
	//nodes.clear();
	//edges.clear();
	
	return 0;
}
