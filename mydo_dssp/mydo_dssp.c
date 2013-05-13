#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "string2.h"
#include "strdb.h"
#include "macros.h"
#include "smalloc.h"
#include "mshift.h"
#include "statutil.h"
#include "copyrite.h"
#include "pdbio.h"
#include "gmx_fatal.h"
#include "xvgr.h"
#include "matio.h"
#include "index.h"
#include "gstat.h"
#include "tpxio.h"
#include "viewit.h"

static void strip_dssp(char *dsspfile,int nres,
		       bool bPhobres[],real t,
		       real *acc,FILE *fTArea,
		       t_matrix *mat,int average_area[])
{
  static bool bFirst=TRUE;
  static char *ssbuf;
  FILE *tapeout;
  static int xsize,frame;
  char buf[STRLEN+1];
  char SSTP;
  int  i,nr,iacc;
  real iaccf,iaccb;
  t_xpmelmt c;
  
  tapeout=ffopen(dsspfile,"r");
  
  /* Skip header */
  do {
    fgets2(buf,STRLEN,tapeout);
  } while (strstr(buf,"KAPPA") == NULL);
  if (bFirst) {
    snew(ssbuf,nres+10);
  }
  
  iaccb=iaccf=0;
  for(nr=0; (fgets2(buf,STRLEN,tapeout) != NULL); nr++) {
    SSTP=buf[16]==' ' ? '~' : buf[16];
    ssbuf[nr]=SSTP;
    
    buf[39]='\0';
    sscanf(&(buf[34]),"%d",&iacc);
    acc[nr]=iacc;
    average_area[nr]+=iacc;
    if (bPhobres[nr])
      iaccb+=iacc;
    else
      iaccf+=iacc;
  }
  ssbuf[nr]='\0';
  
  if (bFirst) {
    sprintf(mat->title,"Secondary structure");
    mat->legend[0]=0;
    sprintf(mat->label_x,"%s",time_label());
    sprintf(mat->label_y,"Residue");
    mat->bDiscrete=TRUE;
    mat->ny=nr;
    snew(mat->axis_y,nr);
    for(i=0; i<nr; i++)
      mat->axis_y[i]=i+1;
    mat->axis_x=NULL;
    mat->matrix=NULL;
    xsize=0;
    frame=0;
    bFirst=FALSE;
  }
  if (frame>=xsize) {
    xsize+=10;
    srenew(mat->axis_x,xsize);
    srenew(mat->matrix,xsize);
  }
  mat->axis_x[frame]=t;
  snew(mat->matrix[frame],nr);
  c.c2=0;
  for(i=0; i<nr; i++) {
    c.c1=ssbuf[i];
    mat->matrix[frame][i] = max(0,searchcmap(mat->nmap,mat->map,c));
  }
  frame++;
  mat->nx=frame;
  
  if (fTArea)
    fprintf(fTArea,"%10g  %10g  %10g\n",t,0.01*iaccb,0.01*iaccf);
  fclose(tapeout);
}

bool *bPhobics(t_atoms *atoms)
{
  int  i,nb;
  char **cb;
  bool *bb;
  
  nb=get_strings("phbres.dat",&cb);
  snew(bb,atoms->nres);
  
  for(i=0; (i<atoms->nres); i++) {
    if (search_str(nb,cb,*atoms->resname[i]) != -1)
      bb[i]=TRUE;
  }
  return bb;
}
 
static void check_oo(t_atoms *atoms)
{
  static char *OOO="O";
  int i;
  
  for(i=0; (i<atoms->nr); i++) {
    if (strcmp(*(atoms->atomname[i]),"OXT")==0)
      atoms->atomname[i]=&OOO;
    else if (strcmp(*(atoms->atomname[i]),"O1")==0)
      atoms->atomname[i]=&OOO;
  }
}

static void norm_acc(t_atoms *atoms, int nres, 
		     real av_area[], real norm_av_area[])
{
  int i,n,n_surf;
  
  char surffn[]="surface.dat";
  char **surf_res, **surf_lines;
  double *surf;
  
  n_surf = get_lines(surffn, &surf_lines);
  snew(surf, n_surf);
  snew(surf_res, n_surf);
  for(i=0; (i<n_surf); i++) {
    snew(surf_res[i], 5);
    sscanf(surf_lines[i],"%s %lf",surf_res[i],&surf[i]);
  }
  
  for(i=0; (i<nres); i++) {
    n = search_str(n_surf,surf_res,*atoms->resname[i]);
    if ( n != -1)
      norm_av_area[i] = av_area[i] / surf[n];
    else 
      fprintf(stderr,"Residue %s not found in surface database (%s)\n",
	      *atoms->resname[i],surffn);
  }
}

void prune_ss_legend(t_matrix *mat)
{
  bool *present;
  int  *newnum;
  int i,r,f,newnmap;
  t_mapping *newmap;
  
  snew(present,mat->nmap);
  snew(newnum,mat->nmap);

  for(f=0; f<mat->nx; f++)
    for(r=0; r<mat->ny; r++)
      present[mat->matrix[f][r]]=TRUE;
      
  newnmap=0;
  newmap=NULL;
  for (i=0; i<mat->nmap; i++) {
    newnum[i]=-1;
    if (present[i]) {
      newnum[i]=newnmap;
      newnmap++;
      srenew(newmap,newnmap);
      newmap[newnmap-1]=mat->map[i];
    }
  }
  if (newnmap!=mat->nmap) {
    mat->nmap=newnmap;
    mat->map=newmap;
    for(f=0; f<mat->nx; f++)
      for(r=0; r<mat->ny; r++)
	mat->matrix[f][r]=newnum[mat->matrix[f][r]];
  }
}

void write_sas_mat(char *fn,real **accr,int nframe,int nres,t_matrix *mat)
{
  real lo,hi;
  int i,j,nlev;
  t_rgb rlo={1,1,1}, rhi={0,0,0};
  FILE *fp;
  
  if(fn) {
    hi=lo=accr[0][0];
    for(i=0; i<nframe; i++)
      for(j=0; j<nres; j++) {
	lo=min(lo,accr[i][j]);
	hi=max(hi,accr[i][j]);
      }
    fp=ffopen(fn,"w");
    nlev=hi-lo+1;
    write_xpm(fp,0,"Solvent Accessible Surface","Surface (A^2)",
	      "Time","Residue Index",nframe,nres,
	      mat->axis_x,mat->axis_y,accr,lo,hi,rlo,rhi,&nlev);
    ffclose(fp);
  }
}

void analyse_ss(char *outfile, t_matrix *mat, char *ss_string)
{
  FILE *fp;
  t_mapping *map;
  int s,f,r,*count,ss_count;
  char **leg;
  
  map=mat->map;
  snew(count,mat->nmap);
  snew(leg,mat->nmap+1);
  leg[0]="Structure";
  for(s=0; s<mat->nmap; s++)
    leg[s+1]=map[s].desc;
  
  fp=xvgropen(outfile,"Secondary Structure",
	      xvgr_tlabel(),"Number of Residues");
  if (bPrintXvgrCodes())
    fprintf(fp,"@ subtitle \"Structure = ");
  for(s=0; s<strlen(ss_string); s++) {
    if (s>0)
      fprintf(fp," + ");
    for(f=0; f<mat->nmap; f++)
      if (ss_string[s]==map[f].code.c1)
	fprintf(fp,"%s",map[f].desc);
  }
  fprintf(fp,"\"\n");
  xvgr_legend(fp,mat->nmap+1,leg);
  
  for(f=0; f<mat->nx; f++) {
    ss_count=0;
    for(s=0; s<mat->nmap; s++)
      count[s]=0;
    for(r=0; r<mat->ny; r++)
      count[mat->matrix[f][r]]++;
    for(s=0; s<mat->nmap; s++) {
      if (strchr(ss_string,map[s].code.c1))
	ss_count+=count[s];
    }
    fprintf(fp,"%8g %5d",mat->axis_x[f],ss_count);
    for(s=0; s<mat->nmap; s++)
      fprintf(fp," %5d",count[s]);
    fprintf(fp,"\n");
  }
  
  fclose(fp);
  sfree(leg);
  sfree(count);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "mydo_dssp, modified from do_dssp by zyxue\n",
    "you need\n    export DSSP=/opt/dssp/bin/dssp\nor something like that\n",
    "only -sc options is working, other options are disabled."
  };
  static bool *bVerbose=FALSE;
  static char *ss_string="HEBT"; 

  t_pargs pa[] = {
    { "-v", FALSE, etBOOL, {&bVerbose},
      /* "HIDDENGenerate miles of useless information" }, */
      "Generate miles of useless information" },
    { "-sss", FALSE, etSTR, {&ss_string},
      "Secondary structures for structure count"},
  };
  
  int        status;
  FILE       *tapein;
  FILE       *ss,*acc,*fTArea,*tmpf;
  char       *fnSCount,*fnArea,*fnTArea,*fnAArea;
  char       *leg[] = { "Phobic", "Phylic" };
  t_topology top;
  int        ePBC;
  t_atoms    *atoms;
  t_matrix   mat;
  int        nres,nr0,naccr;
  bool       *bPhbres,bDoAccSurf;
  real       t;
  int        i,j,natoms,nframe=0;
  matrix     box;
  int        gnx;
  char       *grpnm,*ss_str;
  atom_id    *index;
  rvec       *xp,*x;
  int        *average_area;
  real       **accr,*av_area, *norm_av_area;
  /* char       pdbfile[32],tmpfile[32],title[256]; */
  char       pdbfile[128],tmpfile[128],title[256];
  char       dssp[256],*dptr;
  
  t_filenm   fnm[] = {
    { efTRX, "-f",   NULL,      ffREAD },
    { efTPS, NULL,   NULL,      ffREAD },
    { efNDX, NULL,   NULL,      ffOPTRD },
    { efDAT, "-ssdump", "ssdump", ffOPTWR },
    { efMAP, "-map", "ss",      ffLIBRD },
    { efXPM, "-o",   "ss",      ffWRITE },
    { efXVG, "-sc",  "scount",  ffWRITE },
    { efXPM, "-a",   "area",    ffOPTWR },
    { efXVG, "-ta",  "totarea", ffOPTWR },
    { efXVG, "-aa",  "averarea",ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_TIME_UNIT | PCA_BE_NICE ,
		    NFILE,fnm, asize(pa),pa, asize(desc),desc,0,NULL);

  /* added by zyxue */
  char *ftrx;
  ftrx = opt2fn("-f", NFILE, fnm);
  /* printf("%s\n", "########################################"); */
  /* printf("%s\n", ftrx); */
  /* printf("%s\n", "########################################"); */

  fnSCount= opt2fn("-sc",NFILE,fnm);
  fnArea  = opt2fn_null("-a", NFILE,fnm);
  fnTArea = opt2fn_null("-ta",NFILE,fnm);
  fnAArea = opt2fn_null("-aa",NFILE,fnm);
  bDoAccSurf=(fnArea || fnTArea || fnAArea);
  
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xp,NULL,box,FALSE);
  atoms=&(top.atoms);
  check_oo(atoms);
  bPhbres=bPhobics(atoms);
  
  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpnm);
  nres=0;
  nr0=-1;
  for(i=0; (i<gnx); i++) {
    if (atoms->atom[index[i]].resnr != nr0) {
      nr0=atoms->atom[index[i]].resnr;
      nres++;
    }
  }
  fprintf(stderr,"There are %d residues in your selected group\n",nres);

  /* too avoid name collision when multiple processes of dssp are running at the
     same time */
  strcpy(pdbfile, ftrx);
  strcat(pdbfile, ".zyxue.pdbfile.ddXXXXXX");

  int fd;
  printf("%s\n", "NO_gmx_tmpnap");
  fd = mkstemp(pdbfile);
  close(fd);

  /* gmx_tmpnam(pdbfile); */
  if ((tmpf = fopen(pdbfile,"w")) == NULL) {
    sprintf(pdbfile,"%ctmp%cfilterXXXXXX",DIR_SEPARATOR,DIR_SEPARATOR);
    gmx_tmpnam(pdbfile);
    if ((tmpf = fopen(pdbfile,"w")) == NULL) 
      gmx_fatal(FARGS,"Can not open tmp file %s",pdbfile);
  }
  else
    fclose(tmpf);
    
  strcpy(tmpfile, ftrx);
  strcat(tmpfile, ".zyxue.tmpfile.ddXXXXXX");

  fd = mkstemp(tmpfile);
  close(fd);

  /* gmx_tmpnam(tmpfile); */
  if ((tmpf = fopen(tmpfile,"w")) == NULL) {
    sprintf(tmpfile,"%ctmp%cfilterXXXXXX",DIR_SEPARATOR,DIR_SEPARATOR);
    gmx_tmpnam(tmpfile);
    if ((tmpf = fopen(tmpfile,"w")) == NULL) 
      gmx_fatal(FARGS,"Can not open tmp file %s",tmpfile);
  }
  else
    fclose(tmpf);

  /* upon here, pdbfile and tmpfile have already been created, but dssp program
     won't do backup while ffopen function does, that's why pdbfile always get
     back up once  */
  printf("%s\n", "########################################");
  printf("pdbfile: %s\n", pdbfile);
  printf("tmpfile: %s\n", tmpfile);
  char pdbfile_cmd[100], tmpfile_cmd[100];
  sprintf(pdbfile_cmd, "ls -l %s", pdbfile);
  sprintf(tmpfile_cmd, "ls -l %s", tmpfile);
  system(pdbfile_cmd);
  system(tmpfile_cmd);
  printf("%s\n", "########################################");
  remove(pdbfile);
  remove(tmpfile);

  if ((dptr=getenv("DSSP")) == NULL)
    dptr="/usr/local/bin/dssp";
  if (!fexist(dptr))
    gmx_fatal(FARGS,"DSSP executable (%s) does not exist (use setenv DSSP)",
		dptr);
  sprintf(dssp,"%s %s %s %s > /dev/null %s",
	  dptr,bDoAccSurf?"":"-na",pdbfile,tmpfile,bVerbose?"":"2> /dev/null");
  if (bVerbose)
    fprintf(stderr,"dssp cmd='%s'\n",dssp);

  mat.map=NULL;
  mat.nmap=getcmap(libopen(opt2fn("-map",NFILE,fnm)),
		   opt2fn("-map",NFILE,fnm),&(mat.map));
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if (natoms > atoms->nr) 
    gmx_fatal(FARGS,"\nTrajectory does not match topology!");
  if (gnx > natoms)
    gmx_fatal(FARGS,"\nTrajectory does not match selected group!");
  
  snew(average_area,atoms->nres+10); /* useless to me, but will be used in the following code */
  snew(av_area,atoms->nres+10);
  snew(norm_av_area,atoms->nres+10);

  accr=NULL;
  naccr=0;
  do {
    t = convert_time(t);
    if (nframe>=naccr) {
      naccr+=10;
      srenew(accr,naccr);
      for(i=naccr-10; i<naccr; i++)
	snew(accr[i],atoms->nres+10);
    }
    rm_pbc(&(top.idef),ePBC,natoms,box,x,x);
    tapein=ffopen(pdbfile,"w");

    /* zyxue */
    /* sprintf(pdbfile_cmd, "ls -l %s", pdbfile); */
    /* system(pdbfile_cmd); */
    /* zyxue */

    write_pdbfile_indexed(tapein,NULL,atoms,x,ePBC,box,0,-1,gnx,index);
    fclose(tapein);

#ifdef GMX_NO_SYSTEM
    printf("Warning-- No calls to system(3) supported on this platform.");
    printf("Warning-- Skipping execution of 'system(\"%s\")'.", dssp);
    exit(1);
#else
    if(0 != system(dssp))
    {
	gmx_fatal(FARGS,"Failed to execute command: %s",dssp);
    }
#endif

    strip_dssp(tmpfile,nres,bPhbres,t,
	       accr[nframe],fTArea,&mat,average_area);

    /* zyxue */
    /* sprintf(tmpfile_cmd, "ls -l %s", tmpfile); */
    /* system(tmpfile_cmd); */
    /* zyxue */

    remove(tmpfile);
    remove(pdbfile);
    nframe++;
  } while(read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);
  
  prune_ss_legend(&mat);

  analyse_ss(fnSCount,&mat,ss_string);

  view_all(NFILE, fnm);

  /* thanx(stderr); */
  printf("# Instead of print thanx, PLEASE be noticed this has been modified by zyxue!\n");

  return 0;
}
