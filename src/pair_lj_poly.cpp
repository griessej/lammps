/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Jan Griesser
   Reference: Lerner, E. (2019), Journal of Non-Crystalline Solids, 522, 119570.
------------------------------------------------------------------------- */

#include "pair_lj_poly.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLjPoly::PairLjPoly(LAMMPS *lmp) : Pair(lmp) {
  writedata = 1;
  centroidstressflag = 1;
}

/* ---------------------------------------------------------------------- */

PairLjPoly::~PairLjPoly()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(epsilon);
  }
}

/* ---------------------------------------------------------------------- */

void PairLjPoly::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double isize,jsize,ijsize2,ijsize10,ijsize2inv,ijsize4inv,ijsize6inv;
  double dipl,ddipl,dddipl;
  double rsq,r4,r6,r2inv,r10inv,forceipl;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *size = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    isize = size[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      jsize = size[j];
      ijsize2 = pow(0.5*(isize + jsize)*(1 - 0.1*std::abs(isize-jsize)),2.0);

      if (rsq <= xcsq*ijsize2) {
         r4 = rsq*rsq;
         r6 = r4*rsq;
         r2inv = 1.0/rsq;
         r10inv = 1.0/(r4*r6);
         ijsize10 = pow(ijsize2,5.0);
         ijsize2inv = 1.0/ijsize2; 
         ijsize4inv = ijsize2inv*ijsize2inv;
         ijsize6inv = ijsize4inv*ijsize2inv;

         // 
         dipl = c2*rsq*ijsize2inv;
         ddipl = c4*r4*ijsize4inv;
         dddipl = c6*r6*ijsize6inv;
         forceipl = 2*epsilon[itype][jtype]*(-5*ijsize10*r10inv + dipl + 2*ddipl + 3*dddipl);
         fpair = -forceipl*r2inv;

         f[i][0] += delx*fpair;
         f[i][1] += dely*fpair;
         f[i][2] += delz*fpair;
         if (newton_pair || j < nlocal) {
           f[j][0] -= delx*fpair;
           f[j][1] -= dely*fpair;
           f[j][2] -= delz*fpair;
         }   
        
         if (eflag){
            evdwl = epsilon[itype][jtype]*(ijsize10*r10inv + c0 + dipl + ddipl + dddipl);
         }     
         
         if (evflag){
           ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
         }

      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLjPoly::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLjPoly::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLjPoly::coeff(int narg, char **arg)
{
  if (narg != 3 && narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);

  double cut_one = cut_global;
  if (narg == 4) {
    cut_one = utils::numeric(FLERR,arg[3],false,lmp);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLjPoly::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/poly requires atom attribute q. Notice that q contains the size of the particle!");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLjPoly::init_one(int i, int j)
{
  if (setflag[i][j] == 0){
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut[j][i] = cut[i][j];
  epsilon[j][i] = epsilon[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLjPoly::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]){
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLjPoly::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0){
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error);  
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLjPoly::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLjPoly::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLjPoly::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],cut[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLjPoly::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g\n",i,j,
              epsilon[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLjPoly::single(int i, int j, int /*itype*/, int /*jtype*/,
                           double rsq, double &fforce)
{
  double r2inv,rinv,r4,r6,r10inv,forceipl,evdwl;
  double ijsize2,ijsize2inv,ijsize10,ijsize4inv,ijsize6inv;
  double dipl,ddipl,dddipl;

  // coefficients for smooting q=3
  double c0 = -1.9360103298820364;
  double c2 = 2.4694009309719855;
  double c4 = -1.079912943573755;
  double c6 = 0.1607013308889516;
  double xcsq = 1.96;
  //
  ijsize2 = pow(0.5*(atom->q[i] + atom->q[j])*(1 - 0.1*std::abs(atom->q[i]-atom->q[j])),2.0);
  
  if (rsq <= xcsq*ijsize2) {
    r4 = rsq*rsq;
    r6 = r4*rsq;
    rinv = sqrt(r2inv);
    r2inv = 1.0/rsq;
    r10inv = 1.0/(r4*r6);
    ijsize10 = pow(ijsize2,5.0);
    ijsize2inv = 1.0/ijsize2; 
    ijsize4inv = ijsize2inv*ijsize2inv;
    ijsize6inv = ijsize4inv*ijsize2inv;
    dipl = c2*rsq*ijsize2inv;
    ddipl = c4*r4*ijsize4inv;
    dddipl = c6*r6*ijsize6inv;

    forceipl = 2*epsilon[atom->type[i]][atom->type[j]]*(-5*ijsize10*r10inv + dipl + 2*ddipl + 3*dddipl);
    fforce = forceipl * r2inv;
    evdwl = epsilon[atom->type[i]][atom->type[j]]*(ijsize10*r10inv + c0 + dipl + ddipl + dddipl);
    return evdwl;
  }
  else {
    fforce = 0;
    evdwl = 0;
    return evdwl;
  }
}

/* ---------------------------------------------------------------------- */

void *PairLjPoly::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut") == 0) return (void *) &cut;
  return NULL;
}
