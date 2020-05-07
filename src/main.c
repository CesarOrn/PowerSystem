#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <cairo.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif



#define Ith(v,i)    NV_Ith_S(v,i-1)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   16               /* number of equations  */
#define Y1    RCONST(0.0)      /* initial y components */
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-8)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-14)
#define ATOL3 RCONST(1.0e-6)
#define T0    RCONST(0.0)      /* initia01 time           */
#define T1    RCONST(0.1)      /* first output time      */
#define TMULT RCONST(0.01)     /* output time factor     */
#define NOUT  10000           /* number of output times */

#define THRES 0.8
#define LINES 8
#define KOUPLE 170

#define ZERO  RCONST(0.0)

typedef struct {
  realtype dampingCo;
  realtype **couplingMatrix;
  realtype * powertype;
  realtype **power;
  realtype **dataOut;
  realtype threshold;
  int iteration;
}simData;

typedef struct {
  float r;
  float g;
  float b;
}color;

/* Private functions to output results */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int flin(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);


/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

void drawProducer(cairo_t *cr,double x, double y, double width, double height, color rgb );
void drawConsumer(cairo_t *cr,double x, double y, double radius, color rgb);
void drawLine(cairo_t *cr,double x, double y, double dx, double dy,color rgb);





int main(int argc, char const *argv[]) {
  realtype reltol, t, tout;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval, retvalr, iout;
  FILE *fp;



  realtype *coup[NEQ];
  for(int i = 0; i < NEQ; i++){
    coup[i] = (realtype*)malloc (NEQ*sizeof(realtype));
  }
  coup[0][4]=KOUPLE;
  coup[4][0]=KOUPLE;

  coup[2][4]=KOUPLE;
  coup[4][2]=KOUPLE;

  coup[6][4]=KOUPLE;
  coup[4][6]=KOUPLE;

  coup[4][10]=KOUPLE;
  coup[10][4]=KOUPLE;

  coup[10][6]=KOUPLE;
  coup[6][10]=KOUPLE;

  coup[8][10]=KOUPLE;
  coup[10][8]=KOUPLE;

  coup[12][10]=KOUPLE;
  coup[10][12]=KOUPLE;

  coup[12][14]=KOUPLE;
  coup[14][12]=KOUPLE;


  realtype *results[NOUT];
  for(int i=0;i<NOUT;i++){
    results[i] = (realtype*)malloc ((NEQ+1)*sizeof(realtype));
  }

  realtype *rowResults;

  realtype *powerLine[NOUT];
  for(int i=0;i<NOUT;i++){
    powerLine[i] = (realtype*)malloc ((LINES+1)*sizeof(realtype));
  }

  realtype powerArray[NOUT];
  powerArray[1]= RCONST(-10);
  powerArray[3]= RCONST(-10);
  powerArray[5]= RCONST(50);
  powerArray[7]= RCONST(-10);
  powerArray[9]= RCONST(-15);
  powerArray[11]= RCONST(20);
  powerArray[13]= RCONST(-15);
  powerArray[15]= RCONST(-10);

  simData Data;
  Data.dampingCo= RCONST(0.9);
  Data.couplingMatrix =coup;
  Data.powertype=powerArray;
  Data.power=powerLine;
  Data.threshold=RCONST(THRES);
  Data.dataOut =results;



  y = abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  y = N_VNew_Serial(NEQ);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(NEQ);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

  for(int i=1; i <NEQ+1;i++){
    Ith(y,i)  = ZERO;
  }

  reltol = RTOL;

  for(int i=1; i <NEQ+1;i++){
    Ith(abstol,i)  = ATOL3;
  }

  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

  retval = CVodeSetUserData(cvode_mem,&Data);
  if(check_retval(&retval, "FCVMALLOC", 0)) return(1);

  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  LS = SUNLinSol_Dense(y, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  printf(" \n3-species kinetics problem\n\n");


  iout = 0;  tout = T1;
  while(1) {
    Data.iteration= iout;
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    rowResults=*(results+iout);

    rowResults[0]=t;
    for (int i=1;i<NEQ+1;i++){
      rowResults[i]=Ith(y,i);
    //printf("data in %f \n",rowResults[i]);
    }

//store power of each line
    (*(Data.power+iout))[0]=  t;
    (*(Data.power+iout))[1]=  (Data.couplingMatrix[0][4])*sin(rowResults[5]-rowResults[1]);
    (*(Data.power+iout))[2]=  (Data.couplingMatrix[2][4])*sin(rowResults[5]-rowResults[3]);

    (*(Data.power+iout))[3]=  (Data.couplingMatrix[6][4])*sin(rowResults[5]-rowResults[7]);
    (*(Data.power+iout))[4]=  (Data.couplingMatrix[10][4])*sin(rowResults[5]-rowResults[11]);

    (*(Data.power+iout))[5]=  (Data.couplingMatrix[6][10])*sin(rowResults[11]-rowResults[7]);
    (*(Data.power+iout))[6]=  (Data.couplingMatrix[8][10])*sin(rowResults[11]-rowResults[9]);

    (*(Data.power+iout))[7]=  (Data.couplingMatrix[12][10])*sin(rowResults[11]-rowResults[13]);
    (*(Data.power+iout))[8]=  (Data.couplingMatrix[14][12])*sin(rowResults[13]-rowResults[15]);

    // check if line failes due to too much power


    if (Data.power[iout][1]>Data.couplingMatrix[2][0]*Data.threshold){
      //Data.couplingMatrix[2][0]=0;
    //  Data.couplingMatrix[0][2]=0;
    }

    if (Data.power[iout][2]>Data.couplingMatrix[4][0]*Data.threshold){
    //  Data.couplingMatrix[4][0]=0;
    //  Data.couplingMatrix[0][4]=0;
    }

    if (Data.power[iout][3]>Data.couplingMatrix[8][0]*Data.threshold){
    //  Data.couplingMatrix[8][0]=0;
    //  Data.couplingMatrix[0][8]=0;
    }
    if (Data.power[iout][4]>Data.couplingMatrix[4][6]*Data.threshold){
    //  Data.couplingMatrix[4][6]=0;
    //  Data.couplingMatrix[6][4]=0;
    }
    if (Data.power[iout][5]>Data.couplingMatrix[8][6]*Data.threshold){
    //  Data.couplingMatrix[8][6]=0;
    //  Data.couplingMatrix[6][8]=0;
    }
    if (Data.power[iout][6]>Data.couplingMatrix[10][8]*Data.threshold){
    //  Data.couplingMatrix[10][8]=0;
    //  Data.couplingMatrix[8][10]=0;
    }
    // time failure
    if (tout>30){
      Data.couplingMatrix[8][10]=0;
      Data.couplingMatrix[10][8]=0;
    }




    if (check_retval(&retval, "CVode", 1)) break;
    if (retval == CV_SUCCESS) {
      iout++;
      tout += TMULT;
    }

    //update coupling matrix for failure
    //update results
    //update power flow betwen lines. Becuase depends on coupling matrix.
    // update distance of failurs

    if (iout == NOUT) break;
    }
    printf("writing data");

    fp= fopen("/home/cesar/Code/C++/PowerSystem/src/DataOut.txt","w+");
    for (int i = 0; i< NOUT;i++){
      rowResults=*(results+i);
      fprintf(fp, "%2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f\n",rowResults[0], rowResults[1],rowResults[2],rowResults[3],rowResults[4],rowResults[5],rowResults[6],rowResults[7],rowResults[8],rowResults[9],rowResults[10],rowResults[11],rowResults[12]);
    }
    fclose(fp);

    fp= fopen("/home/cesar/Code/C++/PowerSystem/src/LinePower.txt","w+");
    for (int i = 0; i< NOUT;i++){
      rowResults=*(Data.power+i);
      fprintf(fp, "%2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f \n",rowResults[0],rowResults[1],rowResults[2],rowResults[3],rowResults[4],rowResults[5],rowResults[6],rowResults[7],rowResults[8]);

    }
    fclose(fp);

    /* Print some final statistics */
    PrintFinalStats(cvode_mem);

    /* Free y and abstol vectors */

    N_VDestroy(y);
    N_VDestroy(abstol);

    /* Free integrator memory */

    CVodeFree(&cvode_mem);

    /* Free the linear solver memory */

    SUNLinSolFree(LS);

    /* Free the matrix memory */

    SUNMatDestroy(A);


    //update power flow betwen lines
    //create
    printf("closeing file\n" );

    printf("start coupling matrix!!\n");
    for(int i = 0; i < NEQ; i++){
      free(*(coup+i));
    }
    printf("clean dataout!!\n");
    for(int i =0;i<NOUT;i++){
      free(*(results+i));
    }
    printf("clean power lines!!\n");
    for(int i=0;i<NOUT;i++){
      free(*(powerLine+i));
    }

    double xc = 128.0;
    double yc = 128.0;
    double radius = 100.0;
    double angle1 = 45.0  * (M_PI/180.0);  /* angles are specified */
    double angle2 = 180.0 * (M_PI/180.0);  /* in radians           */

    cairo_surface_t *surface;
    cairo_t *cr;

    color consumer;
    consumer.r=0.9;
    consumer.g=0.1;
    consumer.b=0.1;

    color producer;
    producer.r=0.15;
    producer.g=0.9;
    producer.b=0.15;

    color line;
    line.r=0.15;
    line.g=0.15;
    line.b=0.15;


    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 520, 520);
    cr = cairo_create (surface);
    drawLine(cr , 360,360,-160,-160, line);
    drawConsumer(cr, 200,200,20,consumer);
    drawProducer(cr,300,300,120,130,producer);


    //cairo_surface_write_to_png (surface, "/home/cesar/Code/C++/PowerSystem/src/image.png");


    return 0;

}


static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  simData *data=(simData*)user_data;
  //printf("helleo we got %5.2e and %d \n",data->dampingCo,data->iteration);

  realtype damp= data->dampingCo;
  realtype **coupling= data->couplingMatrix;
  realtype *p= data->powertype;

  for (int i =1; i <NEQ+1;i++ ){
    if(i%2==0){
      realtype rowsum =RCONST(0);
      for(int j= 1;j<NEQ+1;j++){
        rowsum= rowsum + coupling[i-2][j-1]*sin(Ith(y,j)-Ith(y,i-1));
      }
      Ith(ydot,i) = p[i-1]- damp*Ith(y,i) + rowsum;
    }else{
      Ith(ydot,i)=Ith(y,i+1);
    }
  }

  return(0);
}

static int flin(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  simData *data=(simData*)user_data;
  //printf("helleo we got %5.2e and %d \n",data->dampingCo,data->iteration);


  realtype damp= data->dampingCo;
  realtype **coupling= data->couplingMatrix;
  realtype *p= data->powertype;

  for (int i =1; i <NEQ+1;i++ ){
    if(i%2==0){
      realtype rowsum =RCONST(0);
      for(int j= 1;j<NEQ+1;j++){
        rowsum= rowsum + coupling[i-2][j-1]*(Ith(y,j)-Ith(y,i-1));
      }
      Ith(ydot,i) = p[i-1]- damp*Ith(y,i) + rowsum;
    }else{
      Ith(ydot,i)=Ith(y,i+1);
    }

  }

  return(0);
}

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4f      y =%3.6f  %3.6f  %3.6f\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);

  return;
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return 0;
}


void drawProducer(cairo_t *cr,double x, double y, double width, double height, color rgb ){
  cairo_set_source_rgb(cr, rgb.r,rgb.g, rgb.b);
  cairo_set_line_width(cr, 1);
  cairo_rectangle(cr, x, y, width, height);
  cairo_stroke_preserve(cr);
  cairo_fill(cr);
  cairo_close_path(cr);
}

void drawConsumer(cairo_t *cr,double x, double y, double radius, color rgb ){
  cairo_set_source_rgb(cr, rgb.r,rgb.g, rgb.b);
  cairo_set_line_width(cr, 1);
  cairo_arc(cr, x, y,radius, 0, 2*M_PI);
  cairo_stroke_preserve(cr);
  cairo_fill(cr);
  cairo_close_path(cr);
}

void drawLine(cairo_t *cr,double x, double y, double dx, double dy,color rgb){
  cairo_set_source_rgb(cr, rgb.r,rgb.g, rgb.b);
  cairo_set_line_width(cr, 5);
  cairo_move_to(cr,x,y);
  cairo_line_to(cr,x+dx,y+dy);
  cairo_stroke_preserve(cr);
  cairo_fill(cr);
  cairo_close_path(cr);

}
