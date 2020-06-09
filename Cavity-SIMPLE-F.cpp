
// Solving Ψ-Ω Formulation in a Lid-Driven Cavity with energy equatoion

      // necessary libraries

#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>

int main()
{
// ----------------------   Variables and discritization parameters  --------------------//
	int i,j,I,J;
    double L=1, H=1, l=0.1, h=0.1, Re=95, U0=1, residual=8.4e-6, error=1 ,b ,residualb;
    
    int xgrid_w=10, ygrid_w=10, Xgrid=60, Ygrid=60, it=0, xgridm, ygridm;

    double dxw , dyw , dxm , dym , x[61] , y[61] , X[62] , Y[62] , u[61][62],v[62][61], p[62][62] , us[61][62],vs[62][61], ps[62][62] , Pc_old[62][62] , Pc_new[62][62] , dx[60] , dy[60] , dX[61] , dY[61] , Fe , Fw, Fn, Fs , aP , aE , aW, aN , aS , lne , lnw , lnn , lns , lnPe , lnPw , lnPn , lnPs ;

    double d_iJ[61][62] , d_Ij[62][61], a_Pce[(60)*(60)], a_Pcw[(60)*(60)] , a_Pcn[(60)*(60)] , a_Pcs[(60)*(60)] , a_Pcp[(60)*(60)];


	//---------------------- Making structured unequal Mesh (Finer mesh near walls)-----------------------//
	
   xgridm = Xgrid-2*xgrid_w;
   ygridm = Ygrid-2*ygrid_w;
   dxw = l/xgrid_w;
   dyw = h/ygrid_w;   
   dxm = (L-2*l)/xgridm;
   dym = (H-2*h)/ygridm;   
   
   
   for (i=0;i<=xgrid_w;i++)
       x[i]=i*dxw;
   for (i=xgrid_w ; i<=(Xgrid-xgrid_w) ; i++)       
       x[i]=l+(i-xgrid_w)*dxm;
   for (i=Xgrid-xgrid_w ; i<=Xgrid ;i++)       
       x[i]=L - l + (i-Xgrid+xgrid_w)*dxw;       
  
  
   for (j=0;j<=ygrid_w;j++)
       y[j]=j*dyw;
   for (j=ygrid_w ; j<=(Ygrid-ygrid_w) ; j++)       
       y[j]=h+(j-ygrid_w)*dym;
   for (j=Ygrid-ygrid_w ; j<=Ygrid ;j++)       
       y[j]=H - h + (j-Ygrid+ygrid_w)*dyw; 
       
       //----------------------  X,Y,dx,dy  -----------------------//
       
   for (I=0;I<=xgrid_w;I++)
       X[I]=(I-0.5)*dxw;
   for (I=xgrid_w+1 ; I<=(Xgrid-xgrid_w) ; I++)       
       X[I]=l+(I-xgrid_w-0.5)*dxm;
   for (I=Xgrid-xgrid_w+1 ; I<=Xgrid+1 ;I++)       
       X[I]=L - l + (I-Xgrid+xgrid_w-0.5)*dxw; 

   for (J=0;J<=ygrid_w;J++)
       Y[J]=(J-0.5)*dyw;
   for (J=ygrid_w+1 ; J<=(Ygrid-ygrid_w) ; J++)       
       Y[J]=h+(J-ygrid_w-0.5)*dym;
   for (J=Ygrid-ygrid_w+1 ; J<=Ygrid+1 ;J++)       
       Y[J]=H - h + (J-Ygrid+ygrid_w-0.5)*dyw;

   for (i=0;i<Xgrid;i++)
       dx[i]=x[i+1]-x[i];
   for (I=0;I<Xgrid+1;I++)             
       dX[I]=X[I+1]-X[I];
   for (j=0;j<Ygrid;j++)
       dy[j]=y[j+1]-y[j];
   for (J=0;J<Ygrid+1;J++)             
       dY[J]=Y[J+1]-Y[J];
   //---------------------- initial guess for P,Ps,U,V  ------------------------//
   
   for (I=0;I<=Xgrid+1;I++)
   {
	   for (J=0;J<=Ygrid+1;J++)
	   {
           ps[I][J]=0;
           p[I][J]=0;
	   }
   }
   for (i=1;i<Xgrid;i++)
   {
       for (J=1;J<=Ygrid;J++)
	   {
           us[i][J]=U0*J/Ygrid;  
           u[i][J]=U0*J/Ygrid;                        
	   }
   }
   for (I=0;I<=Xgrid+1;I++)
   {
       for (j=0;j<=Ygrid;j++)
	   {
           vs[I][j]=0;                              
           v[I][j]=0;
	   }
   }
   //----------------------------------------------------------//   
   for (I=1;I<=Xgrid;I++)
   {
       for (J=1;J<=Ygrid;J++)
	   { 
           a_Pce[(I-1)*Xgrid+(J-1)]=0;
           a_Pcw[(I-1)*Xgrid+(J-1)]=0;
           a_Pcn[(I-1)*Xgrid+(J-1)]=0;                        
           a_Pcs[(I-1)*Xgrid+(J-1)]=0;
           a_Pcp[(I-1)*Xgrid+(J-1)]=0;
	   }
   }

   for (i=0;i<=Xgrid;i++)
       for (J=0;J<=Ygrid+1;J++)
           d_iJ[i][J]=0;
     
   for (I=0;I<Xgrid+2;I++)
       for (j=0;j<Ygrid+1;j++)
           d_Ij[I][j]=0;  

while ( error > residual)
{
   //-------------------------- boundary conditions for U,Us,V,Vs,P,Ps-------------------------//
	 it++;
	 printf("%d\t",it);
     for (J=0;J<Ygrid+2;J++)
	 {
         us[0][J]=0; 
         us[Xgrid][J] = 0;         
         u[0][J]=0; 
         u[Xgrid][J] = 0;           
     }
     for(I=0;I<Xgrid+2;I++)
	 {
         vs[I][0]=0;
         vs[I][Ygrid]=0;                      
         v[I][0]=0;
         v[I][Ygrid]=0;          
     }  
     for (j=0;j<Ygrid+1;j++)
	 {
         vs[0][j]=-vs[1][j];      
         vs[Xgrid+1][j] = -vs[Xgrid][j];           
         v[0][j]=-v[1][j];     
         v[Xgrid+1][j] = -v[Xgrid][j];
     }
     for(i=0;i<Xgrid+1;i++)
	 {
         us[i][0]=-us[i][1]; 
         us[i][Ygrid+1]=U0-us[i][Ygrid];                          
         u[i][0]=-u[i][1];
         u[i][Ygrid+1]=U0-u[i][Ygrid]; 
     }
	 ps[1][1]=0;
     p[1][1]=0;
     for (J=0;J<Ygrid+2;J++){    
         ps[0][J]=ps[1][J] - (us[1][J]-0.5*us[2][J])/(Re*dxw*dxw);
         ps[Xgrid+1][J] = ps[Xgrid][J] + (us[Xgrid][J]-0.5*us[Xgrid-1][J])/(Re*dxw*dxw);     
         p[0][J]=p[1][J]-(u[1][J] - 0.5*u[2][J])/(Re*dxw*dxw);
         p[Xgrid+1][J] = p[Xgrid][J] + (u[Xgrid][J]-0.5*u[Xgrid-1][J])/(Re*dxw*dxw);        
     }
     for(I=0;I<Xgrid+2;I++)
	 {      
         ps[I][0]=ps[I][1] - (vs[I][1]-0.5*vs[I][2])/(Re*dyw*dyw);   
         ps[I][Ygrid+1] = ps[I][Ygrid] + (vs[I][Ygrid]-0.5*vs[I][Ygrid-1])/(Re*dyw*dyw);          
         p[I][0]=p[I][1] - (v[I][1]-0.5*v[I][2])/(Re*dyw*dyw);      
         p[I][Ygrid+1] = p[I][Ygrid] + (v[I][Ygrid]-0.5*v[I][Ygrid-1])/(Re*dyw*dyw);              
     }                 
  //----------------------- u variables and u momnetom equation -----------------------//
     for (i=1;i<Xgrid;i++)
	 { 
         lne = ( X[i+1] - x[i] )/dx[i];          
         lnw = ( X[i] - x[i-1] )/dx[i-1] ;       
         lnn = ( x[i] - X[i] )/dX[i];           
         lns = lnn;           
         for (J=1;J<=Ygrid;J++)
		 {  
             lnPn = ( y[J] - Y[J] )/dY[J];         
             lnPs = ( y[J-1] - Y[J-1] )/dY[J-1]; 
             Fe = ( lne*us[i+1][J] + (1-lne)*us[i][J] )*dy[J-1];          
             Fw = ( lnw*us[i][J] + (1-lnw)*us[i-1][J] )*dy[J-1];
             Fn = ( lnn*vs[i+1][J] + (1-lnn)*vs[i][J] )*dX[i];
             Fs = ( lns*vs[i+1][J-1] + (1-lns)*vs[i][J-1] )*dX[i];  
             aE = -Fe*lne + dy[J-1]/(dx[i]*Re);
             aW = Fw*(1-lnw) + dy[J-1]/(dx[i-1]*Re);
             aN = -Fn*lnPn + dX[i]/(dY[J]*Re);
             aS = Fs*(1-lnPs) + dX[i]/(dY[J-1]*Re) ;    
             aP = aE + aW + aN +aS + Fe - Fw + Fn - Fs;        
             d_iJ[i][J]=dy[J-1]/aP;      
             u[i][J] = 0.5*( -(ps[i+1][J]-ps[i][J])*dy[J-1] + aE*us[i+1][J] + aW*us[i-1][J] + aN*us[i][J+1] + aS*us[i][J-1] )/aP + (1-0.5)*us[i][J];             
         }          
      }        
   //------------------------- v variables and v momnetom equation-----------------------------//
     for (j=1;j<Ygrid;j++)
	 {
         lne = (y[j]-Y[j])/dY[j];
         lnw = lne;
         lnn = (Y[j+1]-y[j])/dy[j];
         lns = (Y[j]-y[j-1])/dy[j-1];
         for (I=1 ; I<Xgrid+1 ; I++)
		 {
             lnPe = ( x[I]-X[I] )/dX[I]; 
             lnPw = ( x[I-1]-X[I-1] )/dX[I-1];             
             Fe = ( lne*us[I][j+1] + (1-lne)*us[I][j] )*dY[j];
             Fw = ( lnw*us[I-1][j+1] + (1-lnw)*us[I-1][j] )*dY[j];
             Fn = ( lnn*vs[I][j+1] + (1-lnn)*vs[I][j] )*dx[I-1];
             Fs = ( lns*vs[I][j] + (1-lns)*vs[I][j-1] )*dx[I-1]; 
             aE = -Fe*lnPe + dY[j]/(dX[I]*Re); 
             aW = Fw*(1-lnPw) + dY[j]/(dX[I-1]*Re) ;             
             aN = -Fn*lnn + dx[I-1]/(dy[j]*Re);             
             aS = Fs*(1-lns) + dx[I-1]/(dy[j-1]*Re) ;
             aP = aE + aW + aN +aS + Fe - Fw + Fn - Fs; 
             d_Ij[I][j]=dx[I-1]/aP;
             v[I][j]=0.5*( -(ps[I][j+1]-ps[I][j])*dx[I-1] + aE*vs[I+1][j] + aW*vs[I-1][j] + aN*vs[I][j+1] + aS*vs[I][j-1]  )/aP + (1-0.5)*vs[I][j];          
         }          
     }  
   //----------------------- pc, P correction equation ----------------------------//
   for (I=1 ; I<=Xgrid ; I++)           //initial//
       for (J=1 ; J<=Ygrid ; J++)
            Pc_old[I][J]=0;    

   for (I=1 ; I<=Xgrid ; I++)
   {
       for (J=1 ; J<=Ygrid ; J++)
	   { 
	    a_Pcs[(I-1)*(Ygrid)+(J-1)]= d_Ij[I][J-1]*dx[I-1];  
           a_Pcw[(I-1)*(Ygrid)+(J-1)]= d_iJ[I-1][J]*dy[J-1];  
           a_Pce[(I-1)*(Ygrid)+(J-1)]= d_iJ[I][J]*dy[J-1];  
           a_Pcn[(I-1)*(Ygrid)+(J-1)]= d_Ij[I][J]*dx[I-1]; 

           if (I==1)
              a_Pcw[(I-1)*(Ygrid)+(J-1)]=0 ;
           else if (I==Xgrid)
              a_Pce[(I-1)*(Ygrid)+(J-1)]=0;
           if (J==1)
              a_Pcs[(I-1)*(Ygrid)+(J-1)]=0; 
           else if(J==Ygrid)
              a_Pcn[(I-1)*(Ygrid)+(J-1)]=0;    
           a_Pcp[(I-1)*(Ygrid)+(J-1)]= a_Pcs[(I-1)*(Ygrid)+(J-1)]+ a_Pcw[(I-1)*(Ygrid)+(J-1)]+ a_Pce[(I-1)*(Ygrid)+(J-1)]+ a_Pcn[(I-1)*(Ygrid)+(J-1)];
       }
   } 
   //---------------------- b --------------------------//
       for (I=1 ; I<=Xgrid ; I++)
	   {
           for (J=1 ; J<=Ygrid ; J++)
		   {                                  
               b = (-u[I][J] + u[I-1][J] )*dy[J-1] + ( -v[I][J] + v[I][J-1] )*dx[I-1] ;              
               Pc_new[I][J] = (a_Pce[(I-1)*(Ygrid)+(J-1)]*Pc_old[I+1][J] + a_Pcw[(I-1)*(Ygrid)+(J-1)]*Pc_new[I-1][J] + a_Pcn[(I-1)*(Ygrid)+(J-1)]*Pc_old[I][J+1] + a_Pcs[(I-1)*(Ygrid)+(J-1)]*Pc_new[I][J-1] + b )/a_Pcp[(I-1)*(Ygrid)+(J-1)];                                                           
           }
       }         
    //-------------------------- update Pc ---------------------------// 
       for (I=1 ; I<=Xgrid ; I++)
           for (J=1 ; J<=Ygrid ; J++)        
               Pc_old[I][J]=0.6*Pc_old[I][J]+0.4*Pc_new[I][J]; 
  //------------------------ update Ps,us,vs ------------------------------//  
   for (I=1 ; I<=Xgrid ; I++)
       for (J=1 ; J<=Ygrid ; J++)
           ps[I][J] = ps[I][J] + 0.5*Pc_old[I][J];        

   for (i=1 ; i<Xgrid ; i++)        
       for (J=1 ; J<=Ygrid ; J++)
           us[i][J] = u[i][J] - d_iJ[i][J]*( Pc_old[i+1][J]-Pc_old[i][J] );
           
   for (I=1 ; I<=Xgrid ; I++)
       for (j=1 ; j<Ygrid ; j++) 
           vs[I][j]=v[I][j] - d_Ij[I][j]*( Pc_old[I][j+1]-Pc_old[I][j] );

   error=0;        
   for (I=1 ; I<=Xgrid ; I++)
       for (J=1 ; J<=Ygrid ; J++)
	   {
          residualb = abs((-us[I][J] + us[I-1][J])*dy[J-1] + (-vs[I][J] + vs[I][J-1])*dx[I-1]);
          if ( residualb > error)
             error=residualb;
       }

printf("residualb= %e\n", error);  
}

  //------------------------ plt file ------------------------------// 
FILE *output;

output = fopen("u.plt" , "w");
fprintf(output, "variables = x , Y , u \n zone i = %d j = %d\n", Xgrid+1, Ygrid+1);
    for(i = 0; i <= Xgrid; i++)
        for(J = 1; J <= Ygrid+1; J++)
              fprintf(output, "%e %e %e \n" , x[i] , Y[J] , us[i][J] );
    fclose(output);   
    
    output = fopen("v.plt" , "w");
    fprintf(output, "variables = X , y , v \n zone i = %d j = %d\n", Xgrid+1, Ygrid+1);
    for(I = 1; I <= Xgrid+1; I++)
        for(j = 0; j <= Ygrid; j++)
              fprintf(output, "%e %e %e \n" , X[I] , y[j] , vs[I][j] );
    fclose(output);  

    output = fopen("p.plt" , "w");
    fprintf(output, "variables = X , Y , p \n zone i = %d j = %d\n", Xgrid, Ygrid);
    for(I = 1; I <= Xgrid; I++)
        for(J = 1; J <= Ygrid; J++)
              fprintf(output, "%e %e %e \n" , X[I] , Y[J] , ps[I][J] );
    fclose(output);     
  

	printf("finished");
return 0;
}
