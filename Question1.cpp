//Question 1
//It will take around 470 seconds
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define Nx 101 
#define Ny 101
#define pi 3.14159265
#define dx 2.0/(Nx-1)
#define dy 1.0/(Ny-1) 
int main(void){
	int i,j,L=2,W=1,n=0,m=0,count=0;
	double k1[Nx+1][Ny+1]={0},k2[Nx+1][Ny+1]={0},k3[Nx+1][Ny+1]={0},k4[Nx+1][Ny+1]={0},phi[Nx+1][Ny+1]={0},phi_new[Nx+1][Ny+1]={0};
	double alpha=1.0,dt=0.00001,dtt=0.1,err=1.0,t;
	double F1[Nx+1][Ny+1], F2[Nx+1][Ny+1], F3[Nx+1][Ny+1],F4[Nx+1][Ny+1];
	double lambda=0, sum1=0.0, sum2=0.0;
 	double theta[Nx+1][Ny+1]={0};
	FILE *fp, *fp1,*fp2;
	fp1=fopen("x_4_y_4.dat","w");
	fp2=fopen("analy_x_4_y_4.dat","w");
	
	 //Analytical Solution
	 count=0;t=0.0;	
    while(t<=4.1){
    	count=count+1;
    	t=count*dtt;
 	  for(i=1;i<=Nx;i++){
        for(j=1;j<=Ny;j++){sum1=0.0;sum2=0.0;
       		for(n=0;n<200;n++){
       			lambda=(2*n+1)*(pi/(2*L));
       			sum1=sum1+ (pow((-1),n)/(lambda*L))*exp(-alpha*pow(lambda,2)*t)*cos(lambda*(i-1)*dx);
               }
			for(m=0;m<200;m++){
       			lambda=(2*m+1)*(pi/(2*W));
       			sum2=sum2+ (pow((-1),m)/(lambda*W))*exp(-alpha*pow(lambda,2)*t)*cos(lambda*(j-1)*dy);
				}
		    theta[i][j]=4*sum1*sum2;
			}}
	        fprintf(fp2,"%lf\t%lf\n",t,theta[25][25]); 
            printf("Analytical:%lf\t%lf\n",t,theta[25][25]); 
            
			//at t=0.1 
            
		if(t==0.1){
		      fp=fopen("zero_one_analy.dat","w");		  
		      for(i=1;i<=Nx;i++){
              for(j=1;j<=Ny;j++){		
	          fprintf(fp,"%lf\t%lf\t%lf\n",(i-1)*dx,(j-1)*dy,theta[i][j]);		 
	        }}
	        fclose(fp);
			fp=fopen("analy_zero_one_mid_x.dat","w");  
              for(j=1;j<=Ny;j++){		
	        fprintf(fp,"%lf\t%lf\n",(j-1)*dy,theta[50][j]);		 
	        } 
	        fclose(fp);  
	    }
	    
	   	//at t=0.2
		   if(t==0.2){
		    fp=fopen("zero_two_analy.dat","w");		  
		    for(i=1;i<=Nx;i++){
            for(j=1;j<=Ny;j++){		
	        fprintf(fp,"%lf\t%lf\t%lf\n",(i-1)*dx,(j-1)*dy,theta[i][j]);		 
	        }}
	        fclose(fp);  
	    }
	    //at t=1.0
	    
	   	if(t==1.0){
		    fp=fopen("one_zero_analy.dat","w");		  
		    for(i=1;i<=Nx;i++){
            for(j=1;j<=Ny;j++){		
	        fprintf(fp,"%lf\t%lf\t%lf\n",(i-1)*dx,(j-1)*dy,theta[i][j]);		 
	        }}
	        fclose(fp);  
	    }
	    //at t=0.5
	    
		if(t==0.5){
	       fp=fopen("analy_zero_five_mid_y.dat","w");  
           for(i=1;i<=Nx;i++){		
	       fprintf(fp,"%lf\t%lf\n",(i-1)*dx,theta[i][50]);		 
	    }
	    fclose(fp);
	    } 
		}
	    fclose(fp2);

	//Initial condition and Boundary conditions
	t=0;
		for(i=1;i<=Nx;i++){
		    for(j=1;j<=Ny;j++){
	            if(i==Nx||j==Ny)
			      phi[i][j]=0.0;
			    else
			      phi[i][j]=1.0;
	            }}
	count=0;
	//Computational Solution
	while(err>1e-6){
	 err=0.0;
	 count=count+1;
	 t=dt*count;
       
        //F1
			for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        		F1[i][j]= alpha*(( (phi[i-1][j]-2*phi[i][j]+phi[i+1][j])/(dx*dx)) 
	                                               +((phi[i][j-1]-2*phi[i][j]+phi[i][j+1])/(dy*dy)));
			}
		} 
	    //K1
	    for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        		k1[i][j]=phi[i][j] +(dt/2)*F1[i][j];
			}
		} 
		//Boundary conditions on K1
		for(i=1;i<=Nx;i++){
		  for(j=1;j<=Ny;j++){
		    	if(i==Nx||j==Ny)
			       k1[i][j]=0.0;
			    else if(i==1)
			      k1[i][j]=k1[i+1][j];
			    else if(j==1)
			       k1[i][j]=k1[i][j+1];
				 }}
		//F2	
			for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){ 
				F2[i][j]= alpha*(( (k1[i-1][j]-2*k1[i][j]+k1[i+1][j])/(dx*dx)) 
	                                             +((k1[i][j-1]-2*k1[i][j]+k1[i][j+1])/(dy*dy))); 
	                    }}
		//K2      
        for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        		k2[i][j]=phi[i][j] +(dt/2)*F2[i][j];
			}
		}
		//Boundary conditions on K2
		for(i=1;i<=Nx;i++){
		  for(j=1;j<=Ny;j++){
		    	if(i==Nx||j==Ny)
			       k2[i][j]=0.0;
			    else if(i==1)
			      k2[i][j]=k2[i+1][j];
			    else if(j==1)
			       k2[i][j]=k2[i][j+1];
				 }} 
		 //F3      
        for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        		F3[i][j]= alpha*(( (k2[i-1][j]-2*k2[i][j]+k2[i+1][j])/(dx*dx)) 
	                                           +((k2[i][j-1]-2*k2[i][j]+k2[i][j+1])/(dy*dy)));
			}
		} 
		
		//K3      
        for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        		k3[i][j]=phi[i][j] +(dt)*F3[i][j];
			}
		} 
		//Boundary conditions on K3
		for(i=1;i<=Nx;i++){
		  for(j=1;j<=Ny;j++){
		    	if(i==Nx||j==Ny)
			       k3[i][j]=0.0;
			    else if(i==1)
			      k3[i][j]=k3[i+1][j];
			    else if(j==1)
			       k3[i][j]=k3[i][j+1];
				 }} 
		//F4
		for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        		F4[i][j]= alpha*(( (k3[i-1][j]-2*k3[i][j]+k3[i+1][j])/(dx*dx)) 
	                                           +((k3[i][j-1]-2*k3[i][j]+k3[i][j+1])/(dy*dy)));
			}
		} 
				 
		//phi(n+1)
		for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        	phi_new[i][j]=phi[i][j] +(dt/6)*(F1[i][j]+2*F2[i][j] +2*F3[i][j] +F4[i][j]); 			
			}}
	
        	
        //Boundary conditions on phi_new
		for(i=1;i<=Nx;i++){
		for(j=1;j<=Ny;j++){
			if(i==Nx||j==Ny)
			    phi_new[i][j]=0.0;
			else if(i==1)
			     phi_new[i][j]=phi_new[i+1][j];
			else if(j==1)
			     phi_new[i][j]=phi_new[i][j+1];
		 	      
		}}
		err=0.0;
		//Error calculation	
        for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
        	     err=err+fabs(phi_new[i][j]-phi[i][j]);//	printf("%lf	%lf\n",err,t);
        	}}
        	if(count%1000==1)
	     	printf("%lf\t%lf\n",err,t);		
	
	 for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){
              phi[i][j] = phi_new[i][j];
        	}}	
/////////////////////////////////////
   /////Printing Results////
////////////////////////////////////
   // at t=0.1 sec
    if(t==0.1){
      	fp=fopen("zero_one_code.dat","w");
	    for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){	
             fprintf(fp,"%lf	%lf	%lf\n",(dx*(i-1)), (dy*(j-1)),phi_new[i][j]);
        	}} 
	    fclose(fp);
	  
	    fp=fopen("zero_one_mid_x.dat","w");  
	    //for(i=1;i<=Nx;i++){
        for(j=1;j<=Ny;j++){		
	        fprintf(fp,"%lf\t%lf\n",(j-1)*dy,phi_new[50][j] );		 
	    }
		//}
	    fclose(fp);
	  }
	 ///////////////////////////// 
	/////////// at t=0.2 sec//////
	//////////////////////////////
    if(t==0.2){
      	fp=fopen("zero_two_code.dat","w");
	    for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){	
             fprintf(fp,"%lf	%lf	%lf\n",(dx*(i-1)), (dy*(j-1)),phi_new[i][j]);
        	}} 
	    fclose(fp);
	   
	  }
	////////////////////////////// 
    /////// at t=1.0 sec//////////
    //////////////////////////////
    if(t==1.0){
      	fp=fopen("one_zero_code.dat","w");
	    for(i=1;i<=Nx;i++){
        	for(j=1;j<=Ny;j++){	
             fprintf(fp,"%lf	%lf	%lf\n",(dx*(i-1)), (dy*(j-1)),phi_new[i][j]);
        	}} 
	    fclose(fp);
	 
	  }
    //////////////////////////////
	/////// at t=0.5 sec//////////
    //////////////////////////////
	if(t==0.5){
	 
	    fp=fopen("zero_five_mid_y.dat","w");  
        for(i=1;i<=Nx;i++){		
	        fprintf(fp,"%lf\t%lf\n",(i-1)*dx,phi_new[i][50]);		 
	    } 
	    fclose(fp);
	}
    //////////////////////////////
	/////at x=L/4 and y =W/4//////
	/////////////////////////////
 
			if(count%1000==1)
		  { // printf("L_4_W_4=%lf\t%lf\n",t,phi_new[25][25]);
	        fprintf(fp1,"%lf\t%lf\n",t,phi_new[25][25]);  
          }
          
     }
     printf("iterations are: %d\n",count);
     fclose(fp1);
} 
