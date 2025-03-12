#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

//x = gX/(2*t)            &     y = gX_0/(2*t)
//mu_tilde = mu/(2*t) - y &     t = 1
//19th February 2024: I am redefining lambda as: lambda = g^2/(K*t) at K = 1
//17th March 2024: mu = 2*mu_tilde - g^2/(Kt)*n = 2*mu_tilde - lambda*n

//I am using this version to plot 2*mu_tilde where 2*mu_tilde = mu_tilde ciuchi del KTV

double sgn(double x);
double energy(double k,double r,double delta,double mu_tilde,double x);
double omega(double delta,double n,double lambda,double mu_tilde,double x,int N_modes);
double derivative_omega_x(double delta,double lambda,double mu_tilde,double x,int N_modes);
double x_solution_local_minima(double delta,double lambda,double mu_tilde,int N_modes);
double particle_density(double delta,double mu_tilde,double x,int N_modes);
double mu_solution_given_n_and_x_local_minima(double delta,double n,double x,int N_modes,int N_mu_tilde);
int check_dimmer_or_supercoducting_phase(double delta,double n,double lambda,double x,double mu_tilde,int N_modes,int N_mu_tilde);

int main()
{
    //Fixed point
    double lambda = 4;
    //Various definitions of variables and its ranges of the variables
    int N_modes = 400;
    int N_delta = 2000;           double delta[N_delta+1],        delta1 = 0,     delta2 = 2;       //index i
    int N_mu_tilde = 2500;      double mu_tilde[N_mu_tilde+1],  mu_tilde1 = 0,  mu_tilde2 = 2.5;     //index j
    int N_x = 1200;             double x[N_x+1],                x1 = 0,         x2 = lambda/2;             //not used in the end...
    double n;                                                                                       //particle density
    
    //Assigning values to the variables
    int i,j,r;
    
    for(i=0;i<N_mu_tilde;i++){  mu_tilde[i] = mu_tilde1 + (mu_tilde2-mu_tilde1)*i/N_mu_tilde;}
    for(i=0;i<N_x;i++){                x[i] = x1 + (x2-x1)*i/N_x;}
    for(i=0;i<N_delta+1;i++){      delta[i] = delta1 + (delta2-delta1)*i/N_delta;}
    
    
    //************************ THE PROGRAM *******************************************
    
    //----dummy variables----
    double x_dummy,n_limit_dimmer;
    int check;
    int dimmer_exists; //0 if it doesn't exist, 1 if it does
    int null_derivative_topo_exists; //0 if it doesn't exist, 1 if it does
    
    double density_of_particles[2*N_mu_tilde+1], mu[2*N_mu_tilde+1];
    int j_star;
    //----important points----
    double n1,mu1, n2,mu2, n3,mu3, n4,mu4;
    double n_maxwell_topo_1, n_maxwell_topo_2, n_maxwell_dimmer1, n_maxwell_dimmer2;
    double mu_maxwell_topo, mu_maxwell_dimmer;
    double A;
    
    //-----file creation begin----
    FILE *maxw_dimhom1 = fopen("maxwell_dimmer_homogeneous_1.txt", "w");
    FILE *maxw_dimhom2 = fopen("maxwell_dimmer_homogeneous_2.txt", "w");
    FILE *maxw_topo1 = fopen("maxwell_nontopo_topo_1.txt", "w");
    FILE *maxw_topo2 = fopen("maxwell_nontopo_topo_2.txt", "w");
    FILE *spinodal_dimhom1 = fopen("spinodal_dimmer_homogeneous_1.txt", "w");
    FILE *spinodal_dimhom2 = fopen("spinodal_dimmer_homogeneous_2.txt", "w");
    FILE *spinodal_topo1 = fopen("spinodal_nontopo_topo_1.txt", "w");
    FILE *spinodal_topo2 = fopen("spinodal_nontopo_topo_2.txt", "w");
    
    char buffer[6],duffer[6];
    snprintf(buffer, sizeof(buffer), "%.3lf", lambda);
    char mu_vs_n[50] = "initial";
    char char_dummy[] = "mu_vs_n_l=";
    //----end file creation----
    
    for(i=1;i<N_delta+1;i++){
        
        printf("I am at delta = %.3lf out of %.3lf. progress = %.2lf percent\n",delta[i],delta[N_delta], delta[i]/delta[N_delta]*100.0);
        
        strcpy(mu_vs_n, char_dummy);
        snprintf(duffer, sizeof(duffer), "%.3lf", delta[i]);
        strcat(mu_vs_n, buffer);
        strcat(mu_vs_n, "_delta=");
        strcat(mu_vs_n, duffer);
        strcat(mu_vs_n, ".txt");
        FILE *plot_file = fopen(mu_vs_n, "w");
        
        dimmer_exists = 0;
        n_limit_dimmer = 0.5; //if there is no dimmer the limit will still be 0.5
        j_star = 0;
        null_derivative_topo_exists = 0;
        
        //----Writing in mu_vs_n file for the dimmerized phase
        for(j=0;j<N_mu_tilde;j++){
            
            check = 0;
            x_dummy = x_solution_local_minima(delta[i],lambda,mu_tilde[j],N_modes);
            if(x_dummy == 0){ j = N_mu_tilde;} //command to exit the loop since the system is not dimmerized anymore
            if( x_dummy != 0){
                
                n = particle_density(delta[i],mu_tilde[j],x_dummy,N_modes);
                check = check_dimmer_or_supercoducting_phase(delta[i],n,lambda,x_dummy,mu_tilde[j],N_modes,N_mu_tilde);
                if(check == 1){
                    dimmer_exists = 1;
                    n_limit_dimmer = n;
                    fprintf(plot_file, "%.12f %.12f %.12f %.12f \n", n,2*mu_tilde[j],delta[i],x_dummy);
                    //assign values to the arrays
                    density_of_particles[j] = n;
                    mu[j] = 2*mu_tilde[j]-lambda*n;
                    j_star = j;
                } //end if for check
            } //if x !=0 condition
        } //index j
        
        //----Writing in mu_vs_n file for the homogeneous phase
        for(j=0;j<N_mu_tilde;j++){
        
            n = particle_density(delta[i],mu_tilde[j],0,N_modes);
            if( n >= n_limit_dimmer){
                fprintf(plot_file, "%.12f %.12f %.12f %.12f \n", n,2*mu_tilde[j],delta[i],0.0);
                density_of_particles[j_star+1] = n;
                mu[j_star+1] = 2*mu_tilde[j]-lambda*n;
                j_star = j_star + 1;
            }
        } //index j
        
        fclose(plot_file); //up untill this point I have made a file 'mu_vs_n_l=x.xxx_d=x.xxx.txt'
        
       // printf("the n_limit_Dimmer is %.5lf\n",n_limit_dimmer);
        
        //Point 1
        if( dimmer_exists == 1 ){
            mu1 = -10000; //setting a starting impossible value
            for( r = 1; r < j_star; r++){
                if( density_of_particles[r] <= n_limit_dimmer){
                    if(mu[r] >= mu1){
                        mu1 = mu[r];
                        n1 = density_of_particles[r];
                    }
                }
            }
        }
        //Point 2
        if( dimmer_exists == 1){
            mu2 = 10000; //setting a starting impossible value
            for( r = 0; r < j_star; r++){
                if( density_of_particles[r] > n_limit_dimmer){
                    if( mu[r] < mu2 ){
                        mu2 = mu[r];
                        n2 = density_of_particles[r];
                    }
                    if( mu[r+1] > mu[r] ){
                        r = j_star + 1;
                    }
                }
            }
        }
        //Is it possible to make the Maxwell construction if there exists a null derivative around the topo - non topo phase
        for( r = 0; r < j_star-1; r++){
            if( density_of_particles[r] > n_limit_dimmer && mu[r+1] < mu[r]){
                null_derivative_topo_exists = 1;
                r = j_star;
                //printf("Pronaso sam jebenu null derivative za mu = %.5lf\n",mu[r]);
            }
        }
        //Point 4
        if( null_derivative_topo_exists == 1){
            for( r = j_star - 1; r > 1; r--){
                if( density_of_particles[r] > n_limit_dimmer){
                    if(mu[r-1] > mu[r] ){
                        mu4 = mu[r];
                        n4 = density_of_particles[r];
                        r = 0; //exit the loop
                    }
                }
            }
        }
        //Point 3
        if( null_derivative_topo_exists == 1){
            for( r = 0; r < j_star; r++){
                if( density_of_particles[r] > n_limit_dimmer && density_of_particles[r] < n4){
                    if(mu[r] > mu[r+1] && mu[r] > mu[r-1]){
                        mu3 = mu[r];
                        n3 = density_of_particles[r];
                        r = j_star + 1; //exit the loop
                    }
                }
            }
        }
        //----End Searching for the points 1,2,3,4
        
        if( dimmer_exists == 1){
            //printf("n1 = %.5lf mu1 = %.5lf\n", n1, mu1);
            fprintf(spinodal_dimhom1, "%.10f %.10f  \n", mu1+lambda*n1,delta[i]);
        }else{
            //printf("Not possible to find points 1 and 2. The dimmer does not exist.\n");
        }
        if( dimmer_exists == 1){
            //printf("n2 = %.5lf mu2 = %.5lf\n", n2, mu2);
            fprintf(spinodal_dimhom2, "%.10f %.10f  \n", mu2+lambda*n2,delta[i]);
        }
        if( null_derivative_topo_exists == 1){
            //printf("n3 = %.5lf mu3 = %.5lf\n", n3, mu3);
            fprintf(spinodal_topo1, "%.10f %.10f  \n", mu3+lambda*n3,delta[i]);
        }else{
            //printf("Not possible to find points 3 and 4. There is no Maxwell Construction.\n");
        }
        if( null_derivative_topo_exists == 1){
            //printf("n4 = %.5lf mu4 = %.5lf\n", n4, mu4);
            fprintf(spinodal_topo2, "%.10f %.10f  \n", mu4+lambda*n4,delta[i]);
        }
        
        //***********************************************
        //--------------------------------
        //----Maxwell Construction----
        
        //----Maxwell Construction for the topo-non topological phase
        if( null_derivative_topo_exists == 1 ){
            for( j = 0; j < j_star; j++){
                if( mu[j] < mu3 && mu[j] > mu4){
                    //find n_maxwell_topo_1
                    for( r = 0; r < j_star; r++){
                        if( density_of_particles[r] > n_limit_dimmer && mu[r] - mu[j] >= 0 ){
                            mu_maxwell_topo = mu[r];
                            n_maxwell_topo_1 = density_of_particles[r];
                            r = j_star;
                        }
                    }
                    //find n_maxwell_topo_2
                    for( r = j_star - 1; r > 1; r--){
                        if( mu[r] - mu[j] <= 0){
                            n_maxwell_topo_2 = density_of_particles[r];
                            r = 0;
                        }
                    }
                    //Calculating the Area bellow the curve
                    A = 0;
                    for( r = 0; r < j_star; r++){
                        if(density_of_particles[r] >= n_maxwell_topo_1 && density_of_particles[r] <= n_maxwell_topo_2){
                            A = A + (mu[r]-mu[j])*(density_of_particles[r]-density_of_particles[r-1]);
                        }
                    }
                    if( A <= 0 ){
                        //printf("n maxwell 1 is %.5lf and n maxwell 2 is %.5lf \n", n_maxwell_topo_1,n_maxwell_topo_2);
                        j = j_star; //exit loop
                        fprintf(maxw_topo1, "%.10f %.10f  \n", mu_maxwell_topo+lambda*n_maxwell_topo_1,delta[i]);
                        fprintf(maxw_topo2, "%.10f %.10f  \n", mu_maxwell_topo+lambda*n_maxwell_topo_2,delta[i]);
                        
                    }
                }
            }
        }
        //----Maxwell construction for the dimmerized-nontopo phase
        if( dimmer_exists == 1 ){
            for( j = 0; j < j_star; j++){
                if( mu[j] < mu1 && mu[j] > mu2){
                    //find n_maxwell_dimmer1
                    for( r = 0; r < j_star; r++){
                        if( mu[r] - mu[j] >= 0 ){
                            mu_maxwell_dimmer = mu[r];
                            n_maxwell_dimmer1 = density_of_particles[r];
                            r = j_star;
                        }
                    }
                    //find n_maxwell_dimmer2
                    for( r = 0; r < j_star -1; r++){
                        if( density_of_particles[r] > n_limit_dimmer && mu[r] - mu[j] >= 0){
                            n_maxwell_dimmer2 = density_of_particles[r];
                            r = j_star;
                        }
                    }
                    //Calculating the Area bellow the curve
                    A = 0;
                    for( r = 0; r < j_star; r++){
                        if(density_of_particles[r] >= n_maxwell_dimmer1 && density_of_particles[r] <= n_maxwell_dimmer2){
                            A = A + (mu[r]-mu[j])*(density_of_particles[r]-density_of_particles[r-1]);
                        }
                    }
                    if( A <= 0 ){
                        //printf("n maxwell 1 is %.5lf and n maxwell 2 is %.5lf \n", n_maxwell_topo_1,n_maxwell_topo_2);
                        j = j_star; //exit loop
                        fprintf(maxw_dimhom1, "%.10f %.10f  \n", mu_maxwell_dimmer+lambda*n_maxwell_dimmer1,delta[i]);
                        fprintf(maxw_dimhom2, "%.10f %.10f  \n", mu_maxwell_dimmer+lambda*n_maxwell_dimmer2,delta[i]);
                        
                    }
                }
            }
        }
        //----End Maxwell Construction----
        //------------------------------------
        //***********************************************
        
    } //index i
    
    fclose(maxw_dimhom1);
    fclose(maxw_dimhom2);
    fclose(maxw_topo1);
    fclose(maxw_topo2);
    fclose(spinodal_dimhom1);
    fclose(spinodal_dimhom2);
    fclose(spinodal_topo1);
    fclose(spinodal_topo2);
    

  return 0;
}

//************************ END PROGRAM *******************************************

//******************* DEFINITIONS OF FUNCTIONS *******************************************

double sgn(double x){
    if(x<0){x=-1;}
    if(x>0){x=1;}
    if(x==0){x=0;}
    return x;
}

double energy(double k,double r,double delta,double mu_tilde,double x){
    double epsilon = 0.00001; //solving a divergence
    double z = epsilon+sqrt( pow(cos(k),2) + pow(mu_tilde,2) + x*x + pow(delta*sin(k),2) + 2*r*sqrt( pow(cos(k)*mu_tilde,2) +x*x*(pow(delta*sin(k),2)+pow(mu_tilde,2)) ) );
    return z;
}

double derivative_omega_x(double delta,double lambda,double mu_tilde,double x,int N_modes){
    double derivative = 0;int i;double k;
    double epsilon = 0.01; //divergence solution for when mu_tilde = 0 notice the else that there are many places where I need to put the epsilon
    if( mu_tilde != 0){ //maybe I can put here an epsilon
        for(i=0;i<N_modes/2;i++){
            k = 2*M_PI*i/N_modes;
            derivative = derivative + ( 1 + (+1)*(pow(delta*sin(k),2)+pow(mu_tilde,2) )/sqrt(pow(cos(k)*mu_tilde,2)+x*x*(pow(delta*sin(k),2)+pow(mu_tilde,2))) )/energy(k,1,delta,mu_tilde,x);
            derivative = derivative + ( 1 + (-1)*( pow(delta*sin(k),2)+pow(mu_tilde,2) )/sqrt(pow(cos(k)*mu_tilde,2)+x*x*(pow(delta*sin(k),2)+pow(mu_tilde,2))) )/energy(k,-1,delta,mu_tilde,x);
        }
    }else{
        for(i=0;i<N_modes/2;i++){
            k = 2*M_PI*i/N_modes;
            derivative = derivative + ( 1 + (+1)*(pow(delta*sin(k),2)+pow(epsilon,2) )/sqrt(pow(cos(k)*(epsilon),2)+x*x*(pow(delta*sin(k),2)+pow(epsilon,2))) )/energy(k,1,delta,mu_tilde+epsilon,x);
            derivative = derivative + ( 1 + (-1)*( pow(delta*sin(k),2)+pow(epsilon,2) )/sqrt(pow(cos(k)*(epsilon),2)+x*x*(pow(delta*sin(k),2)+pow(epsilon,2))) )/energy(k,-1,delta,mu_tilde+epsilon,x);
        }
    }
    derivative = 4*x/lambda-x*derivative/N_modes;
    return derivative;
}

double x_solution_local_minima(double delta,double lambda,double mu_tilde,int N_modes){
    double x_local_minima = 0;int i;
    double x, x1 = lambda*0.4, x2 = 0; //I know analytically that 0=<x<=lambda/4
    double d_x = 0.005;
    int N_x_modified = fabs((x2-x1))/d_x;
    for(i=0;i<N_x_modified;i++){
        x = x1 - d_x*i;
        if(derivative_omega_x(delta,lambda,mu_tilde,x,N_modes)>0 && derivative_omega_x
           (delta,lambda,mu_tilde,x-d_x,N_modes)<0){
            x_local_minima = x-d_x/2;
            i = N_x_modified;
        }
    }
    return x_local_minima;
}

double particle_density(double delta,double mu_tilde,double x,int N_modes){
    double n = 0;int i;double k;
    if( mu_tilde != 0 && x!=0 ){
        for(i=0;i<N_modes/2;i++){
            k = 2*M_PI*i/N_modes;
            n = n + ( 1 + (+1)*( cos(k)*cos(k)+x*x )/sqrt(pow(cos(k)*(mu_tilde),2)+x*x*(pow(delta*sin(k),2)+pow(mu_tilde,2))) )/energy(k,1,delta,mu_tilde,x);
            n = n + ( 1 + (-1)*( cos(k)*cos(k)+x*x )/sqrt(pow(cos(k)*(mu_tilde),2)+x*x*(pow(delta*sin(k),2)+pow(mu_tilde,2))) )/energy(k,-1,delta,mu_tilde,x);
        }
        n = n*(mu_tilde)/(2*N_modes) + 0.5;
    }
    if( mu_tilde != 0 && x == 0 ){
        for(i=0;i<N_modes;i++){
            k = 2*M_PI*i/N_modes;
            n = n + (mu_tilde-cos(k))/sqrt(pow(cos(k)-mu_tilde,2)+pow(delta*sin(k),2));
        }
        n = 0.5 + n/(2*N_modes);
    }
    if( mu_tilde == 0 ){
        n = 0.5;
    }
    return n;
}

double mu_solution_given_n_and_x_local_minima(double delta,double n,double x,int N_modes,int N_mu_tilde){
    double mu_solution;int i;
    double mu, mu1 = 0, mu2 = 10; //this is actually mu_tilde but for some reason I wrote it as mu
    for(i=0;i<N_mu_tilde;i++){
        mu = mu1 + (mu2-mu1)*i/N_mu_tilde;
            if((particle_density(delta, mu,x,N_modes)-n<=0 && particle_density(delta,mu + (mu2-mu1)/N_mu_tilde ,x,N_modes)-n>=0) || particle_density(delta, mu,x,N_modes)-n==0  ){
                mu_solution = mu;
                i = N_mu_tilde;
            }
        }
    return mu_solution;
}

double omega(double delta,double n,double lambda,double mu_tilde,double x,int N_modes){
    double omega = 0;
    int i;
    double k;
    double y = -lambda*0.5*n;
    for(i=0;i<N_modes/2;i++){
        k = 2*M_PI*i/N_modes;
        omega = omega + energy(k,1,delta,mu_tilde,x);
        omega = omega + energy(k,-1,delta,mu_tilde,x);
    }
    omega = -omega/N_modes + 2*( x*x + y*y )/lambda - (mu_tilde);
    return omega;
}


int check_dimmer_or_supercoducting_phase(double delta,double n,double lambda,double x,double mu_tilde,int N_modes,int N_mu_tilde){
    
    double y = -lambda*0.5*n; //re-introducing this parameter
    double mu_tilde_sc = mu_solution_given_n_and_x_local_minima(delta,n,0,N_modes,N_mu_tilde);
    double mu_sc = 2*(mu_tilde_sc + y);
    double mu = 2*(mu_tilde + y);
    
    double F_sc = omega(delta,n,lambda,mu_tilde_sc,0,N_modes) + mu_sc*n;
    double F_dim = omega(delta,n,lambda,mu_tilde,x,N_modes) + mu*n;
    
    int result;
    
    if( F_dim < F_sc ){
        result = 1;// The system is in a dimmerizing phase
    }else{
        result = 0;// The system is in a superconducting phase
    }
    
    return result;
}
