#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

main()
{    FILE *fp1;
     int lattice_size_x,*iup,*idown,i;
     float deltar,nu,ro,rn,*height,*heightold,*temp,*tempold,*flux;
     float simulationlength,outputinterval,timestep,elapsedtime,maxdelt,maxdelh,del;
		 
     fp1=fopen("viscousflownu100","w");
     lattice_size_x=260;
     deltar=0.02;
     nu=100.0;
     ro=1.0;
     rn=1.1;
     iup=ivector(1,lattice_size_x);
     idown=ivector(1,lattice_size_x);
     height=vector(1,lattice_size_x);
     heightold=vector(1,lattice_size_x);
     temp=vector(1,lattice_size_x);
     tempold=vector(1,lattice_size_x);
     flux=vector(1,lattice_size_x);
     simulationlength=10.0; 
     outputinterval=0.1;
     for (i=1;i<=lattice_size_x;i++)
      {temp[i]=exp(-pow(i*deltar/rn,20));
       tempold[i]=temp[i];
       if (i*deltar<=ro) height[i]=pow(4.0/(1+nu)*log(rn/ro),0.25);
        else if (i*deltar<=rn) height[i]=
         pow(4.0/(1+nu)*log(rn/(i*deltar)),0.25); else height[i]=0.001;
       heightold[i]=height[i];
       flux[i]=0;
       iup[i]=i+1;idown[i]=i-1;}
     iup[lattice_size_x]=lattice_size_x;idown[1]=1;
     timestep=0.00001;
     elapsedtime=0.0;
     while (elapsedtime<simulationlength)
      {maxdelt=0;maxdelh=0;
       flux[1]=1;
       for (i=1;i<=lattice_size_x;i++)
        if (i*deltar<ro) flux[i]=1; 
         else flux[i]=i*(1+nu*tempold[i])*pow(0.5*(heightold[idown[i]]+heightold[i]),3)*
          (heightold[idown[i]]-heightold[i]);
       for (i=1;i<=lattice_size_x-1;i++) 
        {del=timestep/(i*deltar*deltar)*(flux[i]-flux[iup[i]]);
         height[i]+=del;
         if (fabs(del)>maxdelh) maxdelh=del;
         del=0;
         if (heightold[i]>0) 
          {del=timestep*(flux[i]*(tempold[idown[i]]-tempold[i])/
            (heightold[i]*i*deltar*deltar)-
              tempold[i]/(heightold[i]*heightold[i]));
	     if (i*deltar<ro) del=0;
         temp[i]+=del;} else del=0;
         if (fabs(del)>maxdelt) maxdelt=fabs(del);
         if (temp[i]<0) temp[i]=0;
         if (temp[i]>1) temp[i]=1;}
       for (i=1;i<=lattice_size_x;i++)
        if (i*deltar<ro) height[i]=height[(int)(ro/deltar)]; 
       elapsedtime+=timestep;
       if ((maxdelt>0.1)||(maxdelh>0.1))
         {elapsedtime-=timestep;
          timestep/=2.0;
          for (i=1;i<=lattice_size_x;i++)
           {height[i]=heightold[i];
            temp[i]=tempold[i];}}
        else 
         if ((maxdelt<0.01)&&(maxdelh<0.01)) timestep*=1.2;
       for (i=1;i<=lattice_size_x;i++)
         {heightold[i]=height[i];
	  tempold[i]=temp[i];}}
       for (i=1;i<=lattice_size_x;i++)
        fprintf(fp1,"%d %f %f\n",i,height[i],temp[i]);
       fclose(fp1);
}  
