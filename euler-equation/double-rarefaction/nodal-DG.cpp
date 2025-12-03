#pragma GCC optimize(3,"Ofast","inline")
#include<bits/stdc++.h>

using namespace std;

int k,N;
double T_end,h,tau;

const double CFL=0.05;
const double PI=acos(-1.0);
const double gamma=7.0/5.0;
double now_time;

double LEFT,RIGHT;

class Matrix{

private:
    int dimension;
    std::vector<std::vector<double> >matrix;

public:

    Matrix(){}

    int getdim() const{return dimension;}

    void Init(int dim){
        dimension=dim;
        matrix.resize(dim);
        for(register int i=0;i<dim;++i)
        matrix[i].resize(dim,0.0);
    }

    double& operator()(int i,int j){return matrix[i][j];}
    const double& operator()(int i,int j)const{return matrix[i][j];}
};

Matrix result_tmp_matrix;

Matrix operator*(const double a,const Matrix I){

    for(register int i=0;i<I.getdim();i++)
    for(register int j=0;j<I.getdim();j++)
    result_tmp_matrix(i,j)=a*I(i,j);
        
    return result_tmp_matrix;
}

class Vector{

private:
    int dimension;
    std::vector<double>vector;

public:

    Vector(){}

    void Init(int dim){
        dimension=dim;
        vector.resize(dim,0.0);
    }

    int get_dim()const{return dimension;}

    double& operator()(int i){return vector[i];}
    const double& operator()(int i)const{return vector[i];}
};

Vector result_tmp_vector;

Vector operator*(const Matrix A,const Vector x){

    for(register int i=0;i<3;++i)
    result_tmp_vector(i)=0.0;

    for(register int i=0;i<3;++i)
    {
        for(register int j=0;j<3;++j)
        result_tmp_vector(i)+=(A(i,j)*x(j));
    }

    return result_tmp_vector;
}

Vector operator+(const Vector A,const Vector B){

    for(register int i=0;i<3;++i)
    result_tmp_vector(i)=A(i)+B(i);

    return result_tmp_vector;
}

Vector operator*(const double a,const Vector x){

    for(register int i=0;i<3;++i)
    result_tmp_vector(i)=a*x(i);

    return result_tmp_vector;
}

class Numerical_Integral{

    private:
        int Number_of_Points;
        vector<pair<double,double> >Gauss_Labatto_Points;
        Matrix Difference_Matrix;
        
    public:
        Numerical_Integral(){}

        int Get_Number_of_Points() const {return Number_of_Points;}

        void Init(int Degree_of_Poly){
            Number_of_Points=Degree_of_Poly+1;
            Difference_Matrix.Init(Degree_of_Poly+1);
            Gauss_Labatto_Points.resize(10,std::make_pair(0.0,0.0));
        }

        void Calculate_Gauss_Labatto_Points(){
            int Poin=Number_of_Points;
            if(Poin==2){
                Gauss_Labatto_Points[0]={-1.00000,1.0000000};
                Gauss_Labatto_Points[1]={1.000000000,1.0000000};
            }
            if(Poin==3){
                Gauss_Labatto_Points[0]={-1.00000,0.33333333333333333333};
                Gauss_Labatto_Points[1]={0.000000,1.33333333333333333333};
                Gauss_Labatto_Points[2]={1.0000000,0.33333333333333333};
            }
            if(Poin==4){
                Gauss_Labatto_Points[0]={-1.0000000,0.166666666666667};
                Gauss_Labatto_Points[1]={-0.447213595499958,0.8333333333333333};
                Gauss_Labatto_Points[2]={0.447213595499958,0.83333333333333333};
                Gauss_Labatto_Points[3]={1.0000000,0.1666666666666667};
            }
            return;
        }

        void Get_Difference_Matrix(){
            for(register int i=0;i<=k;++i)// i-th point
            for(register int j=0;j<=k;++j)// j-th 
            {
                double sum=0.0;
                if(i==j)
                {
                    for(register int l=0;l<=k;++l)
                    if(l==j)continue;
                    else sum+=(1.0/(Gauss_Labatto_Points[i].first-Gauss_Labatto_Points[l].first));
                    Difference_Matrix(i,j)=sum;
                }
                else
                {
                    sum=1.0;
                    double fac=1.0/(Gauss_Labatto_Points[j].first-Gauss_Labatto_Points[i].first);
                    for(register int l=0;l<=k;++l)
                    if(l==j||l==i)continue;
                    else sum*=(Gauss_Labatto_Points[i].first-Gauss_Labatto_Points[l].first)/(Gauss_Labatto_Points[j].first-Gauss_Labatto_Points[l].first);

                    Difference_Matrix(i,j)=fac*sum;
                }
            }
            return;
        }

        vector<pair<double,double> > get_GL(){return Gauss_Labatto_Points;}
        const vector<pair<double,double> >& get_GL()const{return Gauss_Labatto_Points;}
        Matrix get_DF(){return Difference_Matrix;}
        const Matrix get_DF()const{return Difference_Matrix;}
};

class State_of_Fluid{

    private:
        double density,velocity,pressure,momentum,energy;
        double space,time;
        Vector cal_data;
        Vector input_data;
    
    public:
        State_of_Fluid(){}

        void Init(double pho,double u,double p,double x,double t){
            density=pho,velocity=u,pressure=p,momentum=pho*u,energy=0.5*u*u*pho+p/(gamma-1.0),space=x,time=t,
            cal_data.Init(3);cal_data(0)=pho;cal_data(1)=momentum;cal_data(2)=energy;
            input_data.Init(3);input_data(0)=pho;input_data(1)=velocity;input_data(2)=pressure;
            return;
        }

        void Synchronize(){
            density=cal_data(0);momentum=cal_data(1);energy=cal_data(2);
            velocity=momentum/density;
            pressure=(gamma-1.0)*(energy-0.5*density*velocity*velocity);
            return;
        }

        const double& get_density()const{return density;}
        double& get_density(){return density;}

        const double& get_velocity()const{return velocity;}
        double& get_velocity(){return velocity;}
        
        const double& get_pressure()const{return pressure;}
        double& get_pressure(){return pressure;}

        const double& get_momentum()const{return momentum;}
        double& get_momentum(){return momentum;}

        const double& get_energy()const{return energy;}
        double& get_energy(){return energy;}

        Vector& get_cal_data(){return cal_data;}
        const Vector& get_cal_data()const {return cal_data;}
        Vector& get_input_data(){return input_data;}
};

class Nodal_Vector{
    private:
        int number_of_points;
        std::vector<State_of_Fluid>nodal_vector;
    public:
        Nodal_Vector(){};
        void Init(int num){
            number_of_points=num;
            nodal_vector.resize(num);
            return;
        }

        State_of_Fluid& operator()(int i){return nodal_vector[i];}
        const State_of_Fluid& operator()(int i)const{return nodal_vector[i];}
};

std::vector<Nodal_Vector>u_h,u_1,u_2,u_3;
Matrix Iden;
Numerical_Integral GL;
Vector Ini,Solution;

Vector Init_Function(const double x){
    // Ini(0)=1.0+(0.5*sin(x)*sin(x));Ini(1)=1.0;Ini(2)=2.0;

    // if(x<0.0){Ini(0)=1.0;Ini(1)=0.0;Ini(2)=1.0;}
    // else {Ini(0)=0.125,Ini(1)=0.0,Ini(2)=0.1;}

    // if(x<-4.0){Ini(0)=3.857143,Ini(1)=2.629369,Ini(2)=10.3333;}
    // else {Ini(0)=1.0+0.2*sin(5.0*x),Ini(1)=0.0,Ini(2)=1.0;}

    if(x<-0.0){Ini(0)=7.0,Ini(1)=-1.0,Ini(2)=0.2;}
    else {Ini(0)=7.0,Ini(1)=1.0,Ini(2)=0.2;}

    return Ini;
}

inline double Average_Under_Log(const double R,const double L){
    if(fabs(R-L)<1e-12)return (R+L)*0.5;
    return (R-L)/(log(R)-log(L));
}

Vector ZL,ZR,EC_Flux,Z_ave,Z_ave_log;

double betaL,betaR;

Vector Entropy_Conservative_Flux(const State_of_Fluid UL,const State_of_Fluid UR){

    betaL=UL.get_density()/UL.get_pressure()*0.5;
    betaR=UR.get_density()/UR.get_pressure()*0.5;

    EC_Flux(0)=Average_Under_Log(UR.get_density(),UL.get_density())*0.5*(UL.get_velocity()+UR.get_velocity());
    EC_Flux(1)=0.5*(UL.get_density()+UR.get_density())/(betaL+betaR)+0.5*(UL.get_velocity()+UR.get_velocity())*EC_Flux(0);
    EC_Flux(2)=(1.0/(2.0*(gamma-1.0)*Average_Under_Log(betaR,betaL))-0.125*(UL.get_velocity()+UR.get_velocity())*(UL.get_velocity()+UR.get_velocity()))*EC_Flux(0)+0.5*(UL.get_velocity()+UR.get_velocity())*EC_Flux(1);

    return EC_Flux;
}

Vector Phy_Tmp;

Vector Phy_Flux(const State_of_Fluid U){

    Phy_Tmp(0)=U.get_momentum();
    Phy_Tmp(1)=U.get_density()*U.get_velocity()*U.get_velocity()+U.get_pressure();
    Phy_Tmp(2)=(U.get_energy()+U.get_pressure())*U.get_velocity();

    return Phy_Tmp;
}

class Entropy_Stable_Flux{
    private:
        
    public:
        Entropy_Stable_Flux(){};

        double Sound_Speed(double dens,double pres){return sqrt(gamma*pres/dens);}

        double Q(double lef,double rig){
            if(lef<=rig)return 1.0;
            else return sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(lef/rig-1.0));
        }

        double pres_apporximation(const State_of_Fluid UL,const State_of_Fluid UR){
            double aL=Sound_Speed(UL.get_density(),UL.get_pressure());
            double aR=Sound_Speed(UR.get_density(),UR.get_pressure());
            double tmp1=(aL+aR+(gamma-1.0)*(UL.get_velocity()-UR.get_velocity())*0.5);
            double tmp2=aL/pow(UL.get_pressure(),(gamma-1.0)/(2.0*gamma))+aR/pow(UR.get_pressure(),(gamma-1.0)/(2.0*gamma));
            double tmp=pow(tmp1/tmp2,2.0*gamma/(gamma-1.0));
            return tmp;
        }

        double lambda(const State_of_Fluid UL,const State_of_Fluid UR){
            double pres_mid=pres_apporximation(UL,UR);
            double lambda_L=UL.get_velocity()-Sound_Speed(UL.get_density(),UL.get_pressure())*Q(pres_mid,UL.get_pressure());
            double lambda_R=UR.get_velocity()+Sound_Speed(UR.get_density(),UR.get_pressure())*Q(pres_mid,UR.get_pressure());
            double lambd=max(fabs(lambda_L),fabs(lambda_R));
            return lambd;
        }

        Vector Lax_Friedrichs_Flux(int left,int right,const std::vector<Nodal_Vector>& v){
            int lf=left,rig=right;

            if(left==0){lf=1;rig=1;}
            else if(right==N+1){lf=N;rig=N;}

            if(lf==rig&&lf==N){
                double lamb=lambda(v[lf](k),v[rig](k));
                lamb=max(fabs(v[lf](k).get_velocity())+Sound_Speed(v[lf](k).get_density(),v[lf](k).get_pressure()),lamb);
                lamb=max(fabs(v[rig](k).get_velocity())+Sound_Speed(v[rig](k).get_density(),v[rig](k).get_pressure()),lamb);;

                return -0.5*lamb*(v[rig](k).get_cal_data()+((-1.0)*v[lf](k).get_cal_data()))+0.5*(Phy_Flux(v[rig](k))+Phy_Flux(v[lf](k)));
            }

            if(lf==rig&&lf==N){
                double lamb=lambda(v[lf](0),v[rig](0));
                lamb=max(fabs(v[lf](0).get_velocity())+Sound_Speed(v[lf](0).get_density(),v[lf](0).get_pressure()),lamb);
                lamb=max(fabs(v[rig](0).get_velocity())+Sound_Speed(v[rig](0).get_density(),v[rig](0).get_pressure()),lamb);;

                return -0.5*lamb*(v[rig](0).get_cal_data()+((-1.0)*v[lf](0).get_cal_data()))+0.5*(Phy_Flux(v[rig](0))+Phy_Flux(v[lf](0)));
            }

            double lamb=lambda(v[lf](k),v[rig](0));
            lamb=max(fabs(v[lf](k).get_velocity())+Sound_Speed(v[lf](k).get_density(),v[lf](k).get_pressure()),lamb);
            lamb=max(fabs(v[rig](0).get_velocity())+Sound_Speed(v[rig](0).get_density(),v[rig](0).get_pressure()),lamb);;

            return -0.5*lamb*(v[rig](0).get_cal_data()+((-1.0)*v[lf](k).get_cal_data()))+0.5*(Phy_Flux(v[rig](0))+Phy_Flux(v[lf](k)));
        }
};

Entropy_Stable_Flux FL;

Vector RHS_Tmp,RHS_flux;

Vector RHS(int i,int j,const std::vector<Nodal_Vector>& v){

    for(int w=0;w<3;++w)
    RHS_Tmp(w)=0.0,RHS_flux(w)=0;

    for(int l=0;l<=k;++l)
    RHS_Tmp=RHS_Tmp+(-1.0)*( (GL.get_DF()(j,l)*Iden)*(2.0*Entropy_Conservative_Flux(v[i](j),v[i](l))));

    double b=0.0;
    if(j==0)b=-1.0;
    else if(j==k)b=1.0;
    b=b/GL.get_GL()[j].second;

    if(j==0)RHS_flux=FL.Lax_Friedrichs_Flux(i-1,i,v);
    else if(j==k)RHS_flux=FL.Lax_Friedrichs_Flux(i,i+1,v);      

    RHS_flux=(RHS_flux+(-1.0)*Phy_Flux(v[i](j)));
    RHS_flux=b*RHS_flux;

    RHS_Tmp=RHS_Tmp+(-1.0)*RHS_flux;

    return (2.0/h)*RHS_Tmp;
}

inline void SSP(){
    
    now_time=0.0;
    for(;now_time<T_end;){
        if(now_time+tau<=T_end)now_time+=tau;
        else {tau=T_end-now_time;now_time=T_end;}

        // cout<<now_time<<endl;

        for(int i=1;i<=N;++i)for(int j=0;j<=k;++j)u_1[i](j).get_cal_data()=u_h[i](j).get_cal_data()+tau*RHS(i,j,u_h),u_1[i](j).Synchronize();

        for(int i=1;i<=N;++i)for(int j=0;j<=k;++j)u_2[i](j).get_cal_data()=0.75*u_h[i](j).get_cal_data()+0.25*(u_1[i](j).get_cal_data()+tau*RHS(i,j,u_1)),u_2[i](j).Synchronize();

        for(int i=1;i<=N;++i)for(int j=0;j<=k;++j)u_3[i](j).get_cal_data()=0.333333333333333333333333*u_h[i](j).get_cal_data()+0.666666666666666667*(u_2[i](j).get_cal_data()+tau*RHS(i,j,u_2)),u_3[i](j).Synchronize();

        for(int i=1;i<=N;++i)for(int j=0;j<=k;++j)u_h[i](j)=u_3[i](j);

        double maxx=0.0;
        for(register int i=1;i<=N;++i)
        {
            for(register int j=0;j<=k;++j)
            {
                double sound=sqrt(gamma*u_h[i](j).get_pressure()/u_h[i](j).get_density());
                double lamb=u_h[i](j).get_velocity();
                maxx=max(maxx,max(fabs(lamb),max(fabs(lamb+sound),fabs(lamb-sound))));
            }
        }
        maxx=max(1.5,maxx);

        tau=CFL*h/maxx;
    }

    return;
}

inline void init(){

    h=(RIGHT-LEFT)/N;

    GL.Init(k);
    GL.Calculate_Gauss_Labatto_Points();
    GL.Get_Difference_Matrix();

    u_h.resize(N+5);u_1.resize(N+5);u_2.resize(N+5);u_3.resize(N+5);

    for(int i=0;i<=N+1;++i){
        u_h[i].Init(k+1);
        u_1[i].Init(k+1);
        u_2[i].Init(k+1);
        u_3[i].Init(k+1);
    }

    for(int i=1;i<=N;++i)
    {
        double left=LEFT+(i-1)*h;
        double right=LEFT+i*h;

        for(register int j=0;j<=k;++j)
        {
            u_h[i](j).get_input_data()=Init_Function(GL.get_GL()[j].first*h*0.5+(left+right)*0.5);
            u_h[i](j).Init(u_h[i](j).get_input_data()(0),u_h[i](j).get_input_data()(1),u_h[i](j).get_input_data()(2),GL.get_GL()[j].first*h*0.5+(left+right)*0.5,0.0);
        }
    }

    double maxx=0.0;
    for(register int i=1;i<=N;++i)
    {
        for(register int j=0;j<=k;++j)
        {
            double sound=sqrt(gamma*u_h[i](j).get_pressure()/u_h[i](j).get_density());
            double lamb=u_h[i](j).get_velocity();
            maxx=max(maxx,max(fabs(lamb),max(fabs(lamb+sound),fabs(lamb-sound))));
        }
    }
    maxx=max(1.5,maxx);

    tau=CFL*h/maxx;

    return;
}

void Final_Check(){

    for(register int i=1;i<=N;++i)
    cout<<LEFT-0.5*h+h*i<<',';cout<<endl;

    for(register int i=1;i<=N;++i)
    {
        double tmp=0.0;
        for(register int j=0;j<=k;++j)
        tmp+=u_h[i](j).get_density()*GL.get_GL()[j].second;
        
        cout<<tmp*0.5<<",";
    }
    cout<<endl;

    for(register int i=1;i<=N;++i)
    {
        double tmp=0.0;
        for(register int j=0;j<=k;++j)
        tmp+=u_h[i](j).get_momentum()*GL.get_GL()[j].second;
        
        cout<<tmp/2<<",";
    }
    cout<<endl;

    for(register int i=1;i<=N;++i)
    {
        double tmp=0.0;
        for(register int j=0;j<=k;++j)
        tmp+=u_h[i](j).get_energy()*GL.get_GL()[j].second;
        
        cout<<tmp*0.5<<",";
    }
    cout<<endl;

    return;
}

int main(){

    Iden.Init(3);for(register int i=0;i<3;++i)Iden(i,i)=1.0;
    result_tmp_vector.Init(3);
    result_tmp_matrix.Init(3);
    ZL.Init(3);ZR.Init(3);EC_Flux.Init(3);Z_ave.Init(3);Z_ave_log.Init(3);
    RHS_Tmp.Init(3);RHS_flux.Init(3);
    Phy_Tmp.Init(3);
    Ini.Init(3);Solution.Init(3);

    for(int Test=1;Test<=5;++Test)
    {
        cout<<"Test "<<Test<<": "<<endl;
        cout<<"Please Input the Degree of Polynomials:";cin>>k;
        cout<<"Please Input the Time of End:";cin>>T_end;
        cout<<"Please Input the Number of Girds:";cin>>N;
        cout<<"Please Input the Computational Domain:";cin>>LEFT>>RIGHT;
            
        init();

        SSP();

        Final_Check();
    }

    return 0;
}