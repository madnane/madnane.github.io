//////////////////////////////// Question 3:////////////////////////
// Ici calcule la matrice H, n correspond a la dimension maximal des H, x est le vecteur des points ou seront calculés les polynomes



function Herm=Herm(n,x) //genere la  colonne (H_i(x)) pour i de 0 à K-1
    Herm(1)=1;
    Herm(2)=2*x;
    
    for i=2:n
        Herm(i+1)=2*x*Herm(i)-2*i*Herm(i-1);
    end
endfunction

function MatrixHerm=MatrixHerm(n,x) //génère la matrice (H_i(x_j)); 
    s=size(x)
    s=s(2)
    MatrixHerm=Herm(n,x(1));
    for i=2:s
      MatrixHerm=cat(2, MatrixHerm,Herm(n,x(i)));
    end
endfunction
//calcul de T

function T=T(K)
    N=0:1:K-1
    alpha=(sqrt(%pi)*2^N.*factorial(N)).^(-0.5)
    alpha=diag(alpha);
    fileName=strcat(['C:\Users\Yassine\Documents\Bose einstein\',"k",string(K),'.txt'])
    Casio=fscanfMat(fileName);
    omegaFoisExp=Casio(:,2).*exp((Casio(:,1).^2)/2);
    omegaFoisExp=diag(omegaFoisExp);
    T=alpha*MatrixHerm(K-1,Casio(:,1)')*omegaFoisExp;
endfunction

// calcul de Ti

function Ti=Ti(K)
    N=0:1:K-1
    alpha=(sqrt(%pi)*2^N.*factorial(N)).^(-0.5)
    alpha=diag(alpha);
    fileName=strcat(['C:\Users\Yassine\Documents\Bose einstein\',"k",string(K),'.txt'])
    Casio=fscanfMat(fileName);
    Exp=diag(exp(-Casio(:,1).^2/2));
    Ti=alpha*MatrixHerm(K-1,Casio(:,1)')*Exp;
    Ti=Ti';
endfunction


///////////////////////////////Question 4
// calcul de T de la question 4

function T4=T4(K)
    T4=kron(T(K),T(K))
endfunction

//calcul de Ti de la question 4:

function Ti4=Ti4(K)
    Ti4=kron(Ti(K),Ti(K))
endfunction



/////////////////////////////////////Question 5

// calcul de A qui dépend du pas delta_t, de K, de g et de mu


function A1=A1(deltat,K,mu)
    N=0:1:K-1
    A1=exp((-%i*(2*N+1-mu)*deltat)/2);
    A1=diag(A1);
endfunction
    

// calcul de la matrice \phi_2(x_p,t=0) pour p de 0 à K-1:
function phi2=phi2(deltat,K,mu)
   Cas=fscanfMat('C:\Users\Yassine\Documents\Bose einstein\k40.txt')
   C1=(%pi)^(-1/4)*exp((-Cas(:,1).^2)/2);
   phi2=Ti(K)*A1(deltat,K,mu)*T(K)*C1;
    
endfunction

//calcul de diag(exp(-i*deltat*g*|\phi_2|^2)):
function A2=A2(deltat,K,mu,g)
    
    A2=abs(phi2(deltat,K,mu)).^2;
    A2=diag(exp(-%i*deltat*g*A2));
    
endfunction

function B=B(deltat,K,mu)
    B=Ti(K)*A1(deltat,K,mu)*T(K);
    
endfunction

function A=A(deltat,K,mu,g)
      
       A=B(deltat,K,mu)*A2(deltat,K,mu,g)
       A=A*B(deltat,K,mu);  
endfunction


//Calcul de l'energie en un temps n,  var: n, de deltat, de K, de g et de mu

//cette fct renvoie phi(tm) (un vecteur)
function colonnePhi=colonnePhi(m,deltat,K,mu,g)
    Cas=fscanfMat('C:\Users\Yassine\Documents\Bose einstein\k40.txt');
    C2=(%pi)^(-1/4)*exp((-Cas(:,1).^2)/2);
    A=A(deltat,K,mu,g);
    colonnePhi=A^m*C2;
    
endfunction

// cette fct renvoie la matrice de taille K*(m+1) 
function colonnePhi2=colonnePhi2(m,deltat,K,mu,g)
    Cas=fscanfMat('C:\Users\Yassine\Documents\Bose einstein\k40.txt')
    C2=(%pi)^(-1/4)*exp((-Cas(:,1).^2)/2);
    n=m+1;
    colonnePhi2=ones(K,n);
    colonnePhi2(:,1)=C2;
    A=A(deltat,K,mu,g);
    for i=2:n
        colonnePhi2(:,i)=A*colonnePhi2(:,i-1);
    end
    
endfunction

function TiGen1D=TiGen1D(K,x)
    N=0:1:K-1
    alpha=(sqrt(%pi)*2^N.*factorial(N)).^(-0.5)
    alpha=diag(alpha);
    
    Exp=diag(exp(-x.^2/2));
    TiGen1D=alpha*MatrixHerm(K-1,x)*Exp;
    TiGen1D=TiGen1D';

endfunction
 // calcul Norme L2 Q5
 
function NormeL2Q5=NormeL2Q5(n,deltat,K,mu,g,nbpoint)
    
    X=T(K)*colonnePhi(n,deltat,K,mu,g);
    x=linspace(-4,4,nbpoint);
    X=TiGen1D(K)*T(K)*X;
    NormeL2Q5=norm((x(2)-x(1))*(X));
 endfunction

function Energie=Energie(n,deltat,K,mu,g,nbpoint)
     
    X=T(K)*colonnePhi(n,deltat,K,mu,g);
    x=linspace(-4,4,nbpoint);
    X=TiGen1D(K)*T(K)*X;
    Energie=0.5*sum((x(2)-x(1))*abs(X).^4)+sum(((0:K-1)*2+1).*abs((T(K)*colonnePhi(n,deltat,K,mu,g)))'.^2);
endfunction


///////////////////////////////////////////Question 6
    K=40;
    g=1;
    mu=10;
    deltat=0.01;

function phiQ6=phiQ6(n)
    // on prend g=1, mu=10,  les xi vont de -9 à +9 (segment contenants tous les noeuds pour K=40)
    // on prend la meme condition initiale qu'en q°5
    x=linspace(-8.1,8.1,K);
    sx=size(x);
    sx=sx(2);
    deltax=x(2)-x(1);
    phiQ6=ones(K,n);
    phiQ6(:,1)=(%pi)^(-1/4)*exp((-x.^2)/2)';
    ij=[[[1:K-1],[1:K],[2:K]]',[[2:K],[1:K],[1:K-1]]'];
    Rn=abs(phiQ6(:,1)).^2;
    one=ones(1,sx);
    for i=2:n
        Rn=2*abs(phiQ6(:,i-1)).^2-Rn;
        v=[repmat(-%i*deltat/(2*deltax^2),1,K-1),1+%i*deltat/2*(2/(deltax^2)+x.^2-mu+Rn'),repmat(-%i*deltat/(2*deltax^2),1,K-1)]';
        Mn=sparse(ij,v);
        //[hand,rk]=lufact(Mn);
        //disp(rk)
        phi=conjgrad(Mn,phiQ6(:,i-1));
        //ludel(hand);
        phiQ6(:,i)=2*phi-phiQ6(:,i-1);
    end
    
    
endfunction

function normeQ6=normeQ6(n)
      x=linspace(-8.1,8.1,K)
      N=0.0;
      phi=phiQ6(n);
      for i=1:K-1
          N=N+abs(phi(i))^2*(x(i+1)-x(i));
      end
      normeQ6=sqrt(N);
     
endfunction


    K=40;
    g=1;
    mu=10;
    deltat=0.01;

function phiQ6=phiQ6(n)
    // on prend g=1, mu=10,  les xi vont de -9 à +9 (segment contenants tous les noeuds pour K=40)
    // on prend la meme condition initiale qu'en q°5
    x=linspace(-8.1,8.1,K);
    sx=size(x);
    sx=sx(2);
    deltax=x(2)-x(1);
    phiQ6=ones(K,n);
    phiQ6(:,1)=(%pi)^(-1/4)*exp((-x.^2)/2)';
    ij=[[[1:K-1],[1:K],[2:K]]',[[2:K],[1:K],[1:K-1]]'];
    Rn=abs(phiQ6(:,1)).^2;
    one=ones(1,sx);
    for i=2:n
        Rn=2*abs(phiQ6(:,i-1)).^2-Rn;
        v=[repmat(-%i*deltat/(2*deltax^2),1,K-1),1+%i*deltat/2*(2/(deltax^2)+x.^2-mu+Rn'),repmat(-%i*deltat/(2*deltax^2),1,K-1)]';
        Mn=sparse(ij,v);
        //[hand,rk]=lufact(Mn);
        //disp(rk)
        phi=conjgrad(Mn,phiQ6(:,i-1));
        //ludel(hand);
        phiQ6(:,i)=2*phi-phiQ6(:,i-1);
    end
    
    
endfunction

function normeQ6=normeQ6(n)
      x=linspace(-8.1,8.1,K)
      N=0.0;
      phi=phiQ6(n);
      for i=1:K-1
          N=N+abs(phi(i))^2*(x(i+1)-x(i));
      end
      normeQ6=sqrt(N);
     
endfunction

///////////////////////////////////Question 9
//nous avons commencé par un schéma discrétisé en espace par la méthode des éléments finis, l'autre méthode est dans la partie ci apres "9 Bis"
K=30;
x=linspace(-8.1,8.1,K);
y=x;
gama=0.07;
deltax=x(2)-x(1);
deltat=0.01;
omega=1.60;
mu=10;
ex=0.01;
ey=0.01;// les epsilonnes sont ey et ex
phi_init=(%pi)^(-1/4)*exp((-x.^2)/2)'; 
phi_init=kron(phi_init, phi_init);// condition initiale: phi_init(x,y)=psi_0(x)*psi_0(y) avec psi_0 la premiere 

// construction de M_Delta
ij1=[[1:K^2],[1:K^2-1],[2:K^2],[1:K^2-K],[K+1:K^2]]';
ij2=[[1:K^2],[2:K^2],[1:K^2-1],[K+1:K^2],[1:K^2-K]]';
ij=[ij1,ij2];
vDelta=[ones(1,K-1),0];
vDelta=[repmat(-4,1,K^2),repmat(vDelta,1,K-1),ones(1,K-1),repmat(vDelta,1,K-1),ones(1,K-1),ones(1,2*(K^2-K))]';
MDelta=sparse(ij,vDelta);
MDelta=MDelta*1/deltax^2;


//Construciton de M_Lz
ij1=[[1:K^2-1],[2:K^2],[1:K^2-K],[K+1:K^2]]';
ij2=[[2:K^2],[1:K^2-1],[K+1:K^2],[1:K^2-K]]';
ij=[ij1,ij2];
tmp1=[ones(1,K-1),0];
tmp1=kron(x,tmp1);// produit de kroneker de x et tmp1
stmp1=size(tmp1);// taille de tmp
stmp1=stmp1(2);
tmp1=tmp1(1:(stmp1-1));
vL=[tmp1,tmp1,repmat(y,1,K-1),repmat(y,1,K-1)]';
ML=sparse(ij,vL);
ML=-%i/(2*deltax)*ML;


// Construction de M_x,y
  x=linspace(-8.1,8.1,K);
  y=linspace(-8.1,8.1,K);
  i=0:1:K^2-1
  Mxy=(1+ex)*(x(1+floor(i/K))).^2+(1+ey)*(y(1+i-K*floor(i/K))).^2
  Mxy=diag(Mxy)

function phi=phi(n) // retourne la mtrice de n colonne dont la i-ème colonne représente phi(i*delta)
    W=1-(%i-gama)/2*deltat/2*(-MDelta-omega*ML+Mxy-mu);// attention .. deltat/2 pour qu on puisse avoir phi1(tn+1/2) et non pas phi1 (tn+1)
    phi=ones(K^2,n);
    phi(:,1)=phi_init;
    for i=1:n-1
        phi1=phi(:,i);
        tmp=conjgrad(W,phi1,)//"bicg", 1d-8, K^3);
        phi1=2*tmp-phi1;
        phi3=phi1.*exp(deltat*abs(phi1)*(%i-gama)^-1);//=phi2(t_n+1))=phi3(tn) car  phi2(tn)=phi1(tn+1/2)=phi1
        tmp=conjgrad(W,phi3)//,"bicg", 1d-1, K^3);
        phi3=2*tmp-phi3;
        phi(:,i+1)=phi3./NormeL2(phi3);//phi4 de l'énoncé est egale à phi3./abs(phi3)
    end
    // enregistrer la partie reelle et imaginaire:
    fprintfMat('C:\Users\Yassine\Documents\Bose einstein\imaginary(phi).txt',imag(P));
    fprintfMat('C:\Users\Yassine\Documents\Bose einstein\real(phi).txt',real(P));
endfunction

function NormeL2=NormeL2(Phi)
    NormeL2=sqrt(sum(deltax*(abs(Phi).^2)));
endfunction

///////// La partie qui suit concerne la creation de graphique

// converti la matrice phi de dimention K^2*1 en une patrice Phi de dimension K*K ou les colonnes représentent les coordonnées x et les lignes les coordonnées y:
function ConvFromKron=ConvFromKron(Phi)
    s=size(Phi);
    s=s(1);
    s=sqrt(s);
    tmp=ones(s,s)
    for i=1:s
        tmp(i,:)=Phi(((i-1)*s+1):(i*s),1)';
    end
    ConvFromKron=tmp;
endfunction

////plotting
function plotting(n, Phi) //phi suposé contenir les m premieres valeurs de phi (en temps) avec m>=n
    M=Phi(:,n);
    M=ConvFromKron(M);
    
    gris = graycolormap(64)

    xset('colormap',gris(64:-1:1,:))
    grayplot(x,y,(abs(M))^2);
    
endfunction

function plottingSeq(Phi)
    s=size(Phi)
    s=s(1);
    s=sqrt(s);
    for i=1:s
        plotting(i,Phi);
        xpause(15*10^4);// ralenti l'execusion de n microsecondes
    end
endfunction



///////////////////////////////Question 9 bis
 
K=30;
Nx=50;
x=linspace(-4,4,Nx);
y=x;
gama=0.07;
deltax=x(2)-x(1);
deltat=0.01;
omega=1.78;
mu=10;
ex=0.7;
ey=0.7;// les epsilonnes sont ey et ex
//phi_init=(%pi)^(-1/4)*exp((-x.^2)/2)'; 
//phi_init=kron(phi_init, phi_init);// condition initiale: phi_init(x,y)=psi_0(x)*psi_0(y) avec psi_0 la premiere 
phi_init=zeros(1,K^2);
phi_init(1)=1;


N=0:1:K-1;
alpha=(sqrt(%pi)*2^N.*factorial(N)).^(-0.5);

// Matrice Mx2:

ijx2=[[[1:K^2],[(2*K+1):K^2],[1:(K^2-2*K)]];[[1:K^2],[1:(K^2-2*K)][(2*K+1):K^2]]]';
Vx2=[kron([0:K-1]+1/2,ones(1,K)),kron(0.25*alpha(1:K-2)./alpha(3:K),ones(1,K)),kron([2:K-1].*[1:K-2].*alpha(3:K)./alpha(1:K-2),ones(1,K))]';

Mx2=ex*sparse(ijx2, Vx2);

//Matrice My2:
ijy2=[[[1:K^2],[3:K^2],[1:K^2-2]];[[1:K^2],[1:K^2-2],[3:K^2]]]';
Vy2=[repmat(0.5+[0:K-1],1,K),0.25*alpha(1:K-2)./alpha(2:K-1),repmat([0,0,0.25*alpha(1:K-2)./alpha(3:K)],1,K-1),[2:K-1].*[1:K-2].*alpha(3:K)./alpha(1:K-2),repmat([0,0,[2:K-1].*[1:K-2].*alpha(3:K)./alpha(1:K-2)],1,K-1)]';

My2=ey*sparse(ijy2, Vy2);


//Matrice MLz:

ijLz=[[[1:K^2-K],[K+1:K^2]];[[K+1:K^2],[1:K^2-K]]]';
VLz=%i*[kron([alpha(2:K)./alpha(1:K-1)],[0,alpha(1:K-1)./alpha(2:K)]),-kron([alpha(1:K-1)./alpha(2:K)],[0,[1:K-1].*alpha(2:K)./alpha(1:K-1)])]';

MLz=sparse(ijLz, VLz);

// Mtrice MDelta de -Delta+x^2+y^2:

ijDelta= [[1:K^2];[1:K^2]]'
N=0:K^2-1;
VDelta=2*(floor(N/K)+N-K*floor(N/K)+1)';
MDelta=sparse(ijDelta, VDelta);

// La matrice M (on a pris delta /2 au lieu de delta t parce qu'on fait un demi pas de temps)
M=eye()-deltat/(4*(%i-gama))*(MDelta+Mx2+My2-mu*eye()-omega*MLz);
//Calcul des solutions:

// commençons d'abors pas construire une version "généralisée " de la matrice Ti définie en question 3 plus haut
function TiGen=TiGen(K)
    N=0:1:K-1;
    alpha=(sqrt(%pi)*2^N.*factorial(N)).^(-0.5);
    Diagalpha=diag(alpha);
    Exp=diag(exp(-x.^2/2));
    TiGen=Diagalpha*MatrixHerm(K-1,x)*Exp;

    TiGen=kron(TiGen,TiGen)';
endfunction
function phiBis=phiBis(n)
    //Transphibis: Transformée de fourier
    Ti=Ti4(K);
    TransphiBis=ones(K^2,n);
    TransphiBis(:,1)=phi_init';
    for i=1:n-1
        Transphi1=TransphiBis(:,i);
        tmp=conjgrad(M,Transphi1);
        Transphi1=2*tmp-Transphi1;
        Transphi2=Transphi1;
        //if (i==1) then// phi 2 constante en temps
            Tphi2=Ti*Transphi2;
        //end
        Transphi2=diag(1/(%i-gama)*deltat*abs(Tphi2))*Transphi2;
        Transphi3=conjgrad(M,Transphi2);
        TransphiBis(:,i+1)=Transphi3./norm(Transphi3);
    end
    phiBis=TiGen(K)*TransphiBis;
    // enregistrer la partie reelle et imaginaire:
    fprintfMat('C:\Users\Yassine\Documents\Bose einstein\imaginary(phiBis).txt',imag(phiBis));
    fprintfMat('C:\Users\Yassine\Documents\Bose einstein\real(phiBis).txt',real(phiBis));
endfunction

////plotting
function plottingBis(n, Phi) //phi suposé contenir les m premieres valeurs de phi (en temps) avec m>=n
    M=Phi(:,n);
    M=ConvFromKron(M);
    
    gris = graycolormap(128)

    xset('colormap',gris(128:-1:1,:))
    grayplot(x,y,(abs(M)).^2);
    
endfunction

function plottingSeqBis(Phi)
    s=size(Phi)
    s=s(2);
    //s=sqrt(s);
    for i=1:s
        plottingBis(i,Phi);
        xpause(10^4);// ralenti l'execusion de n microsecondes
    end
endfunction
