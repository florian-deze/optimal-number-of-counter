# Projet de processus stochastiques
# Sujet : Optimalite du nombre guichets
# Binome : BRUNO Celine, DEZE Florian
# Les notations dans le code seront les suivantes:
# Mu : Parametre d'une loi exponentielle mesurant le temps de service d'un guichet
# Lambda : Parametre d'intensite d'un processus de Poisson mesurant l'arrivees des clients
# s : le nombre de guichet
# X : le nombre de clients dans le systeme au temps t
# PI0 : correspond a la mesure stationnaire de X au temps t=0
# LambdaN : parametre correspondant au temps de sejour du processus X au temps t>=0, a l'etat n

################################ QUESTION 1 ################################
	#Fonction

	# Le calcul de PI0 peut etre decoupe en deux sous-parties
	PI0Part1 = function(s,lambda,Mu){
  	res=0
  	for(k in 0:(s-1)){
    	res=res+(1/factorial(k))*((lambda/Mu)**k)
  	}
  	return (res)
	}

	PI0Part2 = function(s,lambda,Mu){
  	res=(1/factorial(s))*((lambda/Mu)**s)*((s*Mu)/(s*Mu-lambda))
  	return (res)
	}

	# Calcul de l'esperance de X au temps t
	Esperance= function(lambda,Mu,X,s){
  		(lambda/Mu)+X*(1/factorial(s))*((lambda/Mu)**s)*(lambda/(Mu*s))*(1/(1-(lambda/(Mu*s)))**2)
	}

	# fonction de mesure du nombre de client en fonction du temps dans le systeme avec s guichet
   	filedattente=function(Temp,n,lambda,Mu,s){
     nbC=c()
     time=c()
     nbclient=round(n)
     nbC=c(nbC,nbclient)
     time=c(time,0)
     x=lambda+min(nbclient,s)*Mu
     t=rexp(1,x)
     while(t<Temp){
       U=runif(1,0,1)
       if(U<(lambda/x)){nbclient=nbclient+1}
       else {nbclient=nbclient-1}
       x=lambda+min(nbclient,s)*Mu
       t=t+rexp(1,x)
       nbC=c(nbC,nbclient)
       time=c(time,t)
     }
     return(list("nbclient"=nbC, "temps"=time))
   	}

	#calcul

	# lambdaN=lambda+min(s*Mu,n*Mu)
	# PI0=1/(PI0Part1(s,lambda,Mu)+PI0Part2(s,lambda,Mu))
	# EsperancePI=Esperance(lambda,Mu,PI0,s)

	# 1.1 Etude de cas pertinents des parametre Mu et Lambda

    #CAS 1 : Mu*S > lambda
    lambda=4
    Mu=4
    s=4
    n=5
    lambdaN=lambda+min(s*Mu,n*Mu)
    PI0=1/(PI0Part1(s,lambda,Mu)+PI0Part2(s,lambda,Mu))
    EsperancePI=Esperance(lambda,Mu,PI0,s)
  
    #CAS 2 : s*Mu >> lambda
    lambda=4
    Mu=10
    s=10
    n=5
    lambdaN=lambda+min(s*Mu,n*Mu)
    PI0=1/(PI0Part1(s,lambda,Mu)+PI0Part2(s,lambda,Mu))
    EsperancePI=Esperance(lambda,Mu,PI0,s)
   
    #CAS 3: s*Mu > lambda mais tres proche
    lambda=7
    Mu=4
    s=2
    n=5
    lambdaN=lambda+min(s*Mu,n*Mu)
    PI0=1/(PI0Part1(s,lambda,Mu)+PI0Part2(s,lambda,Mu))
    EsperancePI=Esperance(lambda,Mu,PI0,s)
   
	# 1.2 Etude la stabilite du systeme en fonction du nombre de guichet s
	# Nous prendrons Lambda=5 et Mu=1 pour le reste des exemples de la question 1 et pour la question 2
	# Le nombre minimum de guichet est de 6
   
    lambda=5
    Mu=1
   
    # CAS 1 : 6 guichets
    s=6

    PI0=1/(PI0Part1(s,lambda,Mu)+PI0Part2(s,lambda,Mu))
    n=Esperance(lambda,Mu,PI0,s)
    nb_6=filedattente(500,100,lambda,Mu,s)
   
    plot("y"=nb_6$nbclient,"x"=nb_6$temps,type="l")
   
    # CAS 2 : 12 guichets
    s=12
   
    PI0=1/(PI0Part1(s,lambda,Mu)+PI0Part2(s,lambda,Mu))
    n=Esperance(lambda,Mu,PI0,s)
    nb_12=filedattente(500,100,lambda,Mu,s)
   
    plot("y"=nb_12$nbclient,"x"=nb_12$temps,type="l")
   
    # CAS 3 : 50 guichets
    s=50
   
    PI0=1/(PI0Part1(s,lambda,Mu)+PI0Part2(s,lambda,Mu))
    n=Esperance(lambda,Mu,PI0,s)
    nb_50=filedattente(500,100,lambda,Mu,s)
   
    plot("y"=nb_50$nbclient,"x"=nb_50$temps,type="l")


   
################################ QUESTION 2 ################################
	# FONCTION
	
	# fonction permettant le suivi d'une personne dans une file d'attente
	# Retourne le temps t passe par cette personne
	file_dattente_personne_unique= function(n,Mu,s){
          t = 0
          # temps avant d'arriver au guichet
          t = t + sum(rexp(max(0, n-s) , Mu))
          # temps passe au guichet
          t = t + rexp(1, Mu)
          
          return (t)
        }
    
    # fonction qui simule et calcul le nombre de personne presente dans la file d'attente
    PI_calcul_discret=function(lambda,mu,PI0,s){
        U=runif(1,0,1)
        indice=0
        p=0
        while(p<U){
          if(indice<s)
            p=p+(1/factorial(indice))*((lambda/mu)**indice)*PI0
          else
            p=p+(1/factorial(s))*((lambda/mu)**s)*((lambda/(Mu*s))**(indice-s))*PI0
          indice=indice+1 
        }
        return(indice)
      }

	# But : on cherche pour combien de guichet on a moins de 5% des personnes avec un temps d'attentes de 4 minutes
	#		on veut minimiser le nombre  guichet
	
	# Parametre de la simulation:
	# nb_simulation : nombre de simulation du temps d'attente d'une personne dans la file avec un nombre de guichet fixe
	# temps : correspond a l'incrementation du temps d'attente d'une personne dans une file pour un nombre de guichet fixe
	# temps_global : vecteur regroupant les temps moyens (temps/1000) d'attente en fonction du nombre de guichet s
	# nb_temps_inf_4 : variable comptant le nombre de personnes parmis les nb_simulation qui attendent plus de 4 mminutes
	# nb_temps_inf_4_global : vecteur contenant le pourcentage de personnes attendant plus de 4 minutes en fonction du nombre guichets s

    lambda=5
    mu=1
    nb_simulation=1000
    
    temps_global=c()
    nb_temps_inf_4_global=c()
    # On test avec un nombre de guichet allant de 6 a 20
    for (i in 6:20){
      s=i
      PI0=1/(PI0Part1(s,lambda,mu)+PI0Part2(s,lambda,mu))
      temps=0
      nb_temps_inf_4=0
      for(k in 1:nb_simulation){
        x=PI_calcul_discret(lambda,mu,PI0,s)
        tmp= file_dattente_personne_unique(x,mu,s)
        temps=temps+tmp
        if (tmp>4)
        	nb_temps_inf_4=nb_temps_inf_4+1
      }
      temps_global=c(temps_global,temps/nb_simulation)#moyenne des temps
      nb_temps_inf_4_global=c(nb_temps_inf_4_global,nb_temps_inf_4*100/nb_simulation)
    }
    
    #graphique montrant le temps d'attente moyen d'un client en fonction du nombre de guichet
    plot("y"=temps_global,"x"=6:20,type="l", main="temps d'attente moyen d'un client\n en fonction du nombre de guichet")
    
    #graphique montrant le pourcentage de personne attendant plus de 4 minutes
    #dans une file d'attente en fonction du nombre de guichets s
    plot("y"=nb_temps_inf_4_global,"x"=6:20,type="h", main="pourcentage de personnes attendant plus de 4 minutes\n en fonction du nombre de guichets")

	#le resultat attendu semble etre entre 8 et 10 guichets
    print(nb_temps_inf_4_global[3])
	print(nb_temps_inf_4_global[4])
	print(nb_temps_inf_4_global[5])