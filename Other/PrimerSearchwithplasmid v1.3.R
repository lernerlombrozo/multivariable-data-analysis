library("Biostrings")
source("http://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
setwd("/Users/MacBook/Fasta")

#Funci??n para calcular la TM
Tmf<- function(SequenceA) {
  A=letterFrequency(SequenceA, "A", as.prob=FALSE)
  Ti=letterFrequency(SequenceA, "T", as.prob=FALSE)
  C=letterFrequency(SequenceA, "C", as.prob=FALSE)
  G=letterFrequency(SequenceA, "G", as.prob=FALSE)
  if (length(SequenceA)>13){
    64.9+41*(G+C-16.4)/(G+C+A+Ti)
  }else{
    2*(A+Ti)+4*(G+C) 
  }
}

Tmf2<- function(SequenceA){
  A=letterFrequency(SequenceA, "A", as.prob=FALSE)
  Ti=letterFrequency(SequenceA, "T", as.prob=FALSE)
  C=letterFrequency(SequenceA, "C", as.prob=FALSE)
  G=letterFrequency(SequenceA, "G", as.prob=FALSE)
  81.5+41*((G+C)/(G+C+Ti+A))-675/(G+C+A+Ti) 
}

#Funci??n para buscar los Hairpins:
Hairpins = function(primer){
  DNAprimer = DNAString(primer)
  HPpresence=FALSE
  j=6
  i=j
  ciclo=1
  while(j<=ceiling(length(DNAprimer)-4)/2 & HPpresence==FALSE){
    while(i<=(length(DNAprimer)-(j+4)) & HPpresence==FALSE){
      A=(DNAprimer[ciclo:i])
      B=reverseComplement(DNAprimer[(i+4+1):length(DNAprimer)])
      HairPin <- countPattern(A,B, max.mismatch=(j-6),
                              min.mismatch=(j-6), fixed=TRUE)
      HPpresence= HairPin>0
      i=i+1
      ciclo=ciclo+1
    }
    j=j+1
    i=j
    ciclo=1
  }
  return(HPpresence)
}


#Funci??n para ver si se forman d??meros
Dimers =function(PrimerA,PrimerB){
  DNAPrimerA = DNAString(PrimerA)
  InversePrimerA=DNAPrimerA[length(DNAPrimerA):1]
  RCPrimerA=reverseComplement(DNAPrimerA)
  DNAPrimerB = DNAString(PrimerB)
  RCPrimerB=reverseComplement(DNAPrimerB)
  InversePrimerB=DNAPrimerB[length(DNAPrimerB):1]
  Dimer=matchPattern(DNAPrimerA[(length(DNAPrimerA)-5):length(DNAPrimerA)],
                     RCPrimerA,max.mismatch=1, min.mismatch=0, fixed=TRUE)
  if(is.na(Dimer@ranges@width[1])){
    print("FW,FW")
    print(Dimer@ranges@width[1])
    Dimer=matchPattern(DNAPrimerB[(length(DNAPrimerB)-5):length(DNAPrimerB)],
                       RCPrimerB,max.mismatch=1, min.mismatch=0, fixed=TRUE)
    if(is.na(Dimer@ranges@width[1])){
      print("BW,BW")
      print(Dimer@ranges@width[1])
      Dimer=matchPattern(DNAPrimerA[(length(DNAPrimerA)-5):length(DNAPrimerA)],
                       RCPrimerB[length(RCPrimerB):1],max.mismatch=1, min.mismatch=0, fixed=TRUE)
      if(is.na(Dimer@ranges@width[1])){
        print("FW,BW")
        print(Dimer@ranges@width[1])
        Dimer=matchPattern(DNAPrimerB[(length(DNAPrimerB)-5):length(DNAPrimerB)],
                       RCPrimerA[length(RCPrimerA):1],max.mismatch=1, min.mismatch=0, fixed=TRUE)
        if(is.na(Dimer@ranges@width[1])){
          print("BW,FW")
          print(Dimer@ranges@width[1])
          return(FALSE)
        }else{
          return(TRUE)
        }
      }else{
        return(TRUE)
      }
    }else{
      return(TRUE)
    }
  }else{
    return(TRUE)
  }
}

# Funci??n para Gener matriz de primers (parte del pl??smido)
PosiblesPrimersEnPlasmido<- function(LengthPrimerEnPlasmido,Enzima,Plasmido){
  DNAEnzima=DNAString(Enzima)
  matchP = matchPattern(Enzima,
           Plasmido[[1]], max.mismatch=0,min.mismatch=0, fixed=TRUE)
  ParteDePlasmidoEnPrimerMatriz=matrix(, nrow = 0, ncol = 3)
    for(i in 0:(LengthPrimerEnPlasmido-length(DNAEnzima))){
    CB=Plasmido[[1]][(matchP@ranges@start-(LengthPrimerEnPlasmido-length(DNAEnzima)-i)):(matchP@ranges@start+matchP@ranges@width+i-1)]
    CC=toString(CB)
    CD=c(CC,(matchP@ranges@start-(LengthPrimerEnPlasmido-length(DNAEnzima)-i)),Tmf2(CB))
    ParteDePlasmidoEnPrimerMatriz=rbind(CD,ParteDePlasmidoEnPrimerMatriz)
  }
  colnames(ParteDePlasmidoEnPrimerMatriz)<- c("Primer","Location","Tm")
  return(ParteDePlasmidoEnPrimerMatriz)
}

# Funci??n para Gener matriz de primers (parte del Gen)
PosiblesPrimersEnGen<- function(Repeats,length,DNAbuscar,DNAencontrar,CHR,FW){
  i=1
  PosiblesPrimers=matrix(, nrow = 0, ncol = 3)
  while(i<=(length(DNAbuscar)-length+1)){
    A=(DNAbuscar[i:(i+length-1)])
    RepeticionesEnGenoma <- countPattern(A,DNAencontrar[[CHR]], 
                          max.mismatch=0,min.mismatch=0, fixed=TRUE)  
      print(RepeticionesEnGenoma)
    if(RepeticionesEnGenoma<=Repeats & RepeticionesEnGenoma > 0){
      if ((countPattern(reverseComplement(A),DNAencontrar[[CHR]],max.mismatch=0,min.mismatch=0, fixed=TRUE))==(Repeats-RepeticionesEnGenoma)){
        AB=toString(A)
        if(FW){
          AC=c(AB,(length(DNAbuscar)-i),Tmf2(DNAString(AB)))
        }else{
          AC=c(AB,i-1,Tmf2(DNAString(AB)))
        }
        print(AC)
        PosiblesPrimers=rbind(AC,PosiblesPrimers)
      }
    } 
    j=1
    i=i+1
  }
  colnames(PosiblesPrimers) <- c("Primer", "Location","Tm")
  return(PosiblesPrimers)
}

#Funci??n para ver cuanto se repiten ciertos primers en el genoma
SubProductos<- function(PrimersEnGenFW,PrimersEnGenBW,Genoma){
  MatProductos=matrix(, nrow = 0, ncol = 5)
  f=1
  while(f<=2){
    h=1
    g=1
    while(h<=nrow(PrimersEnGenFW)){
      while(g<=nrow(PrimersEnGenBW)){
        if(h==1){
          XXX=vmatchPattern(PrimersEnGenFW[h,1],Genoma, max.mismatch=0,min.mismatch=0, fixed=TRUE)
          YYY=vmatchPattern(PrimersEnGenBW[g,1],Genoma, max.mismatch=0,min.mismatch=0, fixed=TRUE)
        }else{
          XXX=vmatchPattern(toString(DNAString(PrimersEnGenFW[h,1])),Genoma, max.mismatch=0,min.mismatch=0, fixed=TRUE)
          YYY=vmatchPattern(toString(DNAString(PrimersEnGenBW[g,1])),Genoma, max.mismatch=0,min.mismatch=0, fixed=TRUE)
        }
        i=1
        o=0
        ListPrimersEnGenFW=list()
        while (i<=length(XXX)){
          if(!is.null(XXX[i]@ends[[1]])){
            o=o+1
            ListPrimersEnGenFW[[o]]=c(i, XXX[i]@ends[[1]])
          }
          i=i+1
        }
        
        i=1
        n=0
        ListPrimersEnGenBW=list()
        while (i<=length(YYY)){
          if(!is.null(YYY[i]@ends[[1]])){
            n=n+1
            ListPrimersEnGenBW[[n]]=c(i, YYY[i]@ends[[1]])
          }
          i=i+1
        }
        
        i=1
        j=1
        while (i<=length(ListPrimersEnGenFW)){
          while (j<=length(ListPrimersEnGenBW)){
            if (ListPrimersEnGenFW[[i]][1]==ListPrimersEnGenBW[[j]][1]){
              k=2
              l=2
              while (k<=length(ListPrimersEnGenFW[[i]])){
                while(l<=length(ListPrimersEnGenBW[[j]])){
                  print(paste(h,g,ListPrimersEnGenFW[[i]][1],(-ListPrimersEnGenFW[[i]][k]+ListPrimersEnGenBW[[j]][l])))
                  if(f==1 & (-ListPrimersEnGenFW[[i]][k]+ListPrimersEnGenBW[[j]][l])>0){
                    MatProductos=rbind(c(h,g,ListPrimersEnGenFW[[i]][1],TRUE,(-ListPrimersEnGenFW[[i]][k]+ListPrimersEnGenBW[[j]][l])),MatProductos)
                  }else if(f==2 & (ListPrimersEnGenFW[[i]][k]-ListPrimersEnGenBW[[j]][l])>0){
                    MatProductos=rbind(c(h,g,ListPrimersEnGenFW[[i]][1],FALSE,(ListPrimersEnGenFW[[i]][k]-ListPrimersEnGenBW[[j]][l])),MatProductos)
                  }
                  l=l+1
                }
                l=2
                k=k+1
              }
            }
            j=j+1
          }
          j=1
          i=i+1
        }
        g=g+1
      }
      g=1
      h=h+1
    }
    f=f+1
  }
  MatProductos=MatProductos[order(MatProductos[,1],MatProductos[,2]),]
  colnames(MatProductos) <- c("FW Primer","BW Primer", "CHR","Strand","Product Length")
  return(MatProductos)
}


#Funci??n para Gener matriz de primers (parte del Gen + Pl??smido)
UnionDePlasmidos <- function(ParteDePlasmidoEn5Prima,ParteDePlasmidoEn3Prima,FW){
  PosiblesPlasmidosCompletos=matrix(, nrow = 0, ncol = 5)
  colnames(PosiblesPlasmidosCompletos) <- c("Primer", "Location in Genome", "Tm Genome" ,"Location in Plasmid", "Tm Plasmid")
  i= nrow(ParteDePlasmidoEn5Prima)
  j=nrow(ParteDePlasmidoEn3Prima)
  while(j>0){
    while(i>0){
      ComparandoPlasmido=ParteDePlasmidoEn5Prima[i,1]
      ComparandoGen=ParteDePlasmidoEn3Prima[j,1]
      k=min(nchar(ComparandoPlasmido),nchar(ComparandoGen))
      while(k>=0){
        if(substr(ComparandoPlasmido, start=(nchar(ComparandoPlasmido)-k+1), stop=nchar(ComparandoPlasmido))==substr(ComparandoGen, start=1, stop=k)){
          BA=paste(ComparandoPlasmido, substr(ComparandoGen, start=k+1, stop=nchar(ComparandoGen)),sep="")
          BB=c(BA, ParteDePlasmidoEn3Prima[j,2],ParteDePlasmidoEn3Prima[j,3],ParteDePlasmidoEn5Prima[i,2],ParteDePlasmidoEn5Prima[i,3])
          PosiblesPlasmidosCompletos=rbind(BB,PosiblesPlasmidosCompletos)
        }
        k=k-1
      }
      i=i-1
    }
    i= nrow(ParteDePlasmidoEn5Prima)
    j=j-1
  }
  
  TmVector=c()
  HairpinVector=c()
  fila=1
  while(fila <= nrow(PosiblesPlasmidosCompletos)){
    PrimerActual=DNAString(PosiblesPlasmidosCompletos[fila,1])
    TmVector= c(TmVector, Tmf(PrimerActual))
    HairpinVector= c(HairpinVector,Hairpins(PrimerActual))
    fila=fila+1
  }
  PosiblesPlasmidosCompletos=cbind(TmVector,PosiblesPlasmidosCompletos)
  PosiblesPlasmidosVector=PosiblesPlasmidosCompletos[!HairpinVector]
  PosiblesPlasmidosCompletos=matrix(PosiblesPlasmidosVector,,6)
  if(FW){
    colnames(PosiblesPlasmidosCompletos) <- c("Tm", "Primer", "Location from SNP","Tm Genome", "Location in Plasmid","Tm Plasmid")
  }else{
    colnames(PosiblesPlasmidosCompletos) <- c("Tm", "Primer", "Location in Plasmid", "Tm Plasmid", "Location from SNP", "Tm Genome")
  }
  return(PosiblesPlasmidosCompletos)
}

# Funci??n para hacer las mejores parejas de FW y BW primers
Parejas = function (PrimersCompletosFW,PrimersCompletosBW,LongitudProductoMin,LongitudProductoMax,DiferenciaEnTm,LengthPrimerEnGen){
  i=1
  j=1
  ParejasPrimers=matrix(, nrow = 0, ncol = 8)
   while(i <= nrow(PrimersCompletosFW)){
    while(j <= nrow(PrimersCompletosBW)){
      if (abs((as.numeric(PrimersCompletosFW[i,"Tm Plasmid"]))-as.numeric(PrimersCompletosBW[j,"Tm Plasmid"]))<=DiferenciaEnTm 
          & abs((as.numeric(PrimersCompletosFW[i,"Tm Genome"]))-as.numeric(PrimersCompletosBW[j,"Tm Genome"]))<=DiferenciaEnTm
          & (as.numeric(PrimersCompletosFW[i,"Location from SNP"])+as.numeric(PrimersCompletosBW[j,"Location from SNP"])+LengthPrimerEnGen)>=LongitudProductoMin
          & (as.numeric(PrimersCompletosFW[i,"Location from SNP"])+as.numeric(PrimersCompletosBW[j,"Location from SNP"])+LengthPrimerEnGen)<=LongitudProductoMax
          & !Dimers(PrimersCompletosFW[i,"Primer"],PrimersCompletosBW[j,"Primer"])){
        DA=c(PrimersCompletosFW[i,"Primer"],PrimersCompletosFW[i,"Tm Plasmid"],PrimersCompletosFW[i,"Tm Genome"],
             toString(reverseComplement(DNAString(PrimersCompletosBW[j,"Primer"]))),
             PrimersCompletosBW[j,"Tm Plasmid"],PrimersCompletosBW[j,"Tm Genome"],
             as.numeric(PrimersCompletosFW[i,"Location from SNP"])+as.numeric(PrimersCompletosBW[j,"Location from SNP"])+LengthPrimerEnGen,
             as.numeric(PrimersCompletosFW[i,"Location from SNP"]))
        ParejasPrimers=rbind(DA,ParejasPrimers) 
        print(nrow(ParejasPrimers))
      }
      j=j+1
    }
    j=1
    i=i+1
   }
  colnames(ParejasPrimers) <- c("FW Primer", "FW Tm Plasmid", "FW Tm Genome" , "BW Primer", "BW Tm Plasmid", "BW Tm Genome","Product Length", "Bp Downstream From Primer")  
  return(ParejasPrimers)
}


#Secuencias:
Genoma = readDNAStringSet("GCF_000001405.36_GRCh38.p10_genomic.fna")    #Archivo Fasta Genoma
CHR = 1
Plasmido = readDNAStringSet("pEGFP.fasta") #Archivo Fasta Pl??smido
EnzimaFW = "AGATCT"                #Secuencia Enzima Forward (5'-3')
EnzimaBW = "GGATCC"                #Secuencia Enzima Reverse (5'-3')
abuscarA <- "TAAGAGCAGATCCCTGGACAGGCG"  #Secuencia para encontrar SNP

#Datos de la parte del primer en el genoma
Repeats=1     #Repeticiones en cromosoma
LengthPrimerEnGen=13  #Bases a hibridar en genoma
BasesUpstreamMin=20  #M??nimo de bases antes de SNP
BasesUpstreamMax=200 #M??ximo de bases antes de SNP
BasesDownstreamMin=20  #M??nimo de bases despu??s de SNP
BasesDownstreamMax=400  #M??ximo de bases despu??s de SNP
LongitudProductoMin=390  #L??ngitud de producto m??nima 
LongitudProductoMax=400 #L??ngitud de producto m??xima
LongitudSubproducto=5000 #l??ngitud de subproductos m??xima 
DiferenciaEnTm=4 #Diferencia de Tm para FW y RV Primers
TmPlasmidMin=78 #Incluye la suma de ambos primers que se hibridan en el pl??msido

#Datos de la parte del primer en el vector
LengthPrimerEnPlasmido=20

#Esto ya no se necesita introducir
DNAbuscarA = DNAString(abuscarA)

Localizacion <- matchPattern(abuscarA,
                          Genoma[[CHR]], max.mismatch=0,
                          min.mismatch=0, fixed=TRUE)
STRAND=TRUE
if(is.na(Localizacion@ranges@width[1])){
  Genoma = reverseComplement(Genoma)
  Localizacion <- matchPattern(abuscarA,
                            Genoma[[CHR]], max.mismatch=0,
                            min.mismatch=0, fixed=TRUE)
  STRAND=FALSE
  if(is.na(Localizacion@ranges@width[1])){
    print("Sequence not found")
    stop()
  }
}
SNPlocation=as.double(Localizacion@ranges@start)+length(DNAbuscarA)-1
SNPlocation
DNAbuscarFW = Genoma[[CHR]][(SNPlocation-BasesUpstreamMax+1):(SNPlocation-BasesUpstreamMin)]
DNAbuscarBW = Genoma[[CHR]][(SNPlocation+BasesDownstreamMin):(SNPlocation+BasesDownstreamMax-1)]

#A continuaci??n se corren las funciones
PrimersEnPlasmidoFW=PosiblesPrimersEnPlasmido(LengthPrimerEnPlasmido,EnzimaFW,Plasmido)
PrimersEnPlasmidoBW=PosiblesPrimersEnPlasmido(LengthPrimerEnPlasmido,EnzimaBW,Plasmido)

PrimersEnGenFW=PosiblesPrimersEnGen(Repeats,LengthPrimerEnGen,DNAbuscarFW,Genoma,CHR, TRUE)
PrimersEnGenBW=PosiblesPrimersEnGen(Repeats,LengthPrimerEnGen,DNAbuscarBW,Genoma,CHR, FALSE)

ProductosNoDeseados= SubProductos(PrimersEnGenFW,PrimersEnGenBW, Genoma)

PrimersCompletosFW=UnionDePlasmidos(PrimersEnPlasmidoFW,PrimersEnGenFW, TRUE)
PrimersCompletosBW=UnionDePlasmidos(PrimersEnGenBW,PrimersEnPlasmidoBW, FALSE)

ParejasPrimers=Parejas(PrimersCompletosFW,PrimersCompletosBW,LongitudProductoMin,LongitudProductoMax,DiferenciaEnTm,LengthPrimerEnGen)

write.csv(ParejasPrimers, file=paste("PrimersParaClonaci??n", Sys.time(),".csv", sep=""))


PrimErendira = DNAString("agccttccttcctgggcatg")
PrimErendira2=DNAString("gagcaatgatcttgatcttc")
BuscaPrimerEndira <- vcountPattern(PrimErendira,
                                   Genoma, max.mismatch=0,
                                   min.mismatch=0, fixed=TRUE)
BuscaPrimerEndira2 <- vcountPattern(PrimErendira2,
                                   Genoma, max.mismatch=0,
                                   min.mismatch=0, fixed=TRUE)
BuscaPrimerEndira3 <- vcountPattern(reverseComplement(PrimErendira),
                                   Genoma, max.mismatch=0,
                                   min.mismatch=0, fixed=TRUE)
BuscaPrimerEndira4 <- vcountPattern(reverseComplement(PrimErendira2),
                                    Genoma, max.mismatch=0,
                                    min.mismatch=0, fixed=TRUE)
PrimEreToF = BuscaPrimerEndira > 0 & BuscaPrimerEndira4 > 0
PrimEre2ToF = BuscaPrimerEndira2 > 0 & BuscaPrimerEndira3 > 0

i=1
Cromosomas= c()
while (i<= length(PrimEreToF)){
  if(PrimEreToF[i]){
    Cromosomas = c(Cromosomas,i)
  }
  i=i+1
}

i=1
Cromosomas2= c()
while (i<= length(PrimEreToF)){
  if(PrimEre2ToF[i]){
    Cromosomas2 = c(Cromosomas2,i)
  }
  i=i+1
}

names(Genoma[Cromosomas])
names(Genoma[Cromosomas2])

Tama??os = c()
for (i in 1:length(Cromosomas)){
  EncuentraPrimerEndira <- matchPattern(PrimErendira,
                                        Genoma[[Cromosomas[i]]], max.mismatch=0,
                                        min.mismatch=0, fixed=TRUE)
  EncuentraPrimerEndira2 <- matchPattern(reverseComplement(PrimerEndira2),
                                        Genoma[[Cromosomas[i]]], max.mismatch=0,
                                        min.mismatch=0, fixed=TRUE)
  Longi=EncuentraPrimerEndira2@ranges@start-EncuentraPrimerEndira@ranges@start+EncuentraPrimerEndira[1]@ranges@width
  print(names(Genoma[Cromosomas[i]]))
  print(Longi)
}

for (i in 1:length(Cromosomas2)){
  EncuentraPrimerEndira <- matchPattern(PrimErendira2,
                                        Genoma[[Cromosomas2[i]]], max.mismatch=0,
                                        min.mismatch=0, fixed=TRUE)
  EncuentraPrimerEndira2 <- matchPattern(reverseComplement(PrimErendira),
                                         Genoma[[Cromosomas2[i]]], max.mismatch=0,
                                         min.mismatch=0, fixed=TRUE)
  #print(EncuentraPrimerEndira)
  #print(EncuentraPrimerEndira2)
  Longi=EncuentraPrimerEndira2@ranges@start-EncuentraPrimerEndira@ranges@start+EncuentraPrimerEndira[1]@ranges@width
  print(names(Genoma[Cromosomas2[i]]))
  print(Longi)
}


names(Genoma[31])
EncuentraPrimerEndira[64]