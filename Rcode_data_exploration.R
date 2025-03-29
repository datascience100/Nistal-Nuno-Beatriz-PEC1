

library(Biobase)


library(readxl)
datafile = read_excel("GastricCancer_NMR (final).xlsx", sheet=1) # leo la Hoja de Datos
datafile_matrix= as.matrix(datafile)  # convierto la hoja de excel a matriz
rownames(datafile_matrix)=datafile_matrix[ , 2]  #asigno a los nombres de las filas la columna del número de muestra
datafile_matrix_reduced=datafile_matrix[ ,c(-1,-2,-3,-4)]  # elimino las columnas que contienen los metadatos de cada muestra

# Transposing the matrix
t_datafile_matrix_reduced <- t(datafile_matrix_reduced)
dim(t_datafile_matrix_reduced)
numbers=as.numeric(t_datafile_matrix_reduced[ ,1:140])  # convierto la matriz de concentraciones de metabolitos a números analizables
t_datafile_matrix_reduced= matrix(numbers, nrow=149, byrow=FALSE)

rownames(t_datafile_matrix_reduced)=colnames(datafile_matrix_reduced)  # asigno nombres a las filas
colnames(t_datafile_matrix_reduced)=rownames(datafile_matrix_reduced)  # asigno nombres a las columnas
t_datafile_matrix_reduced[1:5,1:5]


sample_metadata = datafile_matrix[ ,c(1,2,3,4)]  # selecciono los metadatos de las muestras que luego formarán parte del objeto SummarizedExperiment 
rownames(sample_metadata)=sample_metadata[ , 1] # asigno nombres a las filas
sample_metadata= sample_metadata[,-1]  # elimino la primera columna
sample_metadata[1:5,]
dim(sample_metadata)


metabolite_metadata <- read_excel("GastricCancer_NMR (final).xlsx", sheet=2 )  # leo la Hoja de metabolitos
metabolite_metadata=as.matrix(metabolite_metadata)  # convierto la hoja de excel a matriz

rownames(metabolite_metadata)=metabolite_metadata[,1]
metabolite_metadata=metabolite_metadata[,-1]  # elimino la primera columna


sample_metadata= as.data.frame(sample_metadata) # convierto la matriz de metadatos a dataframe, que es el tipo de datos que necesita el SummarizedExperiment
metabolite_metadata= as.data.frame(metabolite_metadata) # convierto la matriz de metabolitos a dataframe, que es el tipo de datos que necesita el SummarizedExperiment



sapply(metabolite_metadata, class)  # compruebo los tipos de datos de las columnas de los metabolitos
metabolite_metadata$Perc_missing=as.numeric(metabolite_metadata$Perc_missing)  # convierto esta columna a numérica
metabolite_metadata$QC_RSD=as.numeric(metabolite_metadata$QC_RSD)  # convierto esta columna a numérica


library("SummarizedExperiment")


stopifnot(rownames(t_datafile_matrix_reduced) == metabolite_metadata$Name) # compruebo que los nombres coinciden
stopifnot(colnames(t_datafile_matrix_reduced) == sample_metadata$SampleID) # compruebo que los nombres coinciden
rownames(sample_metadata)=sample_metadata$SampleID
rownames(metabolite_metadata)=metabolite_metadata$Name
se <- SummarizedExperiment(assays = list(counts = t_datafile_matrix_reduced), colData = sample_metadata, rowData = metabolite_metadata) # construyo el SummarizedExperiment
se  # visualizo el SummarizedExperiment
assay(se)[1:5,1:5]   # visualizo la matriz de concentraciones que forma parte del SummarizedExperiment


se1 <- se[rowData(se)$Perc_missing==0 & rowData(se)$QC_RSD<20,] # elimino los metabolitos según los criterios descritos anteriormente
se1   # visualizo el nuevo SummarizedExperiment

assay(se1)[1:5, 1:8]


apply(assay(se1),2, summary)[1:5,1:8]  # Column-wise summary statistics

opt <- par(mfrow=c(1,3))
hist(assay(se1)[,1]) # visualizo los histogramas
hist(assay(se1)[,2])
hist(assay(se1)[,3])
par(opt)


colores <- c()  # asigno a cada muestra un color distinto para la visualización en los gráficos según el tipo de muestra 
for (i in 1: 140) {
  if (colData(se1)$Class[i] =="QC"){colores<-c(colores,"blue")}
  if (colData(se1)$Class[i] =="GC"){colores<-c(colores,"red")}
  if (colData(se1)$Class[i] =="BN"){colores<-c(colores,"pink")}
  if (colData(se1)$Class[i] =="HE"){colores<-c(colores,"green")}
  
}  
boxplot(assay(se1), col=colores, xlab="samples", ylab="metabolite concentrations ", las=2, cex.axis=0.7, cex.main=0.7)
rownames(colData(se1))=colData(se1)$Class  # cambio los nombres de las filas a los tipos de muestra para diferenciarlas según el tipo de muestra


logX <- log10(assay(se1))  # Log scale (base-10)

boxplot(logX, col=colores, xlab="samples", ylab="metabolite concentrations ", las=2, cex.axis=0.7, cex.main=0.7)



pcX<-prcomp(t(logX), center = TRUE, scale = TRUE) #  escalando y centrando los datos
loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)

  
  
  
xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))

plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colores, main ="Principal components (PCA)")
legend("bottomleft", 
       legend = c("GC", "BN", "HE", "QC"), 
       col = c("red", "pink", "green", "blue"),
       pch = c(1,1), 
       bty = "n", 
       pt.cex = 2, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))



# Realizamos agrupamiento jerárquico 

clust.euclid.average <- hclust(dist(t(logX)),method="average") # Calcula la matriz de distancias y realiza el agrupamiento
plot(clust.euclid.average, hang=-1, xlab = "samples") # Representamos el dendrograma	


manDist <- dist(t(logX))
heatmap (as.matrix(manDist), col=heat.colors(16))

