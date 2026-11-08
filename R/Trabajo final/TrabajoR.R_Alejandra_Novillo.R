##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R       ##################
#####################                                       			##################
##########################################################################################

#Primero miramos si estamos en el directorio de trabajo que queremos
getwd()


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t")
#Tenemos que leer el archivo y tener en cuenta que hay un tabulador separando las columnas.
#Examinanmos el archivo (en este caso se llama datos)
head(datos)     # primeras 6 filas
tail(datos)     # últimas 6 filas
dim(datos)      # número de filas (7) y columnas (619)


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. 

# La primera columna es el nombre, las 2-4 son los controles y de 5-7 las condiciones. 

datos$media_control <- apply(datos[,2:4],1,mean)    # Media control (2-4)
datos$de_control <- apply(datos[,2:4],1,sd)         # DE control (2-4)
datos$media_condicion <- apply(datos[,5:7],1,mean)  # Media condiciones (5-7)
datos$de_condicion <- apply(datos[,5:7],1,sd)       # DE condiciones (5-7)

#Comprobamos (cuidado al cambiar el nombre después de ejecutar porque duplica las columnas)
head(datos)
dim(datos)   # debería dar 619 filas 11 columnas

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

rango <- range(c(datos$media_control, datos$media_condicion)) #Primero ajustamos los rangos (max y mín) para que la escala del eje X sea la misma y se puedan comparar los histogramas

par(mfrow=c(1,2))                           #Ponemos para que salgan los dos histogramas al lado
hist(datos$media_control,                   #Histograma para la media del control
     main = "Control", col = "lightgreen",  #Título y color verde
     xlab = "Media", ylab = "Frecuencia")   #Título del eje (label)

hist(datos$media_condicion,                 #Histograma para la media de la condición
     main = "Condición", col = "pink",      # Título y color rosa
     xlab = "Media", ylab = "Frecuencia")   # Títulos de los ejes

par(mfrow=c(1,1))                           # Volvemos a poner solo un panel, no dos

#Podemos asociar los colores a las variables si queremos, ya que los vamos a usar más adelante
#col_control   <- "lightgreen"
#col_condicion <- "pink"

#Como los valores son muy heterogéneos, hacemos el logartimo para visualizar el histograma más claro.
#Le sumamos uno (+1) por si hay algún valor que da 0, para no perder esos datos al hacer log, ya que el log(0)= -infinito
#los valores dan como que no existen y los perdemos

datos$log_control   <- log2(datos$media_control + 1)
datos$log_condicion <- log2(datos$media_condicion + 1)

# Para comprobar que se han añadido las columnas, ahora debe haber 13
head(datos)
dim(datos)

# Volvemos ver histogramas con el cambio log:
rango <- range(c(datos$media_control, datos$media_condicion))
par(mfrow=c(1,2))                                    
hist(datos$log_control,                          # Histograma log2 media control + 1
     main = "Control", col = "lightgreen",       
     xlab = "log2 (Media)", ylab = "Frecuencia")   

hist(datos$log_condicion,                        # Histograma log2 media condición + 1
     main = "Condición", col = "pink",         
     xlab = "log2 (Media)", ylab = "Frecuencia")   
par(mfrow=c(1,1))                                # Volvemos a poner solo un panel, no dos

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)


boxplot(datos$log_control, datos$log_condicion,                           #boxplot para control y condición
        names = c("Control", "Condición"),                                #nombres de eje X
        col   = c("lightgreen", "pink"),                                  #colores con las variales
        main  = "Boxplot de las medias (log2)",ylab  = "log2(Media + 1)", #título y nombre de eje Y
        ylim  = c(0, 25))                                                 #el eje Y va de 0 a 25
legend("topleft", legend = c("Control", "Condición"),                     #leyenda arriba izq
       fill = c("lightgreen", "pink"))                                    #colores correspondientes

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

shapiro.test(datos$media_control)   #test normalidad media control
shapiro.test(datos$media_condicion) #test normalidad media condición
#nuestra H0 es que los valores no salen significativos y siguen una distribución normal
#nuestra H1 es que los valores salen significativos y no siguen una distribución normal

#como ambos son  p-value < 2.2e-16, H0 no es cierta,
#utilizamos un test no paramétrico para cada gen, como Wilcox

datos$p_value <- apply(datos[, 2:7], 1,                                   # p-value en "datos"
                       function(f) wilcox.test(f[1:3], f[4:6])$p.value)   # Wilcox test

#comprobamos que se ha creado la columna: 
head(datos)
dim(datos)

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

significativo_0.05 <- datos[datos$p_value <= 0.05, ]   # genes diferencialmente expresados con  corte de p-value 0.05
nrow(significativo_0.05)                               # cuantas filas hay en la nueva variable

#el test de wilcox utiliza el orden de las variables, no el valor en sí, por lo que mezcla los 3 de control y los 3 de condición
#luego los ordena de menor a mayor. Con 6 valores en dos grupos de 3, solo hay 20 formas posibles de repartirlos. 
#los dos extremos de esas 20 combinaciones es que salgan los dos grupos separados (controles y condiciones o viceversa)
#el p value es la probabilidad de que salgan esos dos extremos 2/20=0.1. Son las dos colas de nuestra curva.
#si ponemos el p value en 0.05 solo tiene en cuenta uno de los dos extremos, es decir, una de las colas
#podría salir algún valor significativo si ponemos el pvalue en 0.05 porque sería un extremo (una cola), pero aqui tenemos dos direcciones: sobreexpresados y reprimidos
#por lo que necesitamos las dos colas de la curva. Es por eso que con un p value de 0.05 no nos sale ningún valor significativo,

significativo_0.1 <- datos[datos$p_value <= 0.1, ]   # genes diferencialmente expresados con  corte de p-value 0.05
nrow(significativo_0.1)

#obtenemos 429 diferencialmente expresados para un p value de 0.1

#para saber cuántos han salido sobreexpresados y cuantos reprimidos:
sobreexpresados <- significativo_0.1[significativo_0.1$media_condicion > significativo_0.1$media_control, ]
reprimidos <- significativo_0.1[significativo_0.1$media_condicion < significativo_0.1$media_control, ]

nrow(sobreexpresados)    # cuántos sobreexpresados
nrow(reprimidos)         # cuántos reprimidos

# Comprobación: 429 total 
nrow(sobreexpresados) + nrow(reprimidos)
nrow(significativo_0.1)

#conclusión: de los 619 genes hay 429 diferencialmente expresados, de los cuales 69 eran sobreexpresados y 360 reprimidos

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

# matriz con los genes diferencialmente expresados:
filtrado <- as.matrix(significativo_0.1[, 2:7])
rownames(filtrado) <- significativo_0.1[, 1]

# creamos los dendrogramas haciendo un cluster de columnas y genes:
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtrado))))) #Traspone la tabla,y transforma las diferencias en la distancia de las ramas del dendograma
colv = as.dendrogram(hclust(as.dist(1-cor(filtrado))))    #Lo mismo pero sin trasponer porque ya son columnas

#hacemos el heatmap 1
heatmap(filtrado,
        Rowv = rowv,Colv = colv,                         #heatmap del dendograma y la matriz
        cexCol=0.7,                                      #ajustar el tamaño de letra del eje x
        labRow = FALSE,                                  #oculta los nombres de los genes porque 429 no caben
        col=hcl.colors(50, palette = "Blue-yellow 3"),   #paleta de 50 tonos de azul a amarillo
        margins = c(8, 2))                               #ajusta márgenes
     

#hacemos el heatmap 2 más elaborado. Para ello instalamos:´
install.packages("gplots")		        
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

heatmap(filtrado,
        Rowv = rowv,Colv = colv,                     #heatmap del dendograma y la matriz
        cexCol=0.7,                                  #ajustar el tamaño de letra del eje x
        labRow = FALSE,                              #oculta los nombres de los genes porque 429 no caben
        col = rev(redblue(256)),                     #colores (rojo sobreexpresados, azul reprimidos) de 256 tonos
        scale = "row",                               #estandariza los genes por su z score (valor-media/SD). Así el color sale de distinta intensidad,
                                                     #si no usaría los valores absolutos y todo saldría del mismo color
        trace = "none", density.info = "none",       #quita líneas cian y la leyenda que salen por defecto
               margins = c(8, 2))                    #ajusta márgenes

  

