mydata= read.csv(file=file.choose())#Selecciona un archivo .csv 
y=mydata[["x"]] #Cambia Altura por el t??tulo dela columna que qeieras graficar en y
x=mydata[["y"]] #Cambia Peso por el t??tulo dela columna que qeieras graficar en x
fit=lm(y ~ x) # makes linear model
plot(x,y) # plots y vs x
abline(fit) # puts regresion in plot
#plot(fit) # plots the lineae model
