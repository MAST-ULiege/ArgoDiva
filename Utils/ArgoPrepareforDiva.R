# IN 
# pars: a vector of parameters chosen by the user for DIVAND
# initdf : initial data frame to process
#
# OUT
# a text file named ArgoforDiva.txt

ArgoPrepareforDiva<- function(pars, initdf){

df <- subset(initdf, select=pars)  

#remove rows containing NA's
df <- df[complete.cases(df),]

write.table(df, file="ArgoforDiva.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
                            row.names = FALSE, col.names = FALSE)    

}
