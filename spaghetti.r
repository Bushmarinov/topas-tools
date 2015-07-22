#!/usr/bin/Rscript
options(warn=-1)
#Sys.setenv(LANG="EN")
suppressMessages(library(data.table))
suppressMessages(library(gdata))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(optparse))

k.table <- c(5.62619251369954 ,3.59689767105988, 2.99683950381941, 2.9141936236358 , 
3.37232041437195 ,2.93196313616901, 2.72771585462825, 2.70417682172506,
2.86576896853064 ,2.66921102377684, 2.57092705221669, 2.55881443841962,
2.64995441497466 ,2.53075618782533, 2.47791909362932, 2.4741828679544 ,
2.53384546494551 ,2.45983379319381, 2.4212286638301 , 2.4159429171073 ,
2.45784825894735 ,2.40879237590095, 2.38090443345185, 2.38216835038994,
2.41048201093783 ,2.36989342570794, 2.35215022384351, 2.3524695182014 ,
2.37772035649117 ,2.34643406722578, 2.33269558844597, 2.33379314746303,
2.35213581186961 ,2.3293126888992 , 2.31585973838799, 2.32211834728897,
2.33719246258054 ,2.31822219208683, 2.30836925652429, 2.30928647548234,
2.32053355883361 ,2.30880999917026, 2.29983955811199, 2.30081635190135,
2.31427557030017 ,2.29827012906655)



option_list <- list(
    make_option(c("-W", "--width"),type="double", default=15,
                help="Width of the output file, [default %default]"),
    make_option(c("-H", "--height"), type="double",default=10,
                help="Height of the output file, [default %default]"),
    make_option(c("-p", "--Plot"), action="store_true", default=FALSE,
                help="Use existing _rtable.dat for plotting, [default %default]"),
    make_option(c("-e", "--Errors"), action="store_true", default=FALSE,
                help="Use bond errors in outlier detection, [default %default]"),
    make_option(c("-k", "--Kpar"), type="double",default=-1,
                help="Use k in the outlier detection function, negatives for table values, [default %default]"),
    make_option(c("-a", "--Angles"), action="store_true", default=FALSE,
                help="Use angles instead of bonds")			
    )

parser    <- OptionParser(usage = "spaghetti.r [options] directory", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt       <- arguments$options

if (length(arguments$args) != 1)
{
    print("Incorrect number of required positional arguments")
    print_help(parser)
    stop()
}


directory <-   arguments$args


if (opt$Kpar < 0){
    kpar <- function(l) {if (l < 51) k.table[l-4] else k.table[46]}
}else{
    kpar <- function(l){opt$Kpar}
}



if (opt$Angles) {pattern <- "^\\s*Angle_Restrain(_Breakable)?"
                 name.su <- "_angles"
             }else{
                 pattern = "^\\s*Distance_Restrain(_Breakable|_Morse)?"
                 name.su <- "_bonds"}

na.to.0 <- function(x) {if (is.na(x)) 0 else x}

extract.data <- function (filename) ###Индия!
    {
        xs <- readLines(filename)
        xs <- gsub("^\\s*'.*","",xs, perl=TRUE) #kill full-line comments
        k1 <- unlist(strsplit(grep("penalties_weighting_K1",xs, value=TRUE), "[[:space:]]+"))
        k1 <- as.numeric(k1[length(k1)])
		rwp <- gsub("^.*r_wp\\s+([0-9.]+).*$", "\\1", grep("r_wp", xs, perl=TRUE, value=TRUE), perl=TRUE)[1]
		rwp <- as.numeric(rwp)
		# print(rwp)
        ys <- grep(pattern,xs, value=TRUE, perl=TRUE)
        ys <- ys[!grepl("H\\d",ys, perl=TRUE)]
        ys <- strsplit(ys, "([ \t]*,[][ \t]*)|[()][ \t]*")
        l.ys  <- length(ys)
        Delta <- numeric(l.ys)
        Bond  <- numeric(l.ys)
        Error <- numeric(l.ys)
        for (i in 1:l.ys){
            y <- ys[[i]]
            num <- unlist(strsplit(y[4], "(`_|`)"))
            Delta[i] <- as.numeric(num[1]) - as.numeric(y[3])
            Error[i] <- na.to.0(as.numeric(num[2]))
            Bond[i]  <- y[2]
        }
        data.frame(K1=k1, Rwp=rwp, Bond=Bond, Delta=Delta, Error=Error)
    }

calc.outliers <- function(table, k.fun = kpar, use.errors.flag = opt$Errors)
    {
        l       <- length(table$Delta)
        qq      <- quantile(table$Delta, c(0.25, 0.75), names=FALSE)
        table[["Av.Error"]] <- mean(table$Error)
        table[["Low.Q"]]    <- qq[1] - k.fun(l) * diff(qq)
        table[["High.Q"]]   <- qq[2] + k.fun(l) * diff(qq)
        table[["Outlier"]]  <- (table$Delta > (table$High.Q + 2 * table$Error * use.errors.flag)) |
                               (table$Delta < (table$Low.Q  - 2 * table$Error * use.errors.flag))
        table
    }

cat("Reading data..\n")

if (opt$Plot){
    errors.table <- read.table(directory, sep="\t", header=TRUE)
    directory <- sub("[.].{2,3}","",directory)
}else{
    topas.outs <- list.files(directory, full.names=TRUE, pattern="*.(out|OUT)$")
    errors.table <- rbindlist(Map(extract.data, topas.outs))}

cat("Some calculations...\n")
# print(errors.table)
errors.table <- rbindlist(by(errors.table, errors.table$K1,   calc.outliers))
errors.table <- rbindlist(by(errors.table, errors.table$Bond, function(x) {x[["Bad.Bond"]] <-  any(x[["Outlier"]]); return(x)}))

errors.summary <- rbindlist(by(errors.table, errors.table$K1, 
	function(x){ 
		data.frame(	
			K1=x[["K1"]][1], 
			Rwp=x[["Rwp"]][1],
			RMS=round(sqrt(sum((x[["Delta"]])^2/length(x[["Delta"]]))), digits=4),
			Outliers=sum(x[["Outlier"]])
		)
	}
))
counts <- by(errors.table, errors.table$K1, nrow)
cat("IQR multiplier(s) used", sapply(unique(counts), function(x) {sprintf("%.4f", kpar(x))}),"\n")

## Plotting settings
 
p <- ggplot(data=errors.table, aes(x=K1, y=Delta)) + theme_bw()

if ( opt$Errors &  sum(errors.table$Av.Error) >  0) {
	p <- p +
		geom_ribbon(aes(ymin=High.Q-2*Av.Error, ymax=High.Q+2*Av.Error, y=NULL), fill="#999999") +
		geom_ribbon(aes(ymin=Low.Q -2*Av.Error, ymax=Low.Q +2*Av.Error, y=NULL), fill="#999999") 
	}else {
            p <- p +  geom_line(aes(y = High.Q), linetype="dashed") + geom_line(aes(y = Low.Q), linetype="dashed") }
	
 if (sum(errors.table$Outlier) > 0){
     p <- p +
         geom_line(data=errors.table[errors.table$Bad.Bond==FALSE,], aes(group=Bond), color="black") +
             geom_line(data=errors.table[errors.table$Bad.Bond==TRUE,], aes(group=Bond, color = Bond)) +
                 geom_point(data=errors.table[errors.table$Outlier==TRUE,], aes(group=Bond, color = Bond), size=2) +
				 scale_colour_discrete(name="Outlying bonds")
 } else {
	p <- p + geom_line(aes(group=Bond) , colour="black")  }

## Output

cat("Plotting and tables writing...\n")

write.table(format(errors.table[order(errors.table$K1,errors.table$Bond),], digits=6),
            paste(directory, name.su, "_rtable.dat", sep=""), row.names=FALSE, sep="\t", quote = FALSE)	
write.fwf(errors.summary[order(errors.summary$K1)], paste(directory, name.su,"_summary.txt", sep=""))

ggsave(p, file= paste(directory, name.su, "_plot.pdf", sep = ""), width=opt$width, height=opt$height, units="cm")
ggsave(p, file= paste(directory, name.su, "_plot.eps", sep = ""), width=opt$width, height=opt$height, units="cm")

