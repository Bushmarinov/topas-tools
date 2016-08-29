suppressPackageStartupMessages(library("optparse"))

# Functions:
	## prepare.for.plot (multiplies data on 'mul' after 'beg' ):
		prepare.for.plot <- function (good.line, scale, begining, multiplicator)
			{
			good.line[scale>begining] <- good.line[scale>begining]*multiplicator
			return(good.line)
			}

	## hkl.plot (plots hkl-lines):
		hkl.plot <- function(hkl.vector, y.geometry)
	           	  {mapply(function(x)
	                      {lines(c(x,x), y.geometry)
	                       return(NULL)},
	                     hkl.vector)
	              return(NULL)}

 # Command-line options:

option_list <- list(
                    make_option(c("-o", "--output"), default="topasplot.eps",
                                help="Name of output file, [default \"%default\"]"),
                    make_option(c("-W", "--width"), type="integer", default=19,
                                help="Width of the PS plot [default %default]"),
                    make_option(c("-H", "--height"), type="integer", default=10,
                                help="Height of the PS plot  [default %default]"),
                    make_option(c("-b", "--beg"), type="double", default=35,
                                help="Starting value for multiplying, [default %default]. Let it be negative to disable multiplying"),
                    make_option(c("-m", "--mul"), type="double", default=5,
                                help="Multiplier [default %default]"),
                    make_option(c("-c", "--corr"), type="double", default=3,
                                help="Can help with some margins problems, [default %default]"),
                    make_option(c("-l", "--log"), action="store_true", default=FALSE,
                                help="Plot graph in logarithmic coordinates"),
                     make_option(c("-d", "--DisableDiff"), action="store_true", default=FALSE,
                                help="Do not plot the difference curve")
                    )

parser    <- OptionParser(usage = "%prog [options] plot.data.file [peak.data.file]", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt       <- arguments$options

if (length(arguments$args) ==1) {isPeaks <- FALSE} else {isPeaks <- TRUE}


# Reading data:

topas.data <- read.table(arguments$args[1])
if (isPeaks == TRUE) {
peak.data  <- read.table(arguments$args[2])
}

# Modifying data with log  or identity and prepare.for.plot.

if ( opt$log) {
    modify <-   function (x) {return(log(x))}
    mtext.text <- "log(I)"
} else {
    modify <-   function (x) {return(x)}
    mtext.text <-  "Intensity, counts"
}

if (opt$beg < 0){
    topas.data[[2]] <- modify(topas.data[[2]])
    topas.data[[3]] <- modify(topas.data[[3]])
    topas.data[[4]] <- topas.data[[2]] - topas.data[[3]]
} else {
    topas.data[[2]] <- prepare.for.plot(modify(topas.data[[2]]), topas.data[[1]], opt$beg, opt$mul)
    topas.data[[3]] <- prepare.for.plot(modify(topas.data[[3]]), topas.data[[1]], opt$beg, opt$mul)
    topas.data[[4]] <- topas.data[[2]] - topas.data[[3]]
}

# Graphics parameters:

plot.width  <-  opt$width  / 2.54
plot.height <-  opt$height / 2.54
psfile      <-  opt$output


setEPS()
postscript(psfile, horizontal=FALSE, width=plot.width, height=plot.height, paper="a4")
par(las=1, mar=c(3, 5, 1, 1), lwd = 0.3, mgp = c(40, 1, 0))
plot.new()


x.min <- min(topas.data[[1]]) + opt$corr
x.max <- max(topas.data[[1]]) - opt$corr
if (!opt$Disable){
    y.min <- min(topas.data[-1])
    y.max <- max(topas.data[-1])
    legend.legend <- c("Experimental data", "Calculated data", "Difference")
    legend.col <- c("black", "red", "grey30")
} else {
    y.min <- min(topas.data[-c(1, 4)])
    y.max <- max(topas.data[-c(1, 4)])
    legend.legend <- c("Experimental data", "Calculated data")
    legend.col <- c("black", "red")
}

if (isPeaks == TRUE) {
	hkl.y <- (y.max -y.min) * 0.02
	y.min <- y.min - 1.3 * hkl.y
	}

# Plot!
options(scipen=5)
plot.window(c(x.min, x.max), c(y.min, y.max))
axis(1)
axis(2)
box(lwd=1)
mtext("2Theta, degrees", side=1, line=2)
mtext(mtext.text, side=2, line=4, las=0)
legend("topright", legend = legend.legend, lty=1, col = legend.col, merge = TRUE, cex=0.7)

lines(topas.data[[1]], topas.data[[2]], col="black")
lines(topas.data[[1]], topas.data[[3]], col="red")

if ( !opt$DisableDiff){
	lines(topas.data[[1]], topas.data[[4]], col="grey30")}

if (opt$beg >= 0) {
    abline(v = opt$beg, lty="dashed")
    text(x= opt$beg, y = y.max, pos = 4, paste("x", as.character(opt$mul)), cex=0.7)}

if (isPeaks == TRUE) {
    hkl.plot(peak.data[[1]], c(y.min, y.min + hkl.y))
}

dev.off()
