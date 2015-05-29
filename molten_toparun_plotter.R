#!/usr/bin/env Rscript

## Load libraries
suppressMessages(library(ggplot2))
suppressMessages(library(optparse))

## Create options
option_list <- list(
    make_option(c("-W", "--width"),type="double", default=17,
                help="Width of the output image, cm [default %default]"),
    make_option(c("-H", "--height"), type="double",default=10,
                help="Height of the output image, cm [default %default]"),
	make_option(c("-R", "--rows"), type="integer",default=1,
                help="Number of rows in legend [default %default]"),
	make_option(c("-b", "--brewer"), type="character",default="Set1",
                help="ColorBrewer 2.0 pallette to use in brewer output [default %default]"),
    make_option(c("-f", "--fill"), action="store_true", default=FALSE,
                help="Fill density plot, [default %default]"),
	make_option(c("-l", "--linetype"), action="store_true", default=FALSE,
                help="Vary linetype as well as plot color, [default %default]"),
	make_option(c("-L", "--line"), action="store_true", default=FALSE,
                help="Plot lines, not areas [default %default]")
    )

parser    <- OptionParser(usage = "molten_toparun_plotter.r [options] data.txt\n were data.txt is the output of melt_toparun.pl", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt       <- arguments$options

larg <- length(arguments$args)
if (larg != 1)
{
    print_help(parser)
    stop(paste("Wrong number of arguments (", larg, " for 1)", sep=""))
}


## Read data
toparun <- read.table(arguments$args, sep="\t", header=T)
if (!is.factor(toparun$Name)) {
  toparun$Name <- factor(toparun$Name, levels = sort(unique(toparun$Name)), ordered=TRUE)   # convert integer dataset names to ordered factors
}

view <- geom_density(alpha=0.2)

if (opt$line) {
	view <- geom_line(stat="density")
} 

## Generate plot object
p <- ggplot(data=toparun, aes(x=(Upper-Lower)/2, colour=Name, group=Name)) + view +
    theme_bw() + xlab("HUW")+ theme(legend.position = "bottom") + labs(colour=NULL) 

## Modify if fill=TRUE
if (opt$fill) {p <- p + aes(fill = Name)}

if (opt$linetype) {p <- p + aes(linetype = Name) + labs(linetype=NULL)}



if (opt$rows > 1) {
	g <- guide_legend(nrow = opt$rows, byrow = TRUE)
	p <- p + guides(col = g, linetype = g)
}

## Plot
name = sub("[.].{2,3}","",arguments$args)
# ggsave(p, file= paste(name, "_plot.pdf", sep = ""), width=opt$width, height=opt$height, units="cm")
ggsave(p, file= paste(name, "_plot.eps", sep = ""), width=opt$width, height=opt$height, units="cm")
ggsave(p +scale_colour_brewer(palette=opt$brewer) , file= paste(name, "_plot_brewer.eps", sep = ""), width=opt$width, height=opt$height, units="cm")
# ggsave(p, file= paste(name, "_plot.png", sep = ""), width=opt$width, height=opt$height, units="cm")
