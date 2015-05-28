#!/usr/bin/env Rscript

## Load libraries
suppressMessages(library(ggplot2))
suppressMessages(library(optparse))

## Create options
option_list <- list(
    make_option(c("-W", "--width"),type="double", default=17,
                help="Width of the output file, [default %default]"),
    make_option(c("-H", "--height"), type="double",default=10,
                help="Height of the output file, [default %default]"),
    make_option(c("-f", "--fill"), action="store_true", default=FALSE,
                help="Fill density plot, [default %default]")
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

## Generate plot object
p <- ggplot(data=toparun, aes(x=(Upper-Lower)/2, color = Name, group=Name)) +
    geom_density(alpha=.2) + theme_bw() + xlab("HUW")

## Modify if fill=TRUE
if (opt$fill) {p <- p + aes(fill = Name)}

## Plot
name = sub("[.].{2,3}","",arguments$args)
ggsave(p, file= paste(name, "_plot.pdf", sep = ""), width=opt$width, height=opt$height, units="cm")
ggsave(p, file= paste(name, "_plot.png", sep = ""), width=opt$width, height=opt$height, units="cm")
